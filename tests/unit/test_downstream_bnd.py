"""Same-attempt generated BND handoff using synthetic children, not real callers.

Structural BGZF/TBI fixtures do not prove tabix readability or SV accuracy.
The full native launcher/tool gate is separate from these contract regressions.
"""
import base64
import copy
import csv
import gzip
import json
import os
import sys
from pathlib import Path

import pytest

import hifivar.attempt_receipts as receipts
from hifivar.bnd_records import derive_record_stream
from hifivar.breakends import BndPolicy, POLICY, POLICY_SHA256
from hifivar.command import CommandRunner
from hifivar.downstream_bnd import generated_bnd_roles
from hifivar.downstream_handoff import materialize_downstream
from hifivar.exceptions import InputValidationError, OutputValidationError
from hifivar.phase9 import phase9_bnd_paths
from hifivar.validation import compute_sha256
from test_downstream_handoff import handoff_fixture, _seal, _row, _bgzf, _vcf


def seal_bundle(frozen):
    canonical = frozen.root/'results/sv_harmonized/S1'
    native = canonical/'S1.harmonized.sv.vcf.gz'
    outputs={'vcf':native,'index':Path(str(native)+'.tbi'),
             'evidence':canonical/'S1.sv.evidence.tsv','provenance':canonical/'S1.phase9.provenance.json'}
    paths=phase9_bnd_paths(canonical.parent,'S1')
    outputs.update({'bnd_'+k:v for k,v in paths.items() if k!='directory'})
    source=frozen.root/'results/sv/S1.sawfish.sv.vcf.gz'
    pair={'vcf':source,'index':Path(str(source)+'.tbi')}
    with receipts.attempt_step('sv.harmonization',('S1',),outputs=outputs,
                              inputs={'sawfish':pair},input_manifest=frozen.root/'effective.yaml') as receipt:
        text=gzip.decompress(_vcf('S1')).decode().replace('##fileformat=VCFv4.3','##fileformat=VCFv4.4',1)
        text+='chr1\t100\tb1\tA\tA[chr1:500[\t.\tPASS\tSVTYPE=BND\tGT\t1|0\n'
        native_bytes=_bgzf(text.encode())
        data={str(native):base64.b64encode(native_bytes).decode(),
              str(outputs['index']):base64.b64encode(_bgzf(b'TBI'+bytes([1])+b'synthetic-index')).decode(),
              str(outputs['evidence']):base64.b64encode(b'fixture evidence, not scientific confidence').decode(),
              str(outputs['provenance']):base64.b64encode(b'{"synthetic":true}').decode()}
        CommandRunner().run([sys.executable,'-I','-B','-c',
            'import json,base64,sys;from pathlib import Path;'+
            '[(Path(p).parent.mkdir(parents=True,exist_ok=True),Path(p).write_bytes(base64.b64decode(b))) for p,b in json.loads(sys.argv[1]).items()]',
            json.dumps(data)],inherit_env=False,env={'SystemRoot':os.environ.get('SystemRoot','C:/Windows')} if os.name=='nt' else {},timeout=20)
        paths['native_vcf'].parent.mkdir(parents=True)
        for role,original in (('native_vcf',native),('native_index',outputs['index']),('native_evidence',outputs['evidence'])):
            paths[role].write_bytes(original.read_bytes())
        paths['native_raw'].write_text(text,encoding='utf-8')
        paths['native_result'].write_text('{"synthetic":true}',encoding='utf-8')
        policy=BndPolicy.from_config(frozen.config['sv']['harmonization']['bnd_representation'])
        plain=paths['directory']/'derived.plain.vcf'
        work=paths['directory']/'work';work.mkdir()
        derive_record_stream(native,plain,paths['ledger'],paths['junctions'],reference_lengths={'chr1':1000},
            sample_id='S1',native_sha256=compute_sha256(native),policy=policy,work_directory=work)
        paths['vcf'].write_bytes(_bgzf(plain.read_bytes()))
        paths['index'].write_bytes(outputs['index'].read_bytes())
        ref=frozen.context.reference
        before={str(p.resolve()):receipts.artifact_identity(p) for p in (ref.fasta,ref.fai)}
        copies={role:paths['native_'+role].relative_to(paths['directory']).as_posix() for role in ('vcf','index','raw','evidence')}
        artifacts={p.relative_to(paths['directory']).as_posix():receipts.artifact_identity(p)
                   for role,p in paths.items() if role not in ('directory','manifest')}
        manifest={'schema_version':1,'status':'completed','policy':policy.to_dict(),'native_output_role_unchanged':True,
                  'native_binary_copies':copies,'native_inputs_before':before,'native_inputs_after':before,'artifact_identities':artifacts}
        paths['manifest'].write_text(json.dumps(manifest),encoding='utf-8')
        receipt.complete(outputs=outputs,tool_versions={'jasmine':'1.synthetic'},validation={'validated':True},
            metadata={'reference_build':'GRCh38','bnd_representation':policy.to_dict()})
    return outputs


def generated_case(fixture,*,role='bnd_vcf'):
    config=fixture.config
    config['cohort']['enabled']=False
    config['cohort']['sv']['enabled']=False
    config['small']['enabled']=False
    for caller in ('sniffles2','pbsv','cutesv'):
        config['sv'][caller]['enabled']=False
    config['sv']['harmonization']['enabled']=True
    config['sv']['harmonization']['bnd_representation']=BndPolicy(POLICY['policy_id'],1,POLICY_SHA256,POLICY['mode']).to_dict()
    config['annotation'].update(enabled=True,overlap_enabled=True,gene_version='synthetic-1')
    bed=fixture.root.parent/'genes.bed';bed.write_text('chr1\t99\t100\tLOCAL\nchr1\t499\t500\tMATE\n',encoding='utf-8')
    config['annotation']['gene_bed']=str(bed)
    config['reference'].update(fasta=str(fixture.context.reference.fasta),build='GRCh38')
    row=_row('S1',producer='sv.harmonization',role=role);row['variant_ids']=['b1']
    frozen=fixture({'schema_version':1,'annotation':[row]})
    _seal(frozen,'S1',output_path=frozen.root/'results/sv/S1.sawfish.sv.vcf.gz')
    return frozen,seal_bundle(frozen)


def activate_consumer(frozen,record):
    (frozen.attempt/'status.json').write_text(json.dumps({'status':'DOWNSTREAM_DAG_RUNNING','dry_run':False,
        'downstream_materialization_sha256':receipts.artifact_identity(frozen.attempt/'downstream/materialization.json')['sha256']}))


@pytest.mark.parametrize('role',['vcf','bnd_vcf'])
def test_generated_native_or_derived_keeps_explicit_role_and_current_attempt(handoff_fixture,role):
    frozen,outputs=generated_case(handoff_fixture,role=role)
    watched={str(p):compute_sha256(p) for p in outputs.values()}
    effective_path,record=materialize_downstream(frozen.config,frozen.context,frozen.root,frozen.attempt)
    import yaml
    effective=yaml.safe_load(effective_path.read_text())
    assert effective["sv"]["harmonization"]["enabled"] is False
    assert effective["sv"]["harmonization"]["bnd_representation"]=={"enabled":False}
    assert effective["sv"]["harmonization"]["input_identity"]=={"enabled":False}
    assert frozen.config["sv"]["harmonization"]["enabled"] is True
    with (frozen.attempt/'downstream/annotation.inputs.tsv').open() as f: row=next(csv.DictReader(f,delimiter='	'))
    assert row['source_vcf']==str(outputs[role])
    assert bool(row.get('bnd_bundle_manifest'))==(role=='bnd_vcf')
    if role=='bnd_vcf': assert row['bnd_bundle_manifest']==str(outputs['bnd_manifest'])
    assert record['source_receipts'][0]['output_role']==role
    activate_consumer(frozen,record)
    from hifivar.workflow_report import _materialized_controls, _annotation_sources
    controls=json.loads((frozen.root/"launch_manifest.json").read_text())["control_sha256"]
    verified=_materialized_controls(frozen.root,frozen.attempt,frozen.config,controls)
    assert verified[0]==frozen.attempt/"downstream"
    source_rows,pending=_annotation_sources(frozen.root,frozen.attempt,frozen.config,controls,frozen.context.sample_ids)
    assert pending==[] and len(source_rows)==1
    assert bool(source_rows[0].get("bnd_bundle_manifest"))==(role=="bnd_vcf")
    inputs={'S1':{'vcf':outputs[role],'index':outputs['bnd_index' if role=='bnd_vcf' else 'index']}}
    if role=='bnd_vcf': inputs['S1']['bnd_bundle_manifest']=outputs['bnd_manifest']
    out=frozen.root/'results/annotation/S1/sv.overlap.tsv'
    with receipts.attempt_step('annotation.overlap',('S1',),outputs={'tsv':out},unit_id='sv',
            inputs=inputs,input_manifest=frozen.attempt/'downstream/annotation.inputs.tsv') as receipt:
        out.parent.mkdir(parents=True);out.write_text('synthetic overlap table, real algorithm tested separately')
        receipt.complete(outputs={'tsv':out},tool_versions={'hifivar':'unit'},validation={'validated':True},execution_kind='IN_PROCESS')
    result=receipts.read_step_receipt(frozen.root,frozen.attempt.name,'annotation.overlap',('S1',),unit_id='sv')
    binding=result['consumer_binding']
    assert binding['execution_stage']=='DOWNSTREAM_DAG_RUNNING'
    assert bool(binding['junction_bundle_bindings'])==(role=='bnd_vcf')
    assert bool(binding['producer_receipts'][0].get('bnd_bundle_manifest'))==(role=='bnd_vcf')
    assert watched=={str(p):compute_sha256(p) for p in outputs.values()}


@pytest.mark.parametrize('change',['index_role','bundle_selector','manifest_source','independent_producer'])
def test_report_does_not_accept_inconsistent_derived_role_or_reenabled_producer(handoff_fixture,change):
    import yaml
    from hifivar.workflow_report import _materialized_controls
    frozen,_=generated_case(handoff_fixture)
    effective,record=materialize_downstream(frozen.config,frozen.context,frozen.root,frozen.attempt)
    if change=='index_role':record['source_receipts'][0]['index_role']='index'
    elif change=='bundle_selector':record['source_receipts'][0]['bnd_bundle_manifest_role']=None
    elif change=='manifest_source':
        manifest=frozen.attempt/'downstream/annotation.inputs.tsv'
        with manifest.open() as f:data=list(csv.DictReader(f,delimiter='\t'))
        data[0]['bnd_bundle_manifest']=str(frozen.root/'other.json')
        with manifest.open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(data[0]),delimiter='\t');writer.writeheader();writer.writerows(data)
        record['controls']['annotation.inputs.tsv']=receipts.artifact_identity(manifest)
    else:
        config=yaml.safe_load(effective.read_text());config['sv']['harmonization']['enabled']=True
        effective.write_text(yaml.safe_dump(config));record['controls']['effective.yaml']=receipts.artifact_identity(effective)
    (frozen.attempt/'downstream/materialization.json').write_text(json.dumps(record))
    activate_consumer(frozen,record)
    controls=json.loads((frozen.root/'launch_manifest.json').read_text())['control_sha256']
    with pytest.raises(InputValidationError):_materialized_controls(frozen.root,frozen.attempt,frozen.config,controls)


@pytest.mark.parametrize('damage',['missing_member','changed_ledger','missing_receipt','wrong_policy','wrong_manifest_role','wrong_native_copy'])
def test_incomplete_or_substituted_bundle_does_not_become_ready(handoff_fixture,damage):
    frozen,outputs=generated_case(handoff_fixture)
    if damage=='changed_ledger':
        outputs['bnd_ledger'].write_text('changed')
        with pytest.raises(OutputValidationError,match="output bytes differ"): materialize_downstream(frozen.config,frozen.context,frozen.root,frozen.attempt)
        assert not (frozen.attempt/'downstream').exists()
        return
    record=receipts.read_step_receipt(frozen.root,frozen.attempt.name,'sv.harmonization',('S1',))
    record=copy.deepcopy(record)
    if damage=='missing_member': record['outputs'].pop('bnd_junctions')
    elif damage=='missing_receipt': record['track']='sv.sawfish'
    elif damage=='wrong_policy': record['metadata']['bnd_representation']={'enabled':False}
    elif damage=='wrong_manifest_role': record['outputs']['bnd_manifest']=record['outputs']['provenance']
    else: record['outputs']['vcf']=record['outputs']['bnd_vcf']
    with pytest.raises(InputValidationError): generated_bnd_roles(record,frozen.context.reference,'S1')


def test_derived_handoff_refuses_another_attempt_even_when_output_files_exist(handoff_fixture):
    frozen,_=generated_case(handoff_fixture)
    other=frozen.root/'attempts/attempt_other';other.mkdir()
    with pytest.raises(InputValidationError): materialize_downstream(frozen.config,frozen.context,frozen.root,other)
    assert not (other/'downstream').exists()

@pytest.mark.parametrize('change',['index_role','output_role','bnd_bundle_manifest_role'])
def test_consumer_rejects_coherently_rewritten_materialization_role(handoff_fixture,change):
    frozen,outputs=generated_case(handoff_fixture)
    _,record=materialize_downstream(frozen.config,frozen.context,frozen.root,frozen.attempt)
    record['source_receipts'][0][change]={'index_role':'index','output_role':'vcf','bnd_bundle_manifest_role':None}[change]
    material=frozen.attempt/'downstream/materialization.json'
    material.write_text(json.dumps(record),encoding='utf-8')
    activate_consumer(frozen,record)
    with pytest.raises(InputValidationError):
        with receipts.attempt_step('annotation.overlap',('S1',),outputs={'tsv':frozen.root/'results/bad.tsv'},unit_id='sv',
                inputs={'S1':{'vcf':outputs['bnd_vcf'],'index':outputs['bnd_index'],'bnd_bundle_manifest':outputs['bnd_manifest']}},
                input_manifest=frozen.attempt/'downstream/annotation.inputs.tsv'):
            pytest.fail('Changed source role reached execution')


def test_external_plan_freezes_all_bundle_members_not_only_manifest(handoff_fixture):
    import yaml
    from hifivar.downstream_handoff import prepare_downstream_intent
    frozen,outputs=generated_case(handoff_fixture)
    config=copy.deepcopy(frozen.config);config['sv']['enabled']=False
    new_root=frozen.root.parent/'external_consumer'
    source=frozen.root.parent/'external_plan.yaml'
    row={'sample_id':'S1','track':'sv','origin':'EXTERNAL','path':str(outputs['bnd_vcf']),
         'index_path':str(outputs['bnd_index']),'bnd_bundle_manifest':str(outputs['bnd_manifest']),
         'source_tool':'jasmine','source_version':'1.synthetic','reference_build':'GRCh38','variant_ids':['b1']}
    source.write_text(yaml.safe_dump({'schema_version':1,'annotation':[row]}),encoding='utf-8')
    config['workflow']['downstream_source_plan']=str(source)
    controls,inputs,_=prepare_downstream_intent(copy.deepcopy(config),frozen.context,source.parent,new_root)
    assert not new_root.exists()
    assert set(outputs[k] for k in outputs if k.startswith('bnd_')).issubset(set(inputs))
    normalized=yaml.safe_load(controls['downstream_plan.yaml'])
    assert normalized['annotation'][0]['bnd_bundle_manifest']==str(outputs['bnd_manifest'].resolve())
    outputs['bnd_ledger'].write_text('changed after explicit selection',encoding='utf-8')
    with pytest.raises(InputValidationError,match='artifact bytes changed'):
        prepare_downstream_intent(copy.deepcopy(config),frozen.context,source.parent,new_root)
    assert not new_root.exists()
