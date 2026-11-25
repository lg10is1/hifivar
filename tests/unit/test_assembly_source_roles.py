"""Assembly roles are input provenance, not per-variant or biological phasing truth."""
import json
import os
from pathlib import Path
import sys
from types import SimpleNamespace
import pytest
import yaml
from hifivar import launcher
from hifivar.assembly import AssemblyRole, HaplotypeAssemblyArtifact
from hifivar.assembly_source_roles import assembly_role_binding, source_haplotype_binding
from hifivar.assembly_sv import SVEvidenceSource
from hifivar.attempt_receipts import ATTEMPT_CONTEXT_ENV, artifact_identity, attempt_step, create_attempt_context, read_step_receipt
from hifivar.command import CommandRunner
from hifivar.exceptions import InputValidationError
from hifivar.harmonization import EvidenceRunStatus, SVEvidenceSourceArtifact


def native_files(outputs):
    env={k:v for k,v in os.environ.items() if k.upper()=='SYSTEMROOT'}
    CommandRunner().run([sys.executable,'-I','-B','-c',
        'import pathlib,sys;[(pathlib.Path(p).parent.mkdir(parents=True,exist_ok=True),pathlib.Path(p).write_bytes(b">unit\\nACGT\\n")) for p in sys.argv[1:]]',
        *map(str,outputs.values())],env=env,inherit_env=False)


def complete(outputs,receipt,metadata=None):
    native_files(outputs)
    receipt.complete(outputs=outputs,tool_versions={'synthetic-process':'1'},
        validation={'validated':True,'scope':'SYNTHETIC_BYTES_ONLY_NOT_BIOLOGICAL_TOOL_VALIDATION'},metadata=metadata or {})


@pytest.fixture
def owned(tmp_path,monkeypatch):
    # Synthetic controller fixture; actual native acceptance is a separate gate.
    root=tmp_path/'run';root.mkdir()
    config={'assembly':{'enabled':True},'assembly_sv':{'enabled':True,
        'pav':{'enabled':True},'svim_asm':{'enabled':True}},'sv':{'enabled':False,
        'harmonization':{'enabled':True,'input_vcfs':{}}}}
    sheet='sample_id\tinput\nS1\t/context.fastq\n'
    (root/'samples.tsv').write_text(sheet,encoding='utf-8')
    (root/'effective.yaml').write_text(yaml.safe_dump(config),encoding='utf-8')
    runtime={'payload_sha256':'a'*64,'hifivar_version':'synthetic-fixture'}
    monkeypatch.setattr(launcher,'runtime_identity',lambda _:runtime)
    import hifivar.package_resources
    monkeypatch.setattr(hifivar.package_resources,'installed_workflow_root',lambda:root)
    manifest={'identity':{'schema':2,'run_dir':str(root),'config':config,'sample_sheet':sheet,'runtime':runtime,'inputs':[]},
        'control_sha256':{name:artifact_identity(root/name)['sha256'] for name in ('effective.yaml','samples.tsv')}}
    (root/'launch_manifest.json').write_text(json.dumps(manifest),encoding='utf-8')
    attempt=root/'attempts/roles_test';attempt.mkdir(parents=True)
    (attempt/'status.json').write_text(json.dumps({'status':'DAG_RUNNING','dry_run':False}),encoding='utf-8')
    context=create_attempt_context(root,attempt,dry_run=False)
    monkeypatch.setenv(ATTEMPT_CONTEXT_ENV,str(context))
    outputs={role+'_fasta':root/'results/assembly'/('S1.'+role+'.fasta') for role in ('primary','haplotype1','haplotype2')}
    with attempt_step('assembly',('S1',),outputs=outputs) as receipt:complete(outputs,receipt)
    inputs=tuple(HaplotypeAssemblyArtifact('S1',role,outputs[role.value+'_fasta'],outputs[role.value+'_fasta'],
        'synthetic',outputs[role.value+'_fasta'].stat().st_size) for role in (AssemblyRole.HAPLOTYPE1,AssemblyRole.HAPLOTYPE2))
    return SimpleNamespace(root=root,attempt=attempt,context=context,assemblies=inputs)


def test_external_and_unlinked_assembly_roles_remain_unknown():
    for caller in ('pav','svim_asm'):
        value=source_haplotype_binding('S1',caller,external=True,context_path=None)
        assert value['roles']==[] and value['binding'] is None and value['status'].startswith('UNVERIFIED')
        assert source_haplotype_binding('S1',caller,external=False,context_path=None)['roles']==[]
    assert source_haplotype_binding('S1','sawfish',external=True,context_path=None)['status']=='NOT_APPLICABLE_READ_SOURCE'
    with pytest.raises(InputValidationError):source_haplotype_binding('S1','unknown',external=True,context_path=None)


def test_explicit_roles_bind_actual_parent_output_files(owned):
    value=assembly_role_binding('S1',owned.assemblies,owned.context)
    assert value['phase_accuracy']=='NOT_ASSESSED'
    assert [x['role'] for x in value['assemblies']]==['haplotype1','haplotype2']
    assert value['assemblies'][0]['fasta']==artifact_identity(owned.assemblies[0].path)
    assert value['assembly_receipt']['sha256']==artifact_identity(Path(value['assembly_receipt']['path']))['sha256']
    assert assembly_role_binding('S1',owned.assemblies,None) is None


def source_receipt(owned,caller,binding):
    outputs={'vcf':owned.root/'results/assembly_sv/S1'/f'S1.{caller}.assembly.sv.vcf.gz'}
    outputs['index']=Path(str(outputs['vcf'])+'.tbi')
    with attempt_step('assembly_sv.'+caller,('S1',),outputs=outputs) as receipt:
        complete(outputs,receipt,{'assembly_role_binding':binding})
    return outputs


@pytest.mark.parametrize('caller',['pav','svim_asm'])
def test_generated_source_resolves_only_accepted_same_attempt_roles(owned,caller):
    value=assembly_role_binding('S1',owned.assemblies,owned.context)
    source_receipt(owned,caller,value)
    actual=source_haplotype_binding('S1',caller,external=False,context_path=owned.context)
    assert actual['roles']==['haplotype1','haplotype2'] and actual['binding']==value
    assert 'NOT_VARIANT_SUPPORT' in actual['status']
    owned.assemblies[0].path.write_bytes(b'altered')
    with pytest.raises(InputValidationError):source_haplotype_binding('S1',caller,external=False,context_path=owned.context)


def test_single_explicit_input_does_not_invent_second_haplotype(owned):
    value=assembly_role_binding('S1',owned.assemblies[:1],owned.context)
    source_receipt(owned,'pav',value)
    assert source_haplotype_binding('S1','pav',external=False,context_path=owned.context)['roles']==['haplotype1']


@pytest.mark.parametrize('change',['missing_metadata','bool_schema','wrong_sample','swapped_files','duplicate_role','changed_upstream_receipt'])
def test_malformed_or_changed_role_lineage_fails_closed(owned,change):
    value=assembly_role_binding('S1',owned.assemblies,owned.context)
    if change=='bool_schema':value['schema_version']=True
    elif change=='wrong_sample':value['sample_id']='S2'
    elif change=='swapped_files':value['assemblies'][0]['fasta'],value['assemblies'][1]['fasta']=value['assemblies'][1]['fasta'],value['assemblies'][0]['fasta']
    elif change=='duplicate_role':value['assemblies'][1]['role']='haplotype1'
    elif change=='changed_upstream_receipt':Path(value['assembly_receipt']['path']).write_text('{}',encoding='utf-8')
    source_receipt(owned,'pav',None if change=='missing_metadata' else value)
    with pytest.raises(InputValidationError):source_haplotype_binding('S1','pav',external=False,context_path=owned.context)


def test_absent_assembly_sv_producer_cannot_be_inferred(owned):
    with pytest.raises(InputValidationError,match='no current-attempt'):
        source_haplotype_binding('S1','pav',external=False,context_path=owned.context)


def test_explicit_roles_must_refer_to_matching_sample_and_parent(owned):
    with pytest.raises(InputValidationError):assembly_role_binding('S2',owned.assemblies,owned.context)
    x=owned.assemblies[0]
    other=HaplotypeAssemblyArtifact('S1',AssemblyRole.HAPLOTYPE2,x.path,x.path,'synthetic',x.file_size)
    with pytest.raises(InputValidationError):assembly_role_binding('S1',(x,other),owned.context)
    with pytest.raises(InputValidationError):assembly_role_binding('S1',(),owned.context)


def test_source_json_distinguishes_unknown_declared_and_not_applicable():
    def source(kind,roles=()):return SVEvidenceSourceArtifact('S1',kind,'pav' if kind is SVEvidenceSource.ASSEMBLY else 'sawfish',None,None,EvidenceRunStatus.COMPLETED,roles)
    assert source(SVEvidenceSource.ASSEMBLY).to_dict()['haplotype_role_status']=='UNVERIFIED_ASSEMBLY_INPUT_ROLES'
    assert source(SVEvidenceSource.READ).to_dict()['haplotype_role_status']=='NOT_APPLICABLE_READ_SOURCE'
    assert 'NOT_VARIANT_PHASE_ACCURACY' in source(SVEvidenceSource.ASSEMBLY,('haplotype1',)).to_dict()['haplotype_role_status']


def test_consumer_seal_and_readback_recheck_assembly_file_lineage(owned):
    value=assembly_role_binding('S1',owned.assemblies,owned.context)
    inputs={caller:source_receipt(owned,caller,value) for caller in ('pav','svim_asm')}
    outputs={role:owned.root/'results/harmonized'/role for role in ('vcf','index','evidence','provenance')}
    with attempt_step('sv.harmonization',('S1',),outputs=outputs,inputs=inputs,input_manifest=owned.root/'effective.yaml') as receipt:
        complete(outputs,receipt)
    record=read_step_receipt(owned.root,owned.attempt.name,'sv.harmonization',('S1',))
    assert record['consumer_binding']['source_haplotype_roles']['pav']['roles']==['haplotype1','haplotype2']
    assert record['consumer_binding']['source_haplotype_roles']['svim_asm']['binding']==value
    owned.assemblies[1].path.write_bytes(b'changed despite unchanged final VCF')
    with pytest.raises(InputValidationError):read_step_receipt(owned.root,owned.attempt.name,'sv.harmonization',('S1',))


@pytest.mark.parametrize('invalid',[None,'haplotype1',[{}]])
def test_invalid_artifact_objects_report_explicit_errors(owned,invalid):
    with pytest.raises(InputValidationError,match='explicit assembly artifact objects'):
        assembly_role_binding('S1',invalid,owned.context)
