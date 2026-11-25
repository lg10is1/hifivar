"""Synthetic source/receipt contracts, not native BAM QC or human validation."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml

from hifivar.command import CommandRunner
from hifivar.downstream_plan import read_downstream_plan
from hifivar.exceptions import InputValidationError, OutputValidationError, ReferenceError
from hifivar.reference import ReferenceGenome
from hifivar.review import EvidenceClass, VariantClass
from hifivar.review_sources import ReviewSite, read_review_source_plan, resolve_review_targets
from hifivar.sample import InputType
from test_downstream_handoff import _bgzf, _seal, handoff_fixture

DEFAULT = Path(__file__).parents[2] / 'src/hifivar/resources/configs/default.yaml'


def _context(reference, directory, samples=('S1',)):
    entries = []
    for sample in samples:
        bam = directory / (sample + '.bam')
        bam.write_bytes(b'SYNTHETIC_BAM_NOT_NATIVE_QC')
        Path(str(bam)+'.bai').write_bytes(b'SYNTHETIC_BAI_NOT_NATIVE_QC')
        entries.append(SimpleNamespace(sample=SimpleNamespace(sample_id=sample,
            input=SimpleNamespace(input_type=InputType.BAM, files=(bam,)))))
    return SimpleNamespace(reference=reference, sample_ids=samples, samples=tuple(entries), n_samples=len(samples))


@pytest.fixture
def case(tmp_path):
    ref = tmp_path / 'reference.fa'
    ref.write_text('>chr1\n'+'A'*1000+'\n>chr2\n'+'C'*1000+'\n', encoding='utf-8')
    Path(str(ref)+'.fai').write_text('chr1\t1000\t6\t1000\t1001\nchr2\t1000\t1013\t1000\t1001\n', encoding='utf-8')
    context = _context(ReferenceGenome.from_fasta(ref, build='GRCh38'), tmp_path)
    config = yaml.safe_load(DEFAULT.read_text(encoding='utf-8'))
    config['review']['enabled'] = True
    config['sv']['enabled'] = config['small']['enabled'] = True
    return SimpleNamespace(root=tmp_path, context=context, config=config)


def _vcf(case, records, *, sample='S1', build='GRCh38', name='source'):
    path = case.root / (name+'.vcf.gz')
    header = ('##fileformat=VCFv4.4\n##reference='+build+'\n'
        '##contig=<ID=chr1,length=1000>\n##contig=<ID=chr2,length=1000>\n'
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">\n'
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End">\n'
        '##INFO=<ID=TRID,Number=1,Type=String,Description="Locus">\n'
        '##INFO=<ID=MOTIFS,Number=.,Type=String,Description="Motifs">\n'
        '##INFO=<ID=STRUC,Number=1,Type=String,Description="Structure">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+sample+'\n')
    path.write_bytes(_bgzf((header+''.join(records)).encode()))
    Path(str(path)+'.tbi').write_bytes(_bgzf(b'TBI\x01synthetic-structural-index'))
    return path


def _record(identifier='v1', *, pos=100, ref='A', alt='C', info='.'):
    return f'chr1\t{pos}\t{identifier}\t{ref}\t{alt}\t.\tPASS\t{info}\tGT\t0|1\n'


def _candidate(identifier='v1', kind='SNV', review_id='R1', **extras):
    return dict(review_id=review_id, variant_id=identifier, variant_type=kind, **extras)


def _external(path, *, track='small', candidates=None, **extra):
    return dict(sample_id='S1', track=track, origin='EXTERNAL', path=str(path),
        index_path=str(path)+'.tbi', source_tool='explicit-caller', source_version='1.0-declared',
        reference_build='GRCh38', candidates=candidates or [_candidate()], **extra)


def _plan(case, rows, *, text=None):
    path = case.root/'review_sources.yaml'
    path.write_text(text if text is not None else yaml.safe_dump({'schema_version':1,'review':rows}), encoding='utf-8')
    return path


def _read(case, path):
    return read_review_source_plan(path, sample_ids=case.context.sample_ids,
                                  reference_build='GRCh38', config=case.config)


def _resolve(case, rows):
    plan = _read(case, _plan(case, rows))
    return resolve_review_targets(plan, config=case.config, context=case.context,
                                  output_directory=case.root/'review_out')


def test_explicit_small_source_preserves_raw_identity_and_creates_nothing(case, monkeypatch):
    vcf = _vcf(case, [_record()])
    watched = [vcf,Path(str(vcf)+'.tbi'),case.context.reference.fasta,case.context.reference.fai]
    before = {str(p):p.read_bytes() for p in watched}
    monkeypatch.setattr(CommandRunner,'run',lambda *a,**k:pytest.fail('Resolver must not execute a command'))
    result = _resolve(case, [_external(vcf)])
    assert result.targets[0].variant_id == 'v1'
    assert result.targets[0].evidence_class is EvidenceClass.EXPLICIT
    binding = result.source_bindings[0]
    assert binding['source_version']=='1.0-declared' and binding['source_receipt'] is None
    assert binding['version_source']=='USER_DECLARED_NOT_REPROBED'
    assert binding['selected_records'][0]['record_id']=='v1'
    assert binding['selected_records'][0]['selected_record_sha256']==hashlib.sha256(_record().encode()).hexdigest()
    assert before=={str(p):p.read_bytes() for p in watched}
    assert not (case.root/'review_out').exists()
    data=result.to_dict()
    assert data['review_status']=='NOT_REVIEWED' and data['IGV_executed'] is False
    assert data['accuracy']=='NOT_ASSESSED' and not data['caller_support_is_confidence']
    data['source_bindings'][0]['selected_records'][0]['record_id']='CHANGED_EXPORT_ONLY'
    assert result.to_dict()['source_bindings'][0]['selected_records'][0]['record_id']=='v1'


@pytest.mark.parametrize('ref,alt,kind,end', [('A','C','SNV',100),('AA','A','INDEL',101),
    ('A','AT','INDEL',100),('A','C,G','SNV',100),('A','C,AT','INDEL',100)])
def test_record_level_small_alleles_and_reference_span(case, ref,alt,kind,end):
    vcf=_vcf(case,[_record(ref=ref,alt=alt)])
    target=_resolve(case,[_external(vcf,candidates=[_candidate(kind=kind)])]).targets[0]
    assert (target.start,target.end)==(100,end)
    assert target.variant_class.value==kind


@pytest.mark.parametrize('ref,alt,kind,info', [('AA','CC','INDEL','.'),('AA','A','SNV','.'),
    ('A','<NON_REF>','INDEL','END=120'),('A','*','SNV','.'),('A','A','SNV','.'),
    ('AA','A','INDEL','END=102')])
def test_small_class_is_not_coerced(case,ref,alt,kind,info):
    vcf=_vcf(case,[_record(ref=ref,alt=alt,info=info)])
    with pytest.raises(InputValidationError):
        _resolve(case,[_external(vcf,candidates=[_candidate(kind=kind)])])


@pytest.mark.parametrize('kind,end', [('DEL',150),('DUP',160),('INV',180),('INS',100)])
def test_sv_source_type_and_interval(case,kind,end):
    vcf=_vcf(case,[_record(alt='<'+kind+'>',info=f'SVTYPE={kind};END={end}')])
    target=_resolve(case,[_external(vcf,track='sv',candidates=[_candidate(kind=kind)])]).targets[0]
    assert (target.start,target.end)==(100,end)


@pytest.mark.parametrize('info,kind', [('SVTYPE=DEL','DEL'),('SVTYPE=DEL;END=100','DEL'),
    ('SVTYPE=INV;END=150','DUP'),('END=150','DEL')])
def test_sv_missing_span_or_mismatched_class_fails(case,info,kind):
    vcf=_vcf(case,[_record(alt='<DEL>',info=info)])
    with pytest.raises(InputValidationError):
        _resolve(case,[_external(vcf,track='sv',candidates=[_candidate(kind=kind)])])


@pytest.mark.parametrize('alt', ['A[chr2:500[','A]chr2:500]','[chr2:500[A',']chr2:500]A'])
def test_bnd_uses_two_anchors_and_preserves_alt(case,alt):
    vcf=_vcf(case,[_record('b1',alt=alt,info='SVTYPE=BND;END=900')])
    target=_resolve(case,[_external(vcf,track='sv',candidates=[_candidate('b1','BND')])]).targets[0]
    assert [(x.contig,x.variant_start,x.variant_end) for x in target.loci]==[('chr1',100,100),('chr2',500,500)]
    assert target.to_dict()['junction_record']['alt']==alt


@pytest.mark.parametrize('pos,ref,alt,count',[(100,'A','A.',1),(100,'A','.A',1),
    (0,'N','.[chr2:500[',1),(1001,'N',']chr2:0].',0)])
def test_single_virtual_bnd_is_not_guessed(case,pos,ref,alt,count):
    vcf=_vcf(case,[_record('b1',pos=pos,ref=ref,alt=alt,info='SVTYPE=BND')])
    target=_resolve(case,[_external(vcf,track='sv',candidates=[_candidate('b1','BND')])]).targets[0]
    assert len(target.loci)==count and target.unrendered_bnd_endpoints


def test_trid_namespace_and_repeat_interval_exclude_anchor(case):
    vcf=_vcf(case,[_record('.',ref='AAAAA',alt='AAA',info='TRID=L1;MOTIFS=A;STRUC=(A)n;END=104')])
    row=_external(vcf,track='tr',candidates=[_candidate('TRID:L1','TR')],catalog_id='sha256:'+'c'*64)
    row['source_tool']='TRGT'
    target=_resolve(case,[row]).targets[0]
    assert target.variant_id=='TRID:L1' and (target.start,target.end)==(101,104)
    assert target.evidence_class is EvidenceClass.EXPLICIT


@pytest.mark.parametrize('error', ['wrong_tool','ordinary_id','missing_trid','missing_catalog'])
def test_tr_does_not_accept_ordinary_id_or_unbound_catalog(case,error):
    vcf=_vcf(case,[_record('v1',ref='AAAAA',alt='AAA',info='TRID=L1;MOTIFS=A;STRUC=(A)n;END=104')])
    row=_external(vcf,track='tr',candidates=[_candidate('TRID:L1','TR')],catalog_id='sha256:'+'c'*64)
    row['source_tool']='TRGT'
    if error=='wrong_tool':row['source_tool']='other'
    elif error=='ordinary_id':row['candidates'][0]['variant_id']='v1'
    elif error=='missing_trid':row['candidates'][0]['variant_id']='TRID:missing'
    else:row.pop('catalog_id')
    with pytest.raises(InputValidationError):_resolve(case,[row])


@pytest.mark.parametrize('identifier',['.','v1'])
def test_exact_site_selects_source_record_without_id_rewrite(case,identifier):
    vcf=_vcf(case,[_record(identifier)])
    candidate=dict(review_id='SITE1',variant_type='SNV',site=dict(contig='chr1',position=100,ref='A',alt='C'))
    result=_resolve(case,[_external(vcf,candidates=[candidate])])
    record=result.source_bindings[0]['selected_records'][0]
    assert record['raw_vcf_id']==identifier and record['selection_namespace']=='EXACT_SITE_RECORD'
    assert record['selector']==ReviewSite('chr1',100,'A','C').key
    assert result.targets[0].variant_id==record['selector']


@pytest.mark.parametrize('error',['wrong_alt','absent','duplicate','both','neither','bnd','tr','same_record_twice'])
def test_exact_site_ambiguity_and_namespace_fail_closed(case,error):
    records=[_record('v1')]
    site=dict(contig='chr1',position=100,ref='A',alt='C')
    candidate=dict(review_id='SITE1',variant_type='SNV',site=site)
    if error=='wrong_alt':site['alt']='G'
    elif error=='absent':site['position']=101
    elif error=='duplicate':records*=2
    elif error=='both':candidate['variant_id']='v1'
    elif error=='neither':candidate.pop('site')
    elif error in {'bnd','tr'}:candidate['variant_type']=error.upper()
    vcf=_vcf(case,records)
    candidates=[candidate]
    if error=='same_record_twice':candidates.append(_candidate(review_id='ID1'))
    with pytest.raises(InputValidationError):_resolve(case,[_external(vcf,candidates=candidates)])


@pytest.mark.parametrize('error',['sample','build','index','plain_gzip','missing_id','duplicate_id','type_category','disabled','unknown_field','state','variant_ids'])
def test_external_source_contract_negative(case,error):
    vcf=_vcf(case,[_record(),*([_record()] if error=='duplicate_id' else [])],
             sample='OTHER' if error=='sample' else 'S1',build='GRCh37' if error=='build' else 'GRCh38')
    row=_external(vcf)
    if error=='index':Path(str(vcf)+'.tbi').write_bytes(b'NOT_A_TBI')
    elif error=='plain_gzip':
        import gzip
        vcf.write_bytes(gzip.compress(gzip.decompress(vcf.read_bytes())))
    elif error=='missing_id':row['candidates'][0]['variant_id']='absent'
    elif error=='type_category':row['candidates'][0]['variant_type']='DEL'
    elif error=='disabled':case.config['review']['enabled']=False
    elif error=='unknown_field':row['confidence']='high'
    elif error=='state':row['state']='CALLED'
    elif error=='variant_ids':row['variant_ids']=['v1']
    with pytest.raises((InputValidationError,OutputValidationError,ReferenceError)):_resolve(case,[row])


@pytest.mark.parametrize('error',['duplicate_review_id','duplicate_source','changed_plan','wrong_context','missing_bai','raw_alignment','forged_plan'])
def test_resolution_plan_and_context_invariants(case,error):
    vcf=_vcf(case,[_record()]);row=_external(vcf)
    rows=[row]
    if error=='duplicate_review_id':row['candidates'].append(_candidate(review_id='R1',identifier='v2'))
    elif error=='duplicate_source':rows.append(dict(row,candidates=[_candidate(review_id='R2')]))
    path=_plan(case,rows)
    if error in {'duplicate_review_id','duplicate_source'}:
        with pytest.raises(InputValidationError):_read(case,path)
        return
    plan=_read(case,path)
    if error=='changed_plan':path.write_text(path.read_text()+'# changed\n',encoding='utf-8')
    elif error=='wrong_context':case.context.sample_ids=('OTHER',)
    elif error=='missing_bai':Path(str(case.context.samples[0].sample.input.files[0])+'.bai').unlink()
    elif error=='raw_alignment':case.context.samples[0].sample.input.input_type=InputType.FASTQ
    elif error=='forged_plan':
        from dataclasses import replace
        plan=replace(plan,sha256='f'*64)
    with pytest.raises((InputValidationError,OutputValidationError)):
        resolve_review_targets(plan,config=case.config,context=case.context,output_directory=case.root/'out')


@pytest.mark.parametrize('text',['schema_version: 1\nschema_version: 1\nreview: []\n',
    'schema_version: 1\nreview: &rows []\nextra: *rows\n','schema_version: true\nreview: []\n',
    'schema_version: 1\nreview: []\n','schema_version: 1\nreview: {}\n'])
def test_bounded_strict_yaml(case,text):
    with pytest.raises(InputValidationError):_read(case,_plan(case,[],text=text))


def test_generic_downstream_plan_still_rejects_review(case):
    vcf=_vcf(case,[_record()]);path=_plan(case,[_external(vcf)])
    with pytest.raises(InputValidationError,match='unknown keys'):
        read_downstream_plan(path,sample_ids=('S1',),reference_build='GRCh38',config=case.config)


@pytest.mark.parametrize('error',[None,'missing_receipt','wrong_attempt','changed_source','missing_attempt','gvcf','disabled_producer','future_path'])
def test_generated_source_requires_exact_attempt_native_receipt(handoff_fixture,tmp_path,error):
    frozen=handoff_fixture()
    frozen.config['review']['enabled']=True
    frozen.context.samples=_context(frozen.context.reference,tmp_path,('S2','S1')).samples
    if error!='missing_receipt':outputs=_seal(frozen,'S1',producer='small.deepvariant',version='1.10.0-unit')
    row=dict(sample_id='S1',track='small',origin='GENERATED',producer='small.deepvariant',output_role='vcf',candidates=[_candidate()])
    if error=='gvcf':row['output_role']='gvcf'
    elif error=='disabled_producer':frozen.config['small']['enabled']=False
    elif error=='future_path':row['path']=str(tmp_path/'future.vcf.gz')
    case=SimpleNamespace(root=tmp_path,context=frozen.context,config=frozen.config)
    if error in {'gvcf','disabled_producer','future_path'}:
        with pytest.raises(InputValidationError):_read(case,_plan(case,[row]))
        return
    plan=_read(case,_plan(case,[row]))
    attempt=frozen.attempt
    if error=='wrong_attempt':
        attempt=frozen.root/'attempts'/'attempt_other';attempt.mkdir()
        (attempt/'status.json').write_text(json.dumps({'status':'PRODUCER_DAG_RUNNING','dry_run':False}))
    elif error=='changed_source':outputs['vcf'].write_bytes(outputs['vcf'].read_bytes()+b'x')
    options=dict(run_root=frozen.root,attempt_dir=attempt) if error!='missing_attempt' else {}
    if error:
        with pytest.raises((InputValidationError,OutputValidationError)):
            resolve_review_targets(plan,config=case.config,context=case.context,output_directory=tmp_path/'out',**options)
    else:
        result=resolve_review_targets(plan,config=case.config,context=case.context,output_directory=tmp_path/'out',**options)
        binding=result.source_bindings[0]
        assert binding['source_receipt']['producer_attempt_id']==frozen.attempt.name
        assert binding['version_source']=='SELECTED_ATTEMPT_PRODUCER_RECEIPT'
        assert binding['source_version']=='1.10.0-unit'
        assert result.targets[0].evidence_class is EvidenceClass.SMALL_VARIANT
        assert result.to_dict()['main_launcher_integration']=='NOT_IMPLEMENTED_BY_THIS_LIBRARY'


@pytest.mark.parametrize('position',[0,-1,True,'100',1.5])
def test_exact_site_position_is_strict(case,position):
    vcf=_vcf(case,[_record()])
    candidate=dict(review_id='R1',variant_type='SNV',site=dict(contig='chr1',position=position,ref='A',alt='C'))
    with pytest.raises(InputValidationError):_resolve(case,[_external(vcf,candidates=[candidate])])


@pytest.mark.parametrize('key,value',[('review_id','bad id'),('review_id','../R1'),
    ('variant_id','.'),('variant_id','v1;v2'),('variant_type','MNV')])
def test_candidate_fields_are_explicit_and_not_path_injections(case,key,value):
    vcf=_vcf(case,[_record()]);candidate=_candidate();candidate[key]=value
    with pytest.raises(InputValidationError):_resolve(case,[_external(vcf,candidates=[candidate])])


def test_index_mutation_during_selection_is_not_promoted(case,monkeypatch):
    import hifivar.review_sources as sources
    vcf=_vcf(case,[_record()]);original=sources._fields_for_records
    def changed(path,records):
        result=original(path,records)
        index=Path(str(path)+'.tbi');index.write_bytes(index.read_bytes()+b'mutated synthetic index')
        return result
    monkeypatch.setattr(sources,'_fields_for_records',changed)
    with pytest.raises(OutputValidationError,match='changed'):_resolve(case,[_external(vcf)])
    assert not (case.root/'review_out').exists()


@pytest.mark.parametrize('producer,track,kind',[('sv.sawfish','sv','DEL'),('sv.sniffles2','sv','DEL'),
    ('sv.pbsv','sv','DEL'),('sv.cutesv','sv','DEL'),('sv.harmonization','sv','DEL'),
    ('small.deepvariant','small','SNV'),('phasing.hiphase','small','SNV'),
    ('assembly_sv.pav','sv','DEL'),('assembly_sv.svim_asm','sv','DEL')])
def test_generated_source_intent_selectors_do_not_claim_execution(case,producer,track,kind):
    section,backend=producer.split('.')
    case.config[section]['enabled']=True
    if section=='sv':case.config[section][backend]['enabled']=True
    if section=='assembly_sv':
        case.config['assembly']['enabled']=False
        case.config['assembly_sv']['input_manifest']='explicit_ready_assemblies.tsv'
        case.config['assembly_sv'][backend]['enabled']=True
    row=dict(sample_id='S1',track=track,origin='GENERATED',producer=producer,output_role='vcf',
             candidates=[_candidate(kind=kind)])
    plan=_read(case,_plan(case,[row]))
    assert plan.rows[0].source.path is None and plan.rows[0].source.source_version is None
    assert not (case.root/'review_out').exists()


def test_tr_catalog_mismatch_fails_before_resolution(case):
    catalog=case.root/'catalog.bed';catalog.write_text('chr1\t100\t104\tL1\n')
    case.config['tr']['catalog']=str(catalog)
    vcf=_vcf(case,[_record('.',ref='AAAAA',alt='AAA',info='TRID=L1;MOTIFS=A;STRUC=(A)n;END=104')])
    row=_external(vcf,track='tr',candidates=[_candidate('TRID:L1','TR')],catalog_id='sha256:'+'f'*64)
    row['source_tool']='TRGT'
    with pytest.raises(InputValidationError,match='catalog identity'):_resolve(case,[row])


def test_source_yaml_limits_are_enforced(case):
    for text in ('# comment\n'*4097, 'x'*(1024*1024+1)):
        with pytest.raises(InputValidationError):_read(case,_plan(case,[],text=text))


def test_review_resolution_cannot_enable_main_launcher_implicitly(case):
    import hifivar.launcher as launcher
    from hifivar.exceptions import ConfigurationError
    assert 'review' not in launcher._UNSUPPORTED
    with pytest.raises(ConfigurationError,match='review.source_plan'):
        launcher._enabled_tracks(case.config)


@pytest.mark.parametrize('damage',[None,'missing_manifest','ledger','site_detach'])
def test_derived_bnd_review_requires_complete_source_bundle(case,damage):
    import gzip
    from test_junction_consumers import make_bundle
    genome,derived,bam,manifest=make_bundle(case.root)
    derived.write_bytes(_bgzf(gzip.decompress(derived.read_bytes())))
    index=Path(str(derived)+'.tbi');index.write_bytes(_bgzf(b'TBI\x01synthetic-structural-index'))
    metadata=json.loads(manifest.read_text())
    for p in (derived,index):
        metadata['artifact_identities'][p.relative_to(manifest.parent).as_posix()].update(
            sha256=hashlib.sha256(p.read_bytes()).hexdigest(),bytes=p.stat().st_size)
    manifest.write_text(json.dumps(metadata))
    case.context=_context(genome,case.root)
    row=_external(derived,track='sv',candidates=[_candidate('b1','BND')],bnd_bundle_manifest=str(manifest))
    if damage=='missing_manifest':row.pop('bnd_bundle_manifest')
    elif damage=='ledger':(manifest.parent/'records.jsonl').write_text('changed fixture')
    elif damage=='site_detach':row['candidates']=[dict(review_id='SITE',variant_type='DEL',
        site=dict(contig='chr1',position=100,ref='A',alt='C'))]
    before={str(p):p.read_bytes() for p in (derived,index,manifest)}
    if damage:
        with pytest.raises(InputValidationError):_resolve(case,[row])
    else:
        result=_resolve(case,[row]);target=result.targets[0]
        binding=result.source_bindings[0]['selected_records'][0]
        assert target.bnd_bundle_manifest==manifest
        assert binding['bundle_binding']['record_round_trip']['status']=='PASS'
        assert binding['source_role']=='EXPLICIT_VERIFIED_BND_DERIVED_WITH_NATIVE_LEDGER'
        assert len(target.loci)==2
    assert before=={str(p):p.read_bytes() for p in (derived,index,manifest)}
