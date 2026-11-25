"""Synthetic BND endpoint contracts, not biological accuracy or human review."""
import csv
import gzip
import hashlib
import json
from pathlib import Path

import pytest
import yaml

from hifivar.annotation import (AnnotationInput, RegionCategory, RegionDatabase,
    VariantCategory, annotate_region_overlaps, read_annotation_inputs, read_selected_variant_loci)
from hifivar.bnd_records import derive_record_stream
from hifivar.breakends import BndPolicy, POLICY, POLICY_SHA256
from hifivar.exceptions import InputValidationError, OutputValidationError
from hifivar.igv import IgvWrapper
from hifivar.junction_consumers import read_selected_junction_records, validate_bnd_consumer_bundle
from hifivar.reference import ReferenceGenome
from hifivar.review import EvidenceClass, ReviewManifest, ReviewResult, ReviewTarget, VariantClass
from hifivar.validation import compute_sha256


def case(tmp_path, records=None):
    tmp_path.mkdir(parents=True, exist_ok=True)
    ref = tmp_path / 'reference.fa'
    ref.write_text('>chr1\n' + 'A'*1000 + '\n>chr2\n' + 'C'*1000 + '\n', encoding='utf-8')
    Path(str(ref)+'.fai').write_text('chr1\t1000\t6\t1000\t1001\nchr2\t1000\t1013\t1000\t1001\n', encoding='utf-8')
    genome = ReferenceGenome.from_fasta(ref, build='GRCh38')
    vcf = tmp_path / 'native.vcf'
    header = ('##fileformat=VCFv4.4\n##reference=GRCh38\n##contig=<ID=chr1,length=1000>\n##contig=<ID=chr2,length=1000>\n'
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">\n'
        '##INFO=<ID=END,Number=1,Type=Integer,Description="Native value">\n'
        '##INFO=<ID=SVLEN,Number=A,Type=Integer,Description="Length">\n'
        '##INFO=<ID=MATEID,Number=.,Type=String,Description="Mate IDs">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n')
    vcf.write_text(header + ''.join(records or [row()]), encoding='utf-8', newline='\n')
    bam = tmp_path / 'S1.bam'; bam.write_bytes(b'synthetic BAM; not native')
    Path(str(bam)+'.bai').write_bytes(b'synthetic index')
    return genome, vcf, bam


def row(identifier='b1', pos=100, alt='A[chr2:500[', info='SVTYPE=BND;END=900;SVLEN=800;MATEID=missing'):
    return f'chr1\t{pos}\t{identifier}\t'+('N' if pos in (0,1001) else 'A')+f'\t{alt}\t.\tPASS\t{info}\tGT\t1|0\n'


def target(tmp_path, genome, vcf, bam, *, variant_id='b1', pos=100, **kwargs):
    return ReviewTarget('R1','S1',variant_id,VariantClass.BND,'chr1',pos,pos,vcf,
        'explicit',EvidenceClass.EXPLICIT,bam,genome.fasta,tmp_path/'out',flank_bp=10,
        reference_build='GRCh38',**kwargs)


@pytest.mark.parametrize('alt,sides', [('A[chr2:500[',('LEFT','RIGHT')),('A]chr2:500]',('LEFT','LEFT')),
    ('[chr2:500[A',('RIGHT','RIGHT')),(']chr2:500]A',('RIGHT','LEFT'))])
def test_four_orientations_source_identity_and_no_span(tmp_path,alt,sides):
    genome,vcf,_=case(tmp_path,[row(alt=alt)])
    before=vcf.read_bytes()
    record=read_selected_junction_records(vcf,('b1',),reference=genome,sample_id='S1')[0]
    assert record.end==record.position==100
    assert (record.junction.local_side,record.junction.mate_side)==sides
    assert record.ref=='A' and record.alt==alt and record.ordinal==1
    assert record.record_sha256==hashlib.sha256(before.splitlines(keepends=True)[-1]).hexdigest()
    assert record.mate_qc['links'][0]['id_resolution']=='UNRESOLVED'
    assert record.to_dict()['phase_accuracy']=='NOT_ASSESSED' and vcf.read_bytes()==before


def test_annotation_queries_two_anchors_never_native_end_span(tmp_path):
    genome,vcf,_=case(tmp_path)
    item=AnnotationInput('S1',VariantCategory.SV,vcf,'caller',genome,('b1',))
    loci=read_selected_variant_loci(item)
    bed=tmp_path/'features.bed'
    bed.write_text('chr1\t99\t100\tLOCAL\nchr1\t100\t800\tNOT_SPAN\nchr2\t499\t500\tREMOTE\nchr2\t500\t501\tAFTER\n')
    result=annotate_region_overlaps(loci,(RegionDatabase(RegionCategory.GENE,bed,'1','GRCh38'),),
        reference=genome,output_tsv=tmp_path/'overlap.tsv')
    assert {x.feature_id for x in result.overlaps}=={'LOCAL','REMOTE'}
    assert {(x.to_dict()['contig'],x.to_dict()['start']) for x in result.overlaps}=={('chr1',100),('chr2',500)}
    assert {x.endpoint for x in result.overlaps}=={'LOCAL','ALT_MATE'}
    assert all(x.variant.end==100 for x in result.overlaps)
    assert [x['status'] for x in result.junction_annotations[0]['endpoints']]==['OVERLAP','OVERLAP']
    assert result.to_dict()['junction_annotation_is_fusion_or_confidence'] is False


@pytest.mark.parametrize('pos,alt,expected_windows,skipped',[(100,'A.', ['chr1:90-110'],['UNKNOWN_SINGLE_BREAKEND']),
    (100,'.A',['chr1:90-110'],['UNKNOWN_SINGLE_BREAKEND']),
    (0,'.[chr2:500[',['chr2:490-510'],['NO_REFERENCE_BASE_VIRTUAL_ANCHOR']),
    (100,'A[chr2:1001[',['chr1:90-110'],['NO_REFERENCE_BASE_VIRTUAL_ANCHOR']),
    (1001,']chr2:0].',[],['NO_REFERENCE_BASE_VIRTUAL_ANCHOR','NO_REFERENCE_BASE_VIRTUAL_ANCHOR'])])
def test_review_single_virtual_and_two_anchor_policy(tmp_path,pos,alt,expected_windows,skipped):
    genome,vcf,bam=case(tmp_path,[row(pos=pos,alt=alt,info='SVTYPE=BND')])
    item=target(tmp_path,genome,vcf,bam,pos=pos)
    assert [x.igv_locus for x in item.loci]==expected_windows
    assert [x['status'] for x in item.unrendered_bnd_endpoints]==skipped
    plan=IgvWrapper().plan((item,),batch_path=tmp_path/'batch',snapshot_directory=tmp_path/'out/screenshots')
    assert plan.batch_text.count('goto ')==len(expected_windows)
    assert len(item.screenshot_paths)==len(expected_windows) and item.start==pos
    assert item.to_dict()['junction_record']['junction']['position']==pos


def test_ordinary_paired_review_mate_in_alt_is_not_resolved_mateid(tmp_path):
    genome,vcf,bam=case(tmp_path)
    item=target(tmp_path,genome,vcf,bam)
    assert [x.igv_locus for x in item.loci]==['chr1:90-110','chr2:490-510']
    assert item.to_dict()['junction_record']['mate_qc']['links'][0]['id_resolution']=='UNRESOLVED'
    with pytest.raises(InputValidationError,match='mate differs'):
        target(tmp_path,genome,vcf,bam,mate_contig='chr2',mate_position=501)
    with pytest.raises(InputValidationError,match='exact selected VCF anchor'):
        ReviewTarget('R','S1','b1',VariantClass.BND,'chr1',100,900,vcf,'caller',EvidenceClass.EXPLICIT,bam,genome.fasta,tmp_path/'bad')


def test_virtual_and_single_annotation_statuses_not_lost(tmp_path):
    genome,vcf,_=case(tmp_path,[row('virtual',pos=0,alt='.[chr2:500[',info='SVTYPE=BND'),row('single',alt='A.',info='SVTYPE=BND')])
    item=AnnotationInput('S1',VariantCategory.SV,vcf,'caller',genome,('virtual','single'))
    bed=tmp_path/'empty_hits.bed';bed.write_text('chr1\t200\t201\tNO_HIT\n')
    result=annotate_region_overlaps(read_selected_variant_loci(item),(RegionDatabase(RegionCategory.REPEAT,bed,'1','GRCh38'),),reference=genome,output_tsv=tmp_path/'result.tsv')
    assert not result.overlaps
    assert len(result.junction_annotations)==2
    assert result.junction_annotations[0]['endpoints'][0]['position']==0
    assert result.junction_annotations[0]['endpoints'][0]['status']=='NO_REFERENCE_BASE_VIRTUAL_ANCHOR'
    assert result.junction_annotations[1]['endpoints'][1]['status']=='UNKNOWN_SINGLE_BREAKEND'


@pytest.mark.parametrize('records,selected',[
    ([row(),row()],('b1',)),([row()],('absent',)),([row(alt='A[2:500[')],('b1',)),
    ([row(alt='A[chr2:1002[')],('b1',)),([row(alt='A[chr2:500[,T')],('b1',)),
    ([row(info='SVTYPE=DEL')],('b1',)),([row(info='SVTYPE=BND;MATEID')],('b1',)),
    ([row()],('.',)),([row()],('b1','b1'))])
def test_invalid_selected_sources_fail_closed(tmp_path,records,selected):
    genome,vcf,_=case(tmp_path,records)
    with pytest.raises(InputValidationError):
        read_selected_junction_records(vcf,selected,reference=genome,sample_id='S1')


def test_duplicate_mate_ids_are_qc_not_guessed_or_dropped(tmp_path):
    records=[row(info='SVTYPE=BND;MATEID=other'),row('other',pos=500),row('other',pos=600)]
    genome,vcf,_=case(tmp_path,records)
    record=read_selected_junction_records(vcf,('b1',),reference=genome,sample_id='S1')[0]
    assert record.mate_qc['links']==[{'mate_id':'other','matching_record_count':2,'id_resolution':'AMBIGUOUS_DUPLICATE_ID'}]


def make_bundle(tmp_path):
    genome,native,bam=case(tmp_path/'inputs')
    bundle=tmp_path/'bundle';bundle.mkdir();work=bundle/'work';work.mkdir();copies=bundle/'native';copies.mkdir()
    for name,data in [('native.vcf',native.read_bytes()),('native.vcf.tbi',b'synthetic-index'),('raw.vcf',native.read_bytes()),('evidence.tsv',b'synthetic-evidence')]:
        (copies/name).write_bytes(data)
    policy=BndPolicy(POLICY['policy_id'],1,POLICY_SHA256,POLICY['mode'])
    plain=bundle/'work/derived.vcf'
    derive_record_stream(native,plain,bundle/'records.jsonl',bundle/'junctions.tsv',reference_lengths={'chr1':1000,'chr2':1000},sample_id='S1',native_sha256=compute_sha256(native),policy=policy,work_directory=work)
    derived=bundle/'S1.bndrepr-v1.sv.vcf.gz';derived.write_bytes(gzip.compress(plain.read_bytes()))
    Path(str(derived)+'.tbi').write_bytes(b'synthetic-index; native tests are separate')
    (bundle/'native_result.json').write_text(json.dumps({'fixture':'synthetic representation; not native calling'}))
    def identity(p):return {'path':str(p.resolve()),'sha256':compute_sha256(p),'bytes':p.stat().st_size}
    before={str(p.resolve()):identity(p) for p in (genome.fasta,genome.fai)}
    artifacts=[derived,Path(str(derived)+'.tbi'),bundle/'records.jsonl',bundle/'junctions.tsv',bundle/'native_result.json',*copies.iterdir()]
    data={'schema_version':1,'status':'completed','policy':policy.to_dict(),'native_output_role_unchanged':True,
        'native_binary_copies':{'vcf':'native/native.vcf','index':'native/native.vcf.tbi','raw':'native/raw.vcf','evidence':'native/evidence.tsv'},
        'native_inputs_before':before,'native_inputs_after':before,'artifact_identities':{p.relative_to(bundle).as_posix():identity(p) for p in artifacts}}
    manifest=bundle/'manifest.json';manifest.write_text(json.dumps(data))
    return genome,derived,bam,manifest


def test_explicit_derived_bundle_roundtrip_and_annotation_selection(tmp_path):
    genome,derived,bam,manifest=make_bundle(tmp_path)
    before=compute_sha256(derived)
    with pytest.raises(InputValidationError,match='explicit verified bundle'):
        read_selected_junction_records(derived,('b1',),reference=genome,sample_id='S1')
    item=AnnotationInput('S1',VariantCategory.SV,derived,'jasmine',genome,('b1',),manifest)
    record=read_selected_variant_loci(item)[0].junction_record
    assert record.source_role=='EXPLICIT_VERIFIED_BND_DERIVED_WITH_NATIVE_LEDGER'
    assert record.bundle_binding['record_round_trip']['status']=='PASS'
    review=target(tmp_path,genome,derived,bam,bnd_bundle_manifest=manifest)
    assert len(review.loci)==2 and compute_sha256(derived)==before
    plan=tmp_path/'annotation.tsv'
    plan.write_text('sample\tvariant_category\tsource_vcf\tsource_tool\tsource_variant_ids\tbnd_bundle_manifest\n'+f'S1\tsv\t{derived}\tjasmine\tb1\t{manifest}\n')
    assert read_annotation_inputs(plan,reference=genome)[0].bnd_bundle_manifest==manifest


@pytest.mark.parametrize('damage',['ledger','manifest_state','missing_role','reference','wrong_vcf'])
def test_derived_bundle_damage_cannot_be_hidden(tmp_path,damage):
    genome,derived,_,manifest=make_bundle(tmp_path)
    data=json.loads(manifest.read_text())
    source=derived
    if damage=='ledger':(manifest.parent/'records.jsonl').write_text('changed')
    elif damage=='manifest_state':data['status']='failed';manifest.write_text(json.dumps(data))
    elif damage=='missing_role':data['artifact_identities'].pop('junctions.tsv');manifest.write_text(json.dumps(data))
    elif damage=='reference':genome.fasta.write_text('>chr1\nCHANGED\n')
    else:source=manifest.parent/'copy.vcf.gz';source.write_bytes(derived.read_bytes());Path(str(source)+'.tbi').write_bytes(b'tbi')
    with pytest.raises(InputValidationError):validate_bnd_consumer_bundle(manifest,source,genome,'S1')


@pytest.mark.parametrize('which',['vcf','reference','fai','bed'])
def test_overlap_cannot_replace_source_even_with_overwrite(tmp_path,which):
    genome,vcf,_=case(tmp_path)
    loci=read_selected_variant_loci(AnnotationInput('S1',VariantCategory.SV,vcf,'caller',genome,('b1',)))
    bed=tmp_path/'features.bed';bed.write_text('chr1\t99\t100\tLOCAL\n')
    dest={'vcf':vcf,'reference':genome.fasta,'fai':genome.fai,'bed':bed}[which]
    before=dest.read_bytes()
    with pytest.raises(OutputValidationError,match='immutable'):
        annotate_region_overlaps(loci,(RegionDatabase(RegionCategory.GENE,bed,'1','GRCh38'),),reference=genome,output_tsv=dest,overwrite=True)
    assert dest.read_bytes()==before


def test_late_source_change_cannot_yield_annotation_or_igv_plan(tmp_path):
    genome,vcf,bam=case(tmp_path)
    item=AnnotationInput('S1',VariantCategory.SV,vcf,'caller',genome,('b1',))
    loci=read_selected_variant_loci(item)
    review=target(tmp_path,genome,vcf,bam)
    vcf.write_text(vcf.read_text().replace('1|0','0|1'))
    with pytest.raises(InputValidationError,match='bytes changed'):
        annotate_region_overlaps(loci,(),reference=genome,output_tsv=tmp_path/'new.tsv')
    assert not (tmp_path/'new.tsv').exists()
    with pytest.raises(InputValidationError,match='bytes changed'):
        IgvWrapper().plan((review,),batch_path=tmp_path/'batch',snapshot_directory=tmp_path/'out/screenshots')


def test_ordinary_sv_and_tr_records_retain_existing_span_behavior(tmp_path):
    genome,vcf,_=case(tmp_path,[row('del',alt='<DEL>',info='SVTYPE=DEL;END=150')])
    for category in (VariantCategory.SV,VariantCategory.TR):
        locus=read_selected_variant_loci(AnnotationInput('S1',category,vcf,'caller',genome,('del',)))[0]
        assert (locus.start,locus.end,locus.junction_record)==(100,150,None)


def test_explicit_bundle_source_roles_in_receipt_manifest(tmp_path):
    from hifivar.attempt_receipts import _manifest_sources, _source_paths
    genome,derived,_,manifest=make_bundle(tmp_path)
    table=tmp_path/'selected.tsv'
    table.write_text('sample\tvariant_category\tsource_vcf\tsource_tool\tsource_variant_ids\tbnd_bundle_manifest\n'+f'S1\tsv\t{derived}\tjasmine\tb1\t{manifest}\n')
    actual=_source_paths({'S1':{'vcf':derived,'index':Path(str(derived)+'.tbi'),'bnd_bundle_manifest':manifest}},('S1',))
    expected=_manifest_sources(table,'annotation.overlap',('S1',),'sv')
    assert actual==expected and actual['S1']['bnd_bundle_manifest']==manifest


def test_isolated_derived_annotation_is_rejected_before_backend_execution(tmp_path):
    genome,derived,_,manifest=make_bundle(tmp_path)
    with pytest.raises(InputValidationError,match='explicit bundle manifest'):
        AnnotationInput('S1',VariantCategory.SV,derived,'jasmine',genome,('b1',))
    assert AnnotationInput('S1',VariantCategory.SV,derived,'jasmine',genome,('b1',),manifest).source_identity['bnd_bundle_binding']['role']=='EXPLICIT_VERIFIED_BND_DERIVED_WITH_NATIVE_LEDGER'


@pytest.mark.parametrize('category,alt,info',[(VariantCategory.SV,'<DEL>,<DUP>','SVTYPE=DEL;END=150'),(VariantCategory.TR,'ACAC,ACACAC','TRID=repeat1;END=150')])
def test_non_bnd_multiallelic_source_is_not_sent_to_single_alt_bnd_parser(tmp_path,category,alt,info):
    genome,vcf,_=case(tmp_path,[row('multi',alt=alt,info=info)])
    before=vcf.read_bytes()
    locus=read_selected_variant_loci(AnnotationInput('S1',category,vcf,'caller',genome,('multi',)))[0]
    assert (locus.start,locus.end,locus.junction_record)==(100,150,None)
    assert vcf.read_bytes()==before
