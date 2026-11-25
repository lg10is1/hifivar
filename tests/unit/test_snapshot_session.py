"""Saved-session identity contracts, not screenshot interpretation or accuracy."""
from pathlib import Path
from dataclasses import replace
import os
import xml.etree.ElementTree as ET

import pytest

from hifivar.exceptions import InputValidationError, OutputValidationError
from hifivar.snapshot_session import validate_igv_session


def make_session(tmp_path,*,reference=None,locus='chr1:90-110',resources=None):
    reference=reference or tmp_path/'ref.fa'
    bam=tmp_path/'sample.bam';vcf=tmp_path/'source.vcf'
    for path in (reference,bam,vcf):
        if not path.is_absolute():path=tmp_path/path
        if not path.exists():path.write_bytes(b'synthetic input')
    session=tmp_path/'image.session.xml'
    root=ET.Element('Session',genome=str(reference),locus=locus,version='8')
    node=ET.SubElement(root,'Resources')
    for path in resources if resources is not None else (bam,vcf):
        ET.SubElement(node,'Resource',path=str(path))
    session.write_bytes(ET.tostring(root,encoding='utf-8',xml_declaration=True))
    return session,dict(reference=reference,alignment=bam,source_vcf=vcf,contig='chr1',start=90,end=110)


def review_target(tmp_path,variant_class,**overrides):
    from test_review import target
    item=target(tmp_path,variant_class,**overrides)
    safe=tmp_path/'reference.fa'
    safe.write_bytes(item.reference_fasta.read_bytes())
    Path(str(safe)+'.fai').write_bytes(Path(str(item.reference_fasta)+'.fai').read_bytes())
    return replace(item,reference_fasta=safe)


def test_exact_reference_resources_locus_returns_explicit_technical_status(tmp_path):
    session,expected=make_session(tmp_path)
    result=validate_igv_session(session,**expected)
    assert result['status']=='REFERENCE_RESOURCES_LOCUS_MATCH'
    assert result['reference_path_match'] is True and result['selected_resources_present'] is True
    assert result['human_review']=='NOT_REVIEWED' and result['variant_accuracy']=='NOT_ASSESSED'
    assert result['data_coverage']=='NOT_ASSESSED' and len(result['sha256'])==64


def test_relative_resources_resolve_against_session_parent(tmp_path):
    session,expected=make_session(tmp_path,reference=Path('ref.fa'),resources=['sample.bam','source.vcf'])
    expected['reference']=tmp_path/'ref.fa'
    assert validate_igv_session(session,**expected)['reference_path_match']


def test_local_file_uri_reference_is_supported(tmp_path):
    session,expected=make_session(tmp_path)
    root=ET.fromstring(session.read_bytes());root.set('genome',expected['reference'].as_uri())
    session.write_bytes(ET.tostring(root))
    assert validate_igv_session(session,**expected)['reference_path_match']


@pytest.mark.parametrize('locus',['All','chr1','chr2:90-110','chr1:0-110','chr1:90-9000',
                                  'chr1:89-113','chr1:110-90',''])
def test_wrong_or_whole_genome_locus_is_not_accepted(tmp_path,locus):
    session,expected=make_session(tmp_path,locus=locus)
    with pytest.raises(OutputValidationError,match='locus'):
        validate_igv_session(session,**expected)


def test_rounding_tolerance_is_small_and_recorded(tmp_path):
    session,expected=make_session(tmp_path,locus='chr1:88-112')
    result=validate_igv_session(session,**expected)
    assert result['coordinate_tolerance_bp']==2 and result['observed_locus']=='chr1:88-112'
    with pytest.raises(OutputValidationError):validate_igv_session(session,**expected,tolerance_bp=0)


def test_missing_selected_track_and_wrong_reference_fail(tmp_path):
    session,expected=make_session(tmp_path,resources=['sample.bam'])
    with pytest.raises(OutputValidationError,match='resource'):
        validate_igv_session(session,**expected)
    session,expected=make_session(tmp_path)
    expected['reference']=tmp_path/'different.fa'
    with pytest.raises(OutputValidationError,match='reference'):
        validate_igv_session(session,**expected)


@pytest.mark.parametrize('content',[b'',b'<Session',b'<Wrong/>',
    b'<!DOCTYPE Session [<!ENTITY x "expansion">]><Session genome="&x;"/>',
    b'<?xml version="1.0"?><!ENTITY x SYSTEM "file:///private"><Session/>',
    '<Session/>'.encode('utf-16')])
def test_corrupt_unsafe_or_wrong_root_xml_is_rejected(tmp_path,content):
    session,expected=make_session(tmp_path);session.write_bytes(content)
    with pytest.raises(OutputValidationError):validate_igv_session(session,**expected)


def test_large_session_is_rejected_before_read(tmp_path,monkeypatch):
    session,expected=make_session(tmp_path);session.write_bytes(b'x'*(2*1024*1024+1))
    original=Path.open
    def guarded(self,*args,**kwargs):
        if self==session:raise AssertionError('oversize session must not be read')
        return original(self,*args,**kwargs)
    monkeypatch.setattr(Path,'open',guarded)
    with pytest.raises(OutputValidationError,match='larger'):
        validate_igv_session(session,**expected)


def test_verified_plan_saves_each_session_before_its_snapshot(tmp_path):
    from hifivar.igv import IgvWrapper
    from hifivar.review import VariantClass
    from test_review import target
    item=review_target(tmp_path,VariantClass.SNV,output_directory=tmp_path/'review')
    plan=IgvWrapper().plan((item,),batch_path=item.output_directory/'batch.txt',
        snapshot_directory=item.output_directory/'screenshots',verify_session=True)
    assert plan.sessions==(plan.screenshots[0].with_suffix('.session.xml'),)
    lines=plan.batch_text.splitlines();i=lines.index('goto '+item.loci[0].igv_locus)
    assert lines[i+1]=='saveSession ./'+plan.sessions[0].name
    assert plan.to_dict()['execution_cwd']==str(plan.snapshot_directory.absolute())
    assert lines[i+2]=='snapshot '+plan.screenshots[0].name
    assert plan.to_dict()['verify_session'] is True


def test_legacy_plan_keeps_session_verification_opt_in(tmp_path):
    from hifivar.igv import IgvWrapper
    from hifivar.review import VariantClass
    from test_review import target
    item=review_target(tmp_path,VariantClass.SNV)
    plan=IgvWrapper().plan((item,),batch_path=item.output_directory/'batch.txt',
        snapshot_directory=item.output_directory/'screenshots')
    assert plan.sessions==() and 'saveSession' not in plan.batch_text


def test_verified_plan_supports_whitespace_output_path_with_owned_cwd(tmp_path):
    from hifivar.igv import IgvWrapper
    from hifivar.review import VariantClass
    from test_review import target
    item=review_target(tmp_path,VariantClass.SNV)
    plan=IgvWrapper().plan((item,),batch_path=item.output_directory/'batch.txt',
        snapshot_directory=item.output_directory/'screenshots',verify_session=True)
    assert 'saveSession ./'+plan.sessions[0].name in plan.batch_text
    assert plan.to_dict()['execution_cwd']==str(plan.snapshot_directory.absolute())


def test_session_write_alias_to_input_is_guarded_before_command(tmp_path):
    from hifivar.igv import IgvWrapper
    from hifivar.review import VariantClass
    from test_review import target,FakeIgvRunner
    item=review_target(tmp_path,VariantClass.SNV,output_directory=tmp_path/'review')
    runner=FakeIgvRunner();wrapper=IgvWrapper(runner=runner)
    plan=wrapper.plan((item,),batch_path=item.output_directory/'batch.txt',
        snapshot_directory=item.output_directory/'screenshots',verify_session=True)
    plan.sessions[0].parent.mkdir(parents=True)
    os.link(item.source_vcf,plan.sessions[0]);before=item.source_vcf.read_bytes()
    with pytest.raises(OutputValidationError,match='immutable input'):
        wrapper.run(plan,dry_run=True,overwrite=True)
    assert not runner.calls and item.source_vcf.read_bytes()==before


def test_session_and_stdout_collision_is_guarded(tmp_path):
    from hifivar.igv import IgvWrapper
    from hifivar.review import VariantClass
    from test_review import target,FakeIgvRunner
    item=review_target(tmp_path,VariantClass.SNV,output_directory=tmp_path/'review')
    runner=FakeIgvRunner();wrapper=IgvWrapper(runner=runner)
    plan=wrapper.plan((item,),batch_path=item.output_directory/'batch.txt',
        snapshot_directory=item.output_directory/'screenshots',verify_session=True)
    with pytest.raises(OutputValidationError,match='aliases collide'):
        wrapper.run(plan,dry_run=True,stdout_path=plan.sessions[0])
    assert not runner.calls


def test_wrapper_validates_saved_native_session_and_exposes_metadata(tmp_path,monkeypatch):
    from hifivar.igv import IgvWrapper
    from hifivar.review import VariantClass
    from test_review import FakeIgvRunner
    # Isolate saved-session orchestration from the separately tested image decoder.
    monkeypatch.setattr('hifivar.igv.require_png_decoder',lambda:None)
    monkeypatch.setattr('hifivar.igv.validate_png_snapshots',lambda paths:())
    item=review_target(tmp_path,VariantClass.SNV,output_directory=tmp_path/'review')
    class SessionRunner(FakeIgvRunner):
        def run(self,command,**kwargs):
            result=super().run(command,**kwargs)
            if '--batch' in command:
                assert kwargs['cwd']==(item.output_directory/'screenshots').absolute()
                for locus,screenshot in zip(item.loci,item.screenshot_paths,strict=True):
                    root=ET.Element('Session',genome=str(item.reference_fasta),locus=locus.igv_locus)
                    resources=ET.SubElement(root,'Resources')
                    for path in (item.alignment_path,item.source_vcf):ET.SubElement(resources,'Resource',path=str(path))
                    screenshot.with_suffix('.session.xml').write_bytes(ET.tostring(root))
            return result
    wrapper=IgvWrapper(runner=SessionRunner())
    plan=wrapper.plan((item,),batch_path=item.output_directory/'batch.txt',
        snapshot_directory=item.output_directory/'screenshots',verify_session=True)
    result=wrapper.run(plan)
    metadata=result.to_dict()['session_validation']
    assert metadata['enabled'] is True and metadata['expected_count']==metadata['validated_count']==1
    assert metadata['sessions'][0]['observed_locus']==item.loci[0].igv_locus
    assert metadata['sessions'][0]['variant_accuracy']=='NOT_ASSESSED'


def test_wrapper_rejects_png_success_without_required_session(tmp_path,monkeypatch):
    from hifivar.igv import IgvWrapper
    from hifivar.review import VariantClass
    from test_review import FakeIgvRunner
    monkeypatch.setattr('hifivar.igv.require_png_decoder',lambda:None)
    monkeypatch.setattr('hifivar.igv.validate_png_snapshots',lambda paths:())
    item=review_target(tmp_path,VariantClass.SNV,output_directory=tmp_path/'review')
    wrapper=IgvWrapper(runner=FakeIgvRunner())
    plan=wrapper.plan((item,),batch_path=item.output_directory/'batch.txt',
        snapshot_directory=item.output_directory/'screenshots',verify_session=True)
    with pytest.raises(OutputValidationError,match='session is missing'):
        wrapper.run(plan)
    assert plan.screenshots[0].exists()  # Image existence is insufficient.
