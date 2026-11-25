"""Direct snapshot entry contracts, using synthetic sources and fake native tools."""
from dataclasses import replace
import copy
import csv
import json
from pathlib import Path
from types import SimpleNamespace
import xml.etree.ElementTree as ET

import pytest
import yaml

from hifivar.exceptions import CommandExecutionError, InputValidationError, OutputValidationError, ReferenceError
from hifivar.reference import ReferenceGenome
from hifivar.review_snapshots import run_review_snapshots
from hifivar.review_cli import import_review_decisions
from test_review import FakeIgvRunner as BasicIgvRunner, PNG
from test_review_sources import case, _vcf, _record, _external, _plan, _candidate


class FakeIgvRunner(BasicIgvRunner):
    """Synthetic native-session fixture; never substitute for native acceptance."""
    def __init__(self, *, session_damage=None):
        super().__init__()
        self.session_damage=session_damage

    def run(self, argv, **kwargs):
        result=super().run(argv,**kwargs)
        if '--batch' not in argv or kwargs.get('dry_run'):
            return result
        genome=(str(argv[argv.index('--genome')+1]) if '--genome' in argv else None)
        resources=[];locus='All'
        for line in Path(argv[-1]).read_text(encoding='utf-8').splitlines():
            if line=='new':resources=[];locus='All'
            # IGV 2.19.8 forwards literal genome quotes; unlike load, do not strip them.
            elif line.startswith('genome '):genome=line.split(' ',1)[1]
            elif line.startswith('load '):resources.append(line.split(' ',1)[1].strip('"'))
            elif line.startswith('goto '):locus=line.split(' ',1)[1]
            elif line.startswith('saveSession '):
                path=Path(line.split(' ',1)[1])
                if not path.is_absolute():path=Path(kwargs.get('cwd',Path.cwd()))/path
                root=ET.Element('Session',genome=genome if self.session_damage!='reference' else str(path.parent/'other.fa'),
                                locus=locus if self.session_damage!='locus' else 'All',version='8')
                node=ET.SubElement(root,'Resources')
                for resource in (resources[:-1] if self.session_damage=='resources' else resources):
                    ET.SubElement(node,'Resource',path=resource)
                path.parent.mkdir(parents=True,exist_ok=True)
                ET.ElementTree(root).write(path,encoding='utf-8',xml_declaration=True)
        return result


@pytest.fixture
def snapshot_case(case, monkeypatch):
    """Bind real source resolver to synthetic context and an explicit fake probe."""
    import hifivar.review_snapshots as module
    user = case.root/'user.yaml'
    user.write_text('review: synthetic-test\n', encoding='utf-8')
    case.context.source_sample_sheet = None
    case.config['review']['igv_executable'] = 'synthetic-igv'
    case.config['review']['selection_file'] = None
    case.config['workflow']['samtools_executable'] = 'synthetic-samtools'
    case.probes = []
    case.probe_failure = False
    monkeypatch.setattr(module.AnalysisContext, 'from_config', lambda config: case.context)
    class Probe:
        def __init__(self, runner, executable):
            self.executable = executable
        def run(self, path, index, sample, reference, logdir):
            case.probes.append((path,index,sample,reference,logdir,self.executable))
            logdir.mkdir(parents=True, exist_ok=True)
            (logdir/'synthetic-preflight.txt').write_text('fixture native probe invoked\n', encoding='utf-8')
            if case.probe_failure:
                raise CommandExecutionError('synthetic samtools failure')
            return {'sample_id':sample, 'header_identity':'PASS', 'quickcheck':'PASS',
                    'indexed_access':'PASS', 'scope':'SYNTHETIC_NOT_NATIVE_BAM_QC'}
    monkeypatch.setattr(module, 'IndexedBamPreflight', Probe)
    case.cli_config = SimpleNamespace(sources={'user':user},to_dict=lambda:copy.deepcopy(case.config))
    case.output = case.root/'snapshots'
    return case


def bind(case, rows=None, *, records=None):
    vcf = _vcf(case, records or [_record()])
    plan = _plan(case, rows or [_external(vcf)])
    case.config['review']['source_plan'] = str(plan)
    return vcf, plan


def test_snapshot_dry_run_resolves_exact_sources_without_command_or_output(snapshot_case):
    case = snapshot_case
    vcf, plan = bind(case)
    runner = FakeIgvRunner()
    before = {p:p.read_bytes() for p in (vcf,Path(str(vcf)+'.tbi'),plan,case.context.reference.fasta)}
    result = run_review_snapshots(case.cli_config,case.output,dry_run=True,runner=runner)
    assert result['status'] == 'PLANNED' and result['screenshots_generated'] is False
    assert result['target_count'] == result['expected_screenshots'] == 1
    assert result['source_bindings'][0]['selected_records'][0]['record_id'] == 'v1'
    assert not runner.calls and not case.probes and not case.output.exists()
    assert before == {p:p.read_bytes() for p in before}


def test_generated_sources_cannot_be_replayed_by_direct_entry(snapshot_case):
    case = snapshot_case
    plan = _plan(case,[dict(sample_id='S1',track='small',origin='GENERATED',
        producer='small.deepvariant',output_role='vcf',candidates=[_candidate()])])
    case.config['review']['source_plan'] = str(plan)
    runner = FakeIgvRunner()
    with pytest.raises(InputValidationError,match='EXTERNAL'):
        run_review_snapshots(case.cli_config,case.output,runner=runner)
    assert not runner.calls and not case.probes and not case.output.exists()


@pytest.mark.parametrize('damage',['index','sample','missing_record','build'])
def test_invalid_source_does_not_reserve_output_or_execute(snapshot_case,damage):
    case = snapshot_case
    vcf = _vcf(case,[_record()],sample='OTHER' if damage=='sample' else 'S1',
               build='GRCh37' if damage=='build' else 'GRCh38')
    row = _external(vcf)
    if damage=='index':Path(str(vcf)+'.tbi').write_bytes(b'not an index')
    if damage=='missing_record':row['candidates'][0]['variant_id']='absent'
    plan = _plan(case,[row])
    case.config['review']['source_plan'] = str(plan)
    runner = FakeIgvRunner()
    with pytest.raises((InputValidationError,OutputValidationError,ReferenceError)):
        run_review_snapshots(case.cli_config,case.output,runner=runner)
    assert not runner.calls and not case.probes and not case.output.exists()


def test_completed_technical_bundle_keeps_decisions_pending_and_sources_immutable(snapshot_case):
    case = snapshot_case
    vcf, plan = bind(case)
    watched = (vcf,Path(str(vcf)+'.tbi'),plan,case.context.reference.fasta,
               case.context.reference.fai,case.context.samples[0].sample.input.files[0])
    before = {p:p.read_bytes() for p in watched}
    runner = FakeIgvRunner()
    result = run_review_snapshots(case.cli_config,case.output,runner=runner)
    assert result['status'] == 'SNAPSHOTS_GENERATED_TECHNICAL_QA'
    assert result['manual_review_status'] == 'NOT_REVIEWED' and result['variant_accuracy'] == 'NOT_ASSESSED'
    assert result['igv']['snapshot_validation']['validated_count'] == 1
    assert result['igv']['session_validation']['enabled'] is True
    assert result['igv']['session_validation']['validated_count'] == 1
    assert len(case.probes)==1 and case.probes[0][2]=='S1'
    assert (case.output/'screenshots/R1.locus01.png').read_bytes()==PNG
    assert 'NOT_REVIEWED' in (case.output/'index.html').read_text(encoding='utf-8')
    with (case.output/'manual_decisions.tsv').open(encoding='utf-8',newline='') as stream:
        decisions=list(csv.DictReader(stream,delimiter='\t'))
    assert decisions[0]['status']=='NOT_REVIEWED' and not decisions[0]['reviewer']
    manifest=json.loads((case.output/'review_manifest.json').read_text(encoding='utf-8'))
    assert manifest['results'][0]['status']=='NOT_REVIEWED'
    assert json.loads((case.output/'snapshot_result.json').read_text(encoding='utf-8'))['screenshots_generated'] is True
    assert not (case.output/'snapshot_failure.json').exists()
    assert before=={p:p.read_bytes() for p in before}


@pytest.mark.parametrize('resource',['vcf','bam','plan','settings'])
def test_changed_input_preserves_failure_but_never_completion_receipt(snapshot_case,resource):
    case = snapshot_case
    vcf, plan = bind(case)
    settings=case.root/'settings.yaml'
    settings.write_text('sleep_interval_ms: 0\n',encoding='utf-8')
    selected={'vcf':vcf,'bam':case.context.samples[0].sample.input.files[0],
              'plan':plan,'settings':settings}[resource]
    class Changed(FakeIgvRunner):
        def run(self,argv,**kwargs):
            result=super().run(argv,**kwargs)
            if '--batch' in argv:selected.write_bytes(selected.read_bytes()+b'\nchanged fixture')
            return result
    with pytest.raises((InputValidationError,OutputValidationError),match='changed'):
        run_review_snapshots(case.cli_config,case.output,settings_file=settings,runner=Changed())
    assert not (case.output/'snapshot_result.json').exists()
    assert (case.output/'snapshot_failure.json').is_file()
    assert (case.output/'screenshots/R1.locus01.png').is_file()


def test_native_preflight_failure_is_preserved_and_igv_is_not_started(snapshot_case):
    case = snapshot_case
    bind(case)
    case.probe_failure=True
    runner=FakeIgvRunner()
    with pytest.raises(CommandExecutionError,match='samtools failure'):
        run_review_snapshots(case.cli_config,case.output,runner=runner)
    assert not runner.calls and len(case.probes)==1
    assert (case.output/'logs/synthetic-preflight.txt').is_file()
    failure=json.loads((case.output/'snapshot_failure.json').read_text(encoding='utf-8'))
    assert failure['status']=='FAILED' and failure['manual_review_status']=='NOT_REVIEWED'
    assert not (case.output/'snapshot_result.json').exists()


def test_corrupt_png_does_not_become_completed_bundle(snapshot_case):
    case=snapshot_case
    bind(case)
    class BadImage(FakeIgvRunner):
        def run(self,argv,**kwargs):
            result=super().run(argv,**kwargs)
            if '--batch' in argv:(case.output/'screenshots/R1.locus01.png').write_bytes(PNG[:24])
            return result
    with pytest.raises(OutputValidationError):
        run_review_snapshots(case.cli_config,case.output,runner=BadImage())
    assert (case.output/'snapshot_failure.json').is_file() and not (case.output/'snapshot_result.json').exists()


@pytest.mark.parametrize('damage',['reference','locus','resources'])
def test_decodable_png_with_wrong_native_session_never_completes(snapshot_case,damage):
    case=snapshot_case
    bind(case)
    with pytest.raises(OutputValidationError,match='session'):
        run_review_snapshots(case.cli_config,case.output,runner=FakeIgvRunner(session_damage=damage))
    assert (case.output/'screenshots/R1.locus01.png').read_bytes()==PNG
    assert (case.output/'snapshot_failure.json').is_file()
    assert not (case.output/'snapshot_result.json').exists()
    assert not (case.output/'review_bundle.json').exists()


def test_scientific_snv_window_is_clipped_to_contig(snapshot_case):
    case=snapshot_case
    bind(case,records=[_record(pos=1)])
    result=run_review_snapshots(case.cli_config,case.output,runner=FakeIgvRunner())
    assert result['expected_screenshots']==1
    batch=(case.output/'review.igv.batch').read_text(encoding='utf-8')
    assert 'goto chr1:1-41' in batch
    manifest=json.loads((case.output/'review_manifest.json').read_text(encoding='utf-8'))
    assert manifest['results'][0]['locus']==['chr1:1-41']


def test_whitespace_reference_fails_before_native_execution_and_output(snapshot_case):
    case=snapshot_case
    original=case.context.reference.fasta
    reference=case.root/'reference with spaces.fa'
    reference.write_bytes(original.read_bytes())
    Path(str(reference)+'.fai').write_bytes(case.context.reference.fai.read_bytes())
    case.context.reference=ReferenceGenome.from_fasta(reference,build='GRCh38')
    case.output=case.root/'snapshots with spaces'
    vcf,plan=bind(case)
    watched=(reference,Path(str(reference)+'.fai'),vcf,Path(str(vcf)+'.tbi'),plan)
    before={path:path.read_bytes() for path in watched}
    runner=FakeIgvRunner()
    with pytest.raises(InputValidationError,match='IGV_REFERENCE_PATH_UNSUPPORTED'):
        run_review_snapshots(case.cli_config,case.output,runner=runner)
    assert not runner.calls and not case.probes and not case.output.exists()
    assert before=={path:path.read_bytes() for path in watched}


def test_whitespace_output_with_normal_reference_keeps_owned_session_paths(snapshot_case):
    case=snapshot_case
    case.output=case.root/'snapshots with spaces'
    bind(case)
    runner=FakeIgvRunner()
    result=run_review_snapshots(case.cli_config,case.output,runner=runner)
    assert '--genome' not in result['igv']['plan']['args']
    batch=(case.output/'review.igv.batch').read_text(encoding='utf-8')
    assert any(line=='genome '+str(case.context.reference.fasta.absolute()) for line in batch.splitlines())
    assert 'saveSession ./R1.locus01.session.xml' in batch
    assert runner.calls[-1][1]['cwd']==case.output/'screenshots'
    assert (case.output/'screenshots/R1.locus01.session.xml').is_file()
    assert result['igv']['session_validation']['validated_count']==1


def test_scientific_long_deletion_has_overview_and_two_breakpoints(snapshot_case):
    case=snapshot_case
    ref=case.context.reference.fasta
    ref.write_text('>chr1\n'+'A'*5000+'\n>chr2\n'+'C'*1000+'\n',encoding='utf-8')
    case.context.reference.fai.write_text('chr1\t5000\t6\t5000\t5001\nchr2\t1000\t5013\t1000\t1001\n',encoding='utf-8')
    case.context.reference=ReferenceGenome.from_fasta(ref,build='GRCh38')
    vcf=_vcf(case,[_record('del1',pos=1000,alt='<DEL>',info='SVTYPE=DEL;END=3000')])
    # _vcf header lengths must match the actual reference selected by this test.
    import gzip
    from test_review_sources import _bgzf
    vcf.write_bytes(_bgzf(gzip.decompress(vcf.read_bytes()).replace(b'ID=chr1,length=1000',b'ID=chr1,length=5000')))
    plan=_plan(case,[_external(vcf,track='sv',candidates=[_candidate('del1','DEL')])])
    case.config['review']['source_plan']=str(plan)
    result=run_review_snapshots(case.cli_config,case.output,runner=FakeIgvRunner())
    assert result['expected_screenshots']==3 and result['igv']['snapshot_validation']['validated_count']==3
    target=result['igv']['evidence'][0]['target']
    assert [x['label'] for x in target['loci']]==['overview','left_breakpoint','right_breakpoint']
    assert len(target['screenshots'])==3 and target['view_mode']=='scientific'
    assert 'left_breakpoint' in (case.output/'index.html').read_text(encoding='utf-8')


@pytest.mark.parametrize('alt',['A[chr2:500[','A]chr2:500]','[chr2:500[A',']chr2:500]A'])
def test_bnd_endpoints_are_not_replaced_by_end_info(snapshot_case,alt):
    case=snapshot_case
    vcf=_vcf(case,[_record('b1',alt=alt,info='SVTYPE=BND;END=900')])
    plan=_plan(case,[_external(vcf,track='sv',candidates=[_candidate('b1','BND')])])
    case.config['review']['source_plan']=str(plan)
    result=run_review_snapshots(case.cli_config,case.output,runner=FakeIgvRunner())
    target=result['igv']['evidence'][0]['target']
    assert [(x['contig'],x['variant_start']) for x in target['loci']]==[('chr1',100),('chr2',500)]
    assert [(x['window_start'],x['window_end']) for x in target['loci']]==[(1,350),(250,750)]
    assert result['expected_screenshots']==2


def test_unrenderable_bnd_cannot_claim_screenshots_generated(snapshot_case):
    case=snapshot_case
    vcf=_vcf(case,[_record('b1',pos=1001,ref='N',alt=']chr2:0].',info='SVTYPE=BND')])
    plan=_plan(case,[_external(vcf,track='sv',candidates=[_candidate('b1','BND')])])
    case.config['review']['source_plan']=str(plan)
    runner=FakeIgvRunner()
    with pytest.raises(InputValidationError,match='renderable'):
        run_review_snapshots(case.cli_config,case.output,runner=runner)
    assert not runner.calls and not case.probes and not case.output.exists()


@pytest.mark.parametrize('text',['view_mode: unknown\n','sleep_interval_ms: true\n','unknown: 1\n',
                               'max_panel_height: 1\n','view_mode: scientific\nview_mode: legacy\n'])
def test_bad_settings_fail_before_native_execution(snapshot_case,text):
    case=snapshot_case
    bind(case)
    path=case.root/'settings.yaml';path.write_text(text,encoding='utf-8')
    runner=FakeIgvRunner()
    with pytest.raises(InputValidationError):run_review_snapshots(case.cli_config,case.output,settings_file=path,runner=runner)
    assert not runner.calls and not case.probes and not case.output.exists()


def test_existing_output_is_not_overwritten(snapshot_case):
    case=snapshot_case
    bind(case)
    case.output.mkdir();old=case.output/'keep.txt';old.write_bytes(b'preserve me')
    with pytest.raises(OutputValidationError):run_review_snapshots(case.cli_config,case.output,runner=FakeIgvRunner())
    assert old.read_bytes()==b'preserve me' and not case.probes


def test_cli_has_explicit_snapshot_action_and_external_settings():
    from hifivar.cli import build_parser
    args=build_parser().parse_args(['review','snapshots','--output-dir','review_out','--settings','capture.yaml'])
    assert args.review_command=='snapshots' and args.output_dir==Path('review_out')


def tree_bytes(directory):
    return {p.relative_to(directory).as_posix():p.read_bytes() for p in directory.rglob('*') if p.is_file()}


def test_snapshot_bundle_pending_template_can_be_imported_without_inventing_review(snapshot_case):
    case=snapshot_case
    bind(case)
    run_review_snapshots(case.cli_config,case.output,runner=FakeIgvRunner())
    before=tree_bytes(case.output)
    dest=case.root/'pending_import'
    result=import_review_decisions(case.output,case.output/'manual_decisions.tsv',dest)
    assert result['status']=='DECISIONS_IMPORTED' and result['reviewed_targets']==0
    assert result['status_counts']=={'NOT_REVIEWED':1}
    manifest=json.loads((dest/'review_manifest.adjudicated.json').read_text(encoding='utf-8'))
    assert manifest['results'][0]['status']=='NOT_REVIEWED' and not manifest['results'][0]['reviewer']
    assert result['status_is_truth'] is False and result['reviewer_identity_verified'] is False
    assert tree_bytes(case.output)==before


def test_snapshot_bundle_supplied_decision_import_is_independent_and_immutable(snapshot_case):
    case=snapshot_case
    bind(case)
    run_review_snapshots(case.cli_config,case.output,runner=FakeIgvRunner())
    before=tree_bytes(case.output)
    table=case.root/'supplied_decisions.tsv'
    table.write_text('review_id\tsample\tvariant_id\tstatus\treviewer\treviewed_at\tcomment\n'
                     'R1\tS1\tv1\tUNCERTAIN\tsynthetic-reviewer\t2026-11-03T10:00:00+08:00\tambiguous fixture\n',encoding='utf-8')
    dest=case.root/'supplied_import'
    result=import_review_decisions(case.output,table,dest)
    assert result['reviewed_targets']==1 and result['status_counts']=={'UNCERTAIN':1}
    result_manifest=json.loads((dest/'review_manifest.adjudicated.json').read_text(encoding='utf-8'))
    row=result_manifest['results'][0]
    assert row['status']=='UNCERTAIN' and row['reviewer']=='synthetic-reviewer'
    assert row['status_is_truth'] is False and row['reviewer_identity_verified'] is False
    assert row['screenshots']==json.loads(before['review_manifest.json'])['results'][0]['screenshots']
    assert tree_bytes(case.output)==before


@pytest.mark.parametrize('name',['screenshots/R1.locus01.png','screenshots/R1.locus01.session.xml','snapshot_result.json'])
def test_snapshot_import_rejects_tampered_image_or_execution_record(snapshot_case,name):
    case=snapshot_case
    bind(case)
    run_review_snapshots(case.cli_config,case.output,runner=FakeIgvRunner())
    selected=case.output/name
    selected.write_bytes(selected.read_bytes()+b'tampered')
    dest=case.root/'tampered_import'
    with pytest.raises(InputValidationError,match='changed'):
        import_review_decisions(case.output,case.output/'manual_decisions.tsv',dest)
    assert not dest.exists()


def test_relocated_snapshot_bundle_fails_with_typed_error_without_retargeting(snapshot_case):
    case=snapshot_case
    bind(case)
    run_review_snapshots(case.cli_config,case.output,runner=FakeIgvRunner())
    before=tree_bytes(case.output)
    import shutil
    relocated=case.root/'copied_snapshot_bundle'
    shutil.copytree(case.output,relocated)
    dest=case.root/'copied_import'
    with pytest.raises(InputValidationError):
        import_review_decisions(relocated,relocated/'manual_decisions.tsv',dest)
    assert not dest.exists() and tree_bytes(case.output)==before
