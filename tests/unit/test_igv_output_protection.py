"""Finite write-path safety; synthetic files never establish native IGV accuracy."""
from pathlib import Path
import os

import pytest

from hifivar.exceptions import OutputValidationError, ToolVersionError
from hifivar.igv import IgvRunStatus, IgvWrapper
from hifivar.review import VariantClass
from test_review import FakeIgvRunner, target, PNG
from test_junction_consumers import make_bundle, target as bnd_target


def plan_for(item, runner=None, *, batch=None, executable='igv.sh'):
    runner = runner or FakeIgvRunner()
    wrapper = IgvWrapper(executable, runner=runner)
    plan = wrapper.plan((item,), batch_path=batch or item.output_directory/'review.igv.batch',
                        snapshot_directory=item.output_directory/'screenshots')
    return wrapper, plan, runner


def reject(wrapper, plan, runner, *, resource, **kwargs):
    before = resource.read_bytes()
    with pytest.raises(OutputValidationError, match='IGV_OUTPUT_PROTECTION'):
        wrapper.run(plan, **kwargs)
    assert resource.read_bytes() == before
    assert not runner.calls


@pytest.mark.parametrize('which', ['vcf', 'bam', 'reference', 'fai'])
@pytest.mark.parametrize('destination', ['batch', 'stdout', 'stderr'])
@pytest.mark.parametrize('overwrite', [False, True])
def test_inputs_never_become_batch_or_logs(tmp_path, which, destination, overwrite):
    item = target(tmp_path, VariantClass.SNV)
    resource = dict(vcf=item.source_vcf, bam=item.alignment_path,
                    reference=item.reference_fasta, fai=Path(str(item.reference_fasta)+'.fai'))[which]
    wrapper, plan, runner = plan_for(item, batch=resource if destination=='batch' else None)
    kwargs = {destination+'_path': resource} if destination!='batch' else {}
    reject(wrapper, plan, runner, resource=resource, overwrite=overwrite, **kwargs)
    assert not item.output_directory.exists()


@pytest.mark.parametrize('suffix,short', [('.bai',False),('.bai',True),('.csi',False),('.csi',True)])
def test_every_existing_alignment_index_alias_is_protected(tmp_path,suffix,short):
    item=target(tmp_path,VariantClass.SNV)
    index=item.alignment_path.with_suffix(suffix) if short else Path(str(item.alignment_path)+suffix)
    if not index.exists(): index.write_bytes(b'synthetic alternate index')
    wrapper,plan,runner=plan_for(item)
    reject(wrapper,plan,runner,resource=index,stderr_path=index,overwrite=True)


@pytest.mark.parametrize('suffix',['.tbi','.csi'])
def test_existing_vcf_index_alias_is_protected(tmp_path,suffix):
    item=target(tmp_path,VariantClass.SNV);index=Path(str(item.source_vcf)+suffix)
    index.write_bytes(b'synthetic VCF index')
    wrapper,plan,runner=plan_for(item)
    reject(wrapper,plan,runner,resource=index,stdout_path=index,overwrite=True)


def linked(source,destination,kind):
    destination.parent.mkdir(parents=True,exist_ok=True)
    if kind=='hardlink': os.link(source,destination)
    else:
        try: destination.symlink_to(source)
        except OSError as error: pytest.skip(f'Host cannot create symbolic link: {error}')


@pytest.mark.parametrize('kind',['hardlink','symlink'])
@pytest.mark.parametrize('destination',['batch','temporary','screenshot','stdout'])
def test_linked_outputs_cannot_alias_inputs(tmp_path,kind,destination):
    item=target(tmp_path,VariantClass.SNV);wrapper,plan,runner=plan_for(item)
    path={'batch':plan.batch_path,'temporary':plan.batch_path.with_name('.'+plan.batch_path.name+'.hifivar.tmp'),
          'screenshot':plan.screenshots[0],'stdout':item.output_directory/'log.txt'}[destination]
    linked(item.source_vcf,path,kind)
    reject(wrapper,plan,runner,resource=item.source_vcf,overwrite=True,
           **({'stdout_path':path} if destination=='stdout' else {}))
    assert path.exists()


def test_symlink_parent_alias_is_protected(tmp_path):
    item=target(tmp_path,VariantClass.SNV);link=tmp_path/'alias'
    try: link.symlink_to(tmp_path,target_is_directory=True)
    except OSError as error: pytest.skip(f'Host cannot create directory symlink: {error}')
    wrapper,plan,runner=plan_for(item)
    reject(wrapper,plan,runner,resource=item.source_vcf,stdout_path=link/item.source_vcf.name,overwrite=True)


def test_implicit_stdout_is_protected_before_stderr_only_execution(tmp_path):
    item=target(tmp_path,VariantClass.SNV)
    # The model retains this explicit source path; stderr itself is not an input.
    source=tmp_path/'input.stdout.bam';item.alignment_path.replace(source)
    Path(str(item.alignment_path)+'.bai').replace(Path(str(source)+'.bai'))
    from dataclasses import replace
    item=replace(item,alignment_path=source)
    wrapper,plan,runner=plan_for(item)
    stderr=tmp_path/'input.bam'
    reject(wrapper,plan,runner,resource=source,stderr_path=stderr,overwrite=True)
    assert not stderr.exists()


@pytest.mark.parametrize('destination',['batch','temporary','screenshot'])
def test_output_and_log_collisions_fail_before_any_command(tmp_path,destination):
    item=target(tmp_path,VariantClass.SNV);wrapper,plan,runner=plan_for(item)
    path={'batch':plan.batch_path,'temporary':plan.batch_path.with_name('.'+plan.batch_path.name+'.hifivar.tmp'),
          'screenshot':plan.screenshots[0]}[destination]
    with pytest.raises(OutputValidationError,match='output aliases collide'):
        wrapper.run(plan,stdout_path=path,overwrite=True)
    assert not runner.calls and not item.output_directory.exists()


def test_implicit_stdout_may_not_alias_batch(tmp_path):
    item=target(tmp_path,VariantClass.SNV);batch=item.output_directory/'tool.stdout.log'
    wrapper,plan,runner=plan_for(item,batch=batch)
    with pytest.raises(OutputValidationError,match='output aliases collide'):
        wrapper.run(plan,stderr_path=item.output_directory/'tool.log',overwrite=True)
    assert not runner.calls and not item.output_directory.exists()


def test_logs_with_distinct_names_but_one_inode_are_rejected(tmp_path):
    item=target(tmp_path,VariantClass.SNV);a=tmp_path/'a.log';a.write_bytes(b'old logs');b=tmp_path/'b.log';os.link(a,b)
    wrapper,plan,runner=plan_for(item)
    with pytest.raises(OutputValidationError,match='output aliases collide'):
        wrapper.run(plan,stdout_path=a,stderr_path=b,overwrite=True)
    assert a.read_bytes()==b'old logs' and b.read_bytes()==b'old logs' and not runner.calls


def test_dangling_log_symlink_fails_closed(tmp_path):
    item=target(tmp_path,VariantClass.SNV);link=tmp_path/'dangling.log'
    try: link.symlink_to(tmp_path/'absent.log')
    except OSError as error: pytest.skip(f'Host cannot create symbolic link: {error}')
    wrapper,plan,runner=plan_for(item)
    with pytest.raises(OutputValidationError,match='dangling symlink'):
        wrapper.run(plan,stdout_path=link,overwrite=True)
    assert not runner.calls and not (tmp_path/'absent.log').exists()


def test_tool_executable_is_not_a_log_destination(tmp_path):
    item=target(tmp_path,VariantClass.SNV);program=tmp_path/'igv';program.write_bytes(b'synthetic executable')
    wrapper,plan,runner=plan_for(item,executable=str(program))
    reject(wrapper,plan,runner,resource=program,stdout_path=program,overwrite=True)


def test_existing_trgt_visualization_is_immutable(tmp_path):
    plot=tmp_path/'trgt.svg';plot.write_text('<svg/>',encoding='utf-8')
    item=target(tmp_path,VariantClass.TR,trgt_visualization_path=plot)
    wrapper,plan,runner=plan_for(item)
    reject(wrapper,plan,runner,resource=plot,stderr_path=plot,overwrite=True)


@pytest.mark.parametrize('member',['manifest.json','records.jsonl','junctions.tsv','native_result.json',
                                  'native/raw.vcf','native/evidence.tsv','native/native.vcf','native/native.vcf.tbi'])
def test_complete_bnd_bundle_members_are_immutable(tmp_path,member):
    genome,derived,bam,manifest=make_bundle(tmp_path)
    item=bnd_target(tmp_path,genome,derived,bam,bnd_bundle_manifest=manifest)
    wrapper,plan,runner=plan_for(item)
    resource=manifest.parent/member
    reject(wrapper,plan,runner,resource=resource,stdout_path=resource,overwrite=True)


def test_version_probe_retargeting_is_rechecked_before_batch_write(tmp_path):
    item=target(tmp_path,VariantClass.SNV);log=tmp_path/'late.log'
    class RetargetingRunner(FakeIgvRunner):
        def run(self,command,**kwargs):
            result=super().run(command,**kwargs)
            if '--version' in command: os.link(item.source_vcf,log)
            return result
    runner=RetargetingRunner();wrapper,plan,_=plan_for(item,runner)
    before=item.source_vcf.read_bytes()
    with pytest.raises(OutputValidationError,match='immutable input'):
        wrapper.run(plan,stdout_path=log,overwrite=True)
    assert item.source_vcf.read_bytes()==before and len(runner.calls)==1
    assert '--version' in runner.calls[0][0] and not plan.batch_path.exists()


def test_version_failure_preserves_existing_batch_and_screenshots(tmp_path):
    item=target(tmp_path,VariantClass.SNV)
    class BadVersion(FakeIgvRunner):
        def run(self,command,**kwargs): raise ToolVersionError('synthetic unavailable version')
    wrapper,plan,runner=plan_for(item,BadVersion())
    plan.batch_path.parent.mkdir(parents=True);plan.batch_path.write_bytes(b'old batch')
    plan.snapshot_directory.mkdir();plan.screenshots[0].write_bytes(b'old screenshot')
    with pytest.raises(ToolVersionError):wrapper.run(plan,overwrite=True)
    assert plan.batch_path.read_bytes()==b'old batch' and plan.screenshots[0].read_bytes()==b'old screenshot'


def test_safe_outputs_preserve_inputs_and_existing_wrapper_contract(tmp_path):
    item=target(tmp_path,VariantClass.SNV);before=item.source_vcf.read_bytes()
    wrapper,plan,runner=plan_for(item);stdout=tmp_path/'logs/out.log';stderr=tmp_path/'logs/err.log'
    result=wrapper.run(plan,stdout_path=stdout,stderr_path=stderr)
    assert result.status is IgvRunStatus.COMPLETED and result.tool_version=='2.19.8'
    assert item.source_vcf.read_bytes()==before
    assert runner.calls[-1][1]=={'stdout_path':stdout,'stderr_path':stderr,'timeout':14400,'cwd':None}
    assert plan.screenshots[0].read_bytes()==PNG  # Synthetic only; no native/image accuracy claim.


def test_dry_run_rejects_input_alias_without_command_or_output(tmp_path):
    item=target(tmp_path,VariantClass.SNV);wrapper,plan,runner=plan_for(item)
    reject(wrapper,plan,runner,resource=item.source_vcf,dry_run=True,stdout_path=item.source_vcf)
    assert not item.output_directory.exists()


@pytest.mark.parametrize('name',['stdout_path','stderr_path'])
@pytest.mark.parametrize('value',[True,17,'   '])
def test_invalid_log_paths_fail_with_typed_error_before_commands(tmp_path,name,value):
    item=target(tmp_path,VariantClass.SNV);wrapper,plan,runner=plan_for(item)
    with pytest.raises(OutputValidationError,match='IGV_OUTPUT_PROTECTION'):
        wrapper.run(plan,**{name:value})
    assert not runner.calls and not item.output_directory.exists()


def test_relative_log_alias_is_resolved_against_current_directory(tmp_path,monkeypatch):
    item=target(tmp_path,VariantClass.SNV);wrapper,plan,runner=plan_for(item)
    monkeypatch.chdir(tmp_path)
    reject(wrapper,plan,runner,resource=item.source_vcf,stdout_path=Path(item.source_vcf.name),overwrite=True)


def test_new_batch_cannot_be_an_ancestor_of_screenshot_path(tmp_path):
    item=target(tmp_path,VariantClass.SNV);wrapper,plan,runner=plan_for(item,batch=item.output_directory/'screenshots')
    with pytest.raises(OutputValidationError,match='output aliases collide'):
        wrapper.run(plan,overwrite=True)
    assert not runner.calls and not item.output_directory.exists()


@pytest.mark.parametrize('name',['review_manifest.json','review_manifest.yaml','review_manifest.tsv','.review_manifest.tsv.hifivar.tmp'])
def test_phase10_manifest_cannot_replace_selected_tr_visualization(tmp_path,name):
    from hifivar.phase10 import run_phase10
    out=tmp_path/'审阅 output';out.mkdir();plot=out/name;plot.write_bytes(b'selected immutable visualization')
    item=target(tmp_path,VariantClass.TR,trgt_visualization_path=plot)
    runner=FakeIgvRunner();before=plot.read_bytes()
    with pytest.raises(OutputValidationError,match='IGV_OUTPUT_PROTECTION'):
        run_phase10((item,),output_directory=out,igv_wrapper=IgvWrapper(runner=runner),overwrite=True)
    assert plot.read_bytes()==before and not runner.calls


def test_phase10_manifest_hardlink_to_bnd_source_is_rejected(tmp_path):
    from hifivar.phase10 import run_phase10
    genome,derived,bam,manifest=make_bundle(tmp_path)
    item=bnd_target(tmp_path,genome,derived,bam,bnd_bundle_manifest=manifest)
    item.output_directory.mkdir();os.link(manifest,item.output_directory/'review_manifest.json')
    runner=FakeIgvRunner();before=manifest.read_bytes()
    with pytest.raises(OutputValidationError,match='IGV_OUTPUT_PROTECTION'):
        run_phase10((item,),output_directory=item.output_directory,igv_wrapper=IgvWrapper(runner=runner),overwrite=True)
    assert manifest.read_bytes()==before and not runner.calls


def test_manifest_alias_created_during_igv_is_rechecked_before_serialization(tmp_path):
    from hifivar.phase10 import run_phase10
    plot=tmp_path/'plot.svg';plot.write_bytes(b'original immutable plot')
    item=target(tmp_path,VariantClass.TR,trgt_visualization_path=plot)
    class LateManifestRunner(FakeIgvRunner):
        def run(self,command,**kwargs):
            result=super().run(command,**kwargs)
            if '--batch' in command:os.link(plot,item.output_directory/'review_manifest.json')
            return result
    runner=LateManifestRunner()
    with pytest.raises(OutputValidationError,match='IGV_OUTPUT_PROTECTION'):
        run_phase10((item,),output_directory=item.output_directory,igv_wrapper=IgvWrapper(runner=runner),overwrite=True)
    assert plot.read_bytes()==b'original immutable plot'
    assert (item.output_directory/'review_manifest.json').read_bytes()==plot.read_bytes()
    assert not (item.output_directory/'review_manifest.yaml').exists()
