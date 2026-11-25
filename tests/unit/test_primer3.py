"""Native-adapter contract with a fake runner; no wet-lab or native claim."""
from pathlib import Path
import json

import pytest

from hifivar.command import CommandResult
from hifivar.exceptions import CommandExecutionError, OutputValidationError, ToolVersionError
from hifivar.primer3 import Primer3Wrapper
from test_primer_design import template, response


class FakePrimer3Runner:
    def __init__(self, item):self.item=item;self.calls=[];self.mode=None
    def require_executable(self, executable):return Path(executable)
    def run(self,command,**kwargs):
        self.calls.append((command,kwargs))
        if '--about' in command:
            out='no version' if self.mode=='version' else 'libprimer3 release 2.6.1\n'
            return CommandResult(tuple(command),0,out,'',0.001,None,True)
        if self.mode=='command':raise CommandExecutionError('synthetic failure')
        tags=response(self.item)
        if self.mode=='identity':tags['SEQUENCE_ID']='wrong'
        if self.mode=='template':tags['SEQUENCE_TEMPLATE']='ACGT'
        if self.mode=='invalid_utf8':content=b'\xff'
        elif self.mode=='empty':content=b''
        else:content=('\n'.join(f'{k}={v}' for k,v in tags.items())+'\n=\n').encode()
        Path(kwargs['stdout_path']).write_bytes(content)
        Path(kwargs['stderr_path']).write_text('')
        return CommandResult(tuple(command),0,None,None,0.01,kwargs.get('cwd'),True,
                             Path(kwargs['stdout_path']),Path(kwargs['stderr_path']))


def test_wrapper_runs_shell_free_and_persists_exact_design(tmp_path):
    item=template();runner=FakePrimer3Runner(item);out=tmp_path/'中文 unit'
    result=Primer3Wrapper(runner=runner).design(item,out)
    assert result['candidate_pair_count']==3
    assert json.loads((out/'primer_design.json').read_text())==result
    command,options=runner.calls[-1]
    assert command==['primer3_core','--strict_tags','--io_version=4',str(out/'primer3.input.txt')]
    assert options['timeout']==120 and options['cwd']==out
    assert result['provenance']['tool_version']=='2.6.1'
    assert result['provenance']['returncode']==0
    assert result['specificity_status']=='NOT_ASSESSED'


def test_dry_run_has_no_tool_probe_or_files(tmp_path):
    item=template();runner=FakePrimer3Runner(item);out=tmp_path/'planned'
    result=Primer3Wrapper(runner=runner).design(item,out,dry_run=True)
    assert result['status']=='PLANNED' and not runner.calls and not out.exists()


def test_existing_outputs_never_overwritten(tmp_path):
    item=template();runner=FakePrimer3Runner(item);out=tmp_path/'existing';out.mkdir()
    original=out/'primer3.input.txt';original.write_bytes(b'keep me')
    with pytest.raises(OutputValidationError,match='already exists'):
        Primer3Wrapper(runner=runner).design(item,out)
    assert original.read_bytes()==b'keep me' and not runner.calls


@pytest.mark.parametrize('mode',['identity','template','empty','invalid_utf8'])
def test_incomplete_or_wrong_native_outputs_fail_preserving_logs(tmp_path,mode):
    item=template();runner=FakePrimer3Runner(item);runner.mode=mode;out=tmp_path/mode
    with pytest.raises(OutputValidationError):Primer3Wrapper(runner=runner).design(item,out)
    assert (out/'primer3.input.txt').is_file() and (out/'primer3.output.txt').is_file()
    assert not (out/'primer_design.json').exists()


def test_tool_failure_preserves_input_and_never_writes_success(tmp_path):
    item=template();runner=FakePrimer3Runner(item);runner.mode='command';out=tmp_path/'failed'
    with pytest.raises(CommandExecutionError):Primer3Wrapper(runner=runner).design(item,out)
    assert (out/'primer3.input.txt').is_file() and not (out/'primer_design.json').exists()


def test_bad_version_does_not_create_output(tmp_path):
    item=template();runner=FakePrimer3Runner(item);runner.mode='version';out=tmp_path/'absent'
    with pytest.raises(ToolVersionError):Primer3Wrapper(runner=runner).design(item,out)
    assert not out.exists()
