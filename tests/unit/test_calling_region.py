"""Bounded DeepVariant request tests; no scientific calls are simulated as truth."""
from dataclasses import replace

import pytest

from hifivar.calling_region import parse_calling_region, validate_calling_region
from hifivar.config import validate_config
from hifivar.deepvariant import DeepVariantExecutionMode, DeepVariantRuntime, DeepVariantWrapper
from hifivar.exceptions import ConfigurationError, InputValidationError, ReferenceError
from hifivar.phase3 import Phase3Settings, run_phase3
from test_deepvariant import make_request, FakeCommandRunner
from test_phase3 import make_context, DryRunner


@pytest.mark.parametrize('value', ['', 'chr1', 'chr1:0-2', 'chr1:2-1', 'chr1:01-4',
    'chr1:1-4 ', ' chr1:1-4', 'chr1:1-4 chr2:1-2', '/tmp/regions.bed',
    'chr1:1,000-2,000', '--regions=chr1:1-4', True, 4, ['chr1:1-4']])
def test_invalid_syntax_is_not_normalized(value):
    with pytest.raises(InputValidationError):
        parse_calling_region(value)
    with pytest.raises(ConfigurationError):
        validate_config({'small': {'region': value}}, require_complete=False)


def test_exact_bounds_and_missing_contig_fail_before_writes(tmp_path):
    req = make_request(tmp_path)
    assert validate_calling_region('chr1:1-4', req.artifact.reference).end == 4
    assert validate_calling_region('chr1:4-4', req.artifact.reference).start == 4
    for value in ('1:1-4', 'chr1:1-5'):
        with pytest.raises(ReferenceError):
            replace(req, region=value)
    assert not req.output_vcf.parent.exists()


@pytest.mark.parametrize('mode', list(DeepVariantExecutionMode))
def test_default_command_unchanged_and_one_region_argument(mode, tmp_path):
    req = make_request(tmp_path)
    runtime = DeepVariantRuntime(mode=mode, image=None if mode.value == 'native' else 'deepvariant:test')
    wrapper = DeepVariantWrapper(runtime=runtime)
    whole = wrapper.build_command(req)
    regional = wrapper.build_command(replace(req, region='chr1:1-4'))
    assert not any(a.startswith('--regions') for a in whole)
    assert regional == whole + ['--regions=chr1:1-4']
    assert not req.output_vcf.parent.exists()
    assert req.to_dict()['calling_region'] is None
    assert req.to_dict()['requested_calling_scope'] == 'whole_reference'


def test_region_request_and_result_provenance(tmp_path):
    req = replace(make_request(tmp_path), region='chr1:1-4')
    runner = FakeCommandRunner()
    result = DeepVariantWrapper(runner=runner).run(req, dry_run=True)
    metadata = result.to_dict()
    assert metadata['request']['calling_region'] == 'chr1:1-4'
    assert metadata['request']['requested_calling_scope'] == 'region_limited'
    assert metadata['command']['args'][-1] == '--regions=chr1:1-4'
    assert not req.output_vcf.parent.exists()


def test_phase3_region_settings_propagate(tmp_path):
    context = make_context(tmp_path)
    settings = replace(Phase3Settings.from_config(context.config), region='chr1:2-4')
    runner = DryRunner()
    report = run_phase3(context, tmp_path/'results', settings=settings,
                        deepvariant_wrapper=DeepVariantWrapper(runner=runner), dry_run=True)
    assert report.to_dict()['settings']['calling_region'] == 'chr1:2-4'
    assert report.sample_results[0].call.request.region == 'chr1:2-4'
    assert '--regions=chr1:2-4' in runner.commands[0]
    config = dict(context.config)
    config['small'] = dict(config['small'], region='chr1:2-4')
    assert Phase3Settings.from_config(config).region == 'chr1:2-4'


def test_region_is_not_silently_ignored_by_phase15():
    with pytest.raises(ConfigurationError, match='Phase 15'):
        validate_config({'small': {'region': 'chr1:1-4'}, 'low_depth': {'enabled': True}}, require_complete=False)
