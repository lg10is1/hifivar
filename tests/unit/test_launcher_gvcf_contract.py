"""Read-only gVCF prefix evidence, separate from full GLnexus validation."""
from __future__ import annotations

import hashlib
import json
import struct
from pathlib import Path
from types import SimpleNamespace
import zlib

import pytest
import yaml

from hifivar import launcher
from hifivar.config import load_config
from hifivar.exceptions import InputValidationError, ReferenceError


DEFAULT = Path(__file__).parents[2] / 'src/hifivar/resources/configs/default.yaml'


def record(alt='<*>', *, position='100', end='150', gt='0/0', pl='0,30,300',
           contig='chr1', ref='A', info=None, format_keys='GT:PL', sample=None):
    fields = [contig, position, '.', ref, alt, '.', '.',
              f'END={end}' if info is None else info, format_keys,
              f'{gt}:{pl}' if sample is None else sample]
    return '\t'.join(fields) + '\n'


def text(*, records=None, alt_declaration=False, sample='S1', contig='chr1', length=1000,
         build='GRCh38', metadata=None):
    default = ['##INFO=<ID=END,Number=1,Type=Integer,Description="Reference block, inclusive end">',
               '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
               '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred likelihoods">']
    lines = ['##fileformat=VCFv4.2', f'##reference={build}',
             f'##contig=<ID={contig},length={length}>', '##DeepVariant_version=1.10.0']
    lines.extend(default if metadata is None else metadata)
    if alt_declaration:
        lines.append('##ALT=<ID=NON_REF,Description="Any non-reference allele">')
    lines.append(f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}')
    return '\n'.join(lines) + '\n' + ''.join([record()] if records is None else records)


def bgzf(data: bytes) -> bytes:
    compressor = zlib.compressobj(wbits=-15)
    payload = compressor.compress(data) + compressor.flush()
    size = 18 + len(payload) + 8
    return (b'\x1f\x8b\x08\x04' + b'\0' * 4 + b'\0\xff\x06\0BC\x02\0'
            + struct.pack('<H', size - 1) + payload + struct.pack('<II', zlib.crc32(data), len(data)))


@pytest.fixture
def context():
    return SimpleNamespace(reference=SimpleNamespace(contigs=(SimpleNamespace(name='chr1', length=1000),)))


def source(tmp_path, content):
    path = tmp_path / 'source.g.vcf.gz'
    path.write_bytes(bgzf(content.encode()))
    # Structural index fixture only; unit tests do not claim tabix query success.
    Path(f'{path}.tbi').write_bytes(bgzf(b'TBI\x01structural-fixture'))
    return path


@pytest.mark.parametrize('alt,declaration', [('<*>', False), ('<NON_REF>', True)])
def test_both_gvcf_reference_alleles_have_explicit_bounded_evidence(tmp_path, context, alt, declaration):
    path = source(tmp_path, text(records=[record('C', info='.'), record(alt)], alt_declaration=declaration))
    before = hashlib.sha256(path.read_bytes()).hexdigest()
    audit = launcher._audit_ready_gvcf(path, context)
    assert audit['scope'] == 'TYPED_HEADER_AND_BOUNDED_RECORD_PREFIX'
    assert audit['record_prefix_limit'] == 4096 and audit['records_examined'] == 2
    assert audit['reached_eof'] is True and audit['reference_blocks_in_prefix'] == 1
    assert audit['first_reference_block']['prefix_record_number'] == 2
    assert audit['first_reference_block']['alt'] == alt
    assert audit['full_record_validation'] == audit['reference_allele_sequence_check'] == 'NOT_RUN'
    assert audit['input_modified'] is False
    assert hashlib.sha256(path.read_bytes()).hexdigest() == before


def test_deepvariant_style_block_needs_no_alt_declaration(tmp_path, context):
    content = text(records=[record(format_keys='GT:GQ:MIN_DP:PL', sample='0/0:50:28:0,90,900')])
    content = content.replace('##FORMAT=<ID=PL,', '##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum block depth">\n##FORMAT=<ID=PL,')
    assert '##ALT=' not in content
    audit = launcher._audit_ready_gvcf(source(tmp_path, content), context)
    assert audit['first_reference_block']['alt'] == '<*>'


def test_metadata_attribute_order_and_quoted_description_do_not_change_types(tmp_path, context):
    content = text().replace(
        'ID=END,Number=1,Type=Integer,Description="Reference block, inclusive end"',
        'Description="Reference block, Type=String is not an attribute",Type=Integer,Number=1,ID=END')
    assert launcher._audit_ready_gvcf(source(tmp_path, content), context)['reference_blocks_in_prefix'] == 1


def test_no_call_reference_block_can_have_missing_likelihoods(tmp_path, context):
    audit = launcher._audit_ready_gvcf(source(tmp_path, text(records=[record(gt='./.', pl='.')])), context)
    assert audit['reference_blocks_in_prefix'] == 1


@pytest.mark.parametrize('content', [
    text(records=[record('C', info='.')]),
    text(records=[record('C', info='.')], alt_declaration=True),
    text(records=[record('<DEL>')], alt_declaration=True),
    text(records=[record('<NON_REF_BAD>')]),
    text(records=[record('*')]),
    text(records=[]),
])
def test_ordinary_or_wrong_marker_vcf_does_not_become_gvcf_from_metadata(tmp_path, context, content):
    with pytest.raises(InputValidationError, match='GVCF_STRUCTURE_UNVERIFIED.*ordinary variant-only'):
        launcher._audit_ready_gvcf(source(tmp_path, content), context)


def test_reference_block_after_cap_is_unverified_not_claimed_globally_invalid(tmp_path, context, monkeypatch):
    monkeypatch.setattr(launcher, '_GVCF_RECORD_PREFIX_LIMIT', 2)
    path = source(tmp_path, text(records=[record('C', info='.'), record('G', info='.'), record()]))
    with pytest.raises(InputValidationError, match=r'GVCF_STRUCTURE_UNVERIFIED.*2/2'):
        launcher._audit_ready_gvcf(path, context)


@pytest.mark.parametrize('metadata', [
    [],
    ['##INFO=<ID=END,Number=1,Type=String,Description="Wrong">',
     '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
     '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="PL">'],
    ['##INFO=<ID=END,Number=1,Type=Integer,Description="END">',
     '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
     '##FORMAT=<ID=PL,Number=1,Type=Integer,Description="Wrong cardinality">'],
    ['##INFO=<ID=END,Number=1,Type=Integer,Description="END">',
     '##FORMAT=<ID=GT,Number=1,Type=Integer,Description="Wrong type">',
     '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="PL">'],
])
def test_typed_end_gt_pl_metadata_is_required(tmp_path, context, metadata):
    with pytest.raises(InputValidationError, match='required declarations'):
        launcher._audit_ready_gvcf(source(tmp_path, text(metadata=metadata)), context)


def test_duplicate_definition_does_not_silently_override_type(tmp_path, context):
    content = text().replace('#CHROM', '##INFO=<ID=END,Number=1,Type=String,Description="Duplicate">\n#CHROM')
    with pytest.raises(InputValidationError, match='Duplicate gVCF'):
        launcher._audit_ready_gvcf(source(tmp_path, content), context)


def test_non_ref_block_without_alt_declaration_fails(tmp_path, context):
    with pytest.raises(InputValidationError, match='NON_REF.*ALT declaration'):
        launcher._audit_ready_gvcf(source(tmp_path, text(records=[record('<NON_REF>')])), context)


@pytest.mark.parametrize('updates,match', [
    ({'end': '99'}, 'END/span'), ({'end': '1001'}, 'END/span'),
    ({'end': '-1'}, 'END/span'), ({'end': 'NaN'}, 'END/span'),
    ({'info': '.'}, 'END/span'), ({'info': 'END=120;END=150'}, 'END/span'),
    ({'ref': 'AA'}, 'END/span'), ({'gt': '0/1'}, 'cannot call'),
    ({'gt': '1/1'}, 'cannot call'), ({'gt': '0/a'}, 'cannot call'),
    ({'pl': '0,bad,30'}, 'PL values'), ({'pl': '.'}, 'PL values'),
    ({'format_keys': 'GT', 'sample': '0/0'}, 'GT and PL'),
    ({'format_keys': 'GT:PL:PL', 'sample': '0/0:0,30,300:0,30,300'}, 'GT and PL'),
    ({'sample': '0/0'}, 'GT and PL'), ({'contig': '1'}, 'REFERENCE_CONTIG_MISMATCH'),
    ({'position': '0'}, 'record position'), ({'position': '1001'}, 'record position'),
])
def test_malformed_reference_blocks_fail_without_mutating_input(tmp_path, context, updates, match):
    path = source(tmp_path, text(records=[record(**updates)]))
    before = path.read_bytes()
    with pytest.raises(InputValidationError, match=match):
        launcher._audit_ready_gvcf(path, context)
    assert path.read_bytes() == before


def test_line_bound_is_explicit(tmp_path, context, monkeypatch):
    monkeypatch.setattr(launcher, '_GVCF_LINE_CHARACTER_LIMIT', 50)
    with pytest.raises(InputValidationError, match='bounded line limit'):
        launcher._audit_ready_gvcf(source(tmp_path, text()), context)


@pytest.fixture
def ready_config(tmp_path, monkeypatch):
    (tmp_path / 'ref.fa').write_text('>chr1\n' + 'A' * 1000 + '\n')
    (tmp_path / 'ref.fa.fai').write_text('chr1\t1000\t6\t1000\t1001\n')
    (tmp_path / 'S1.bam').write_bytes(b'unused-metadata-context-unit-fixture')
    (tmp_path / 'samples.tsv').write_text('sample_id\tinput\nS1\tS1.bam\n')
    path = source(tmp_path, text())
    fields = 'sample\ttrack\tstate\tsource_path\tindex_path\tsource_tool\tsource_version\treference_build\tcatalog_id\n'
    (tmp_path / 'cohort.tsv').write_text(fields + f'S1\tsmall_variants\tCALLED\t{path.name}\t{path.name}.tbi\tDeepVariant\t1.10.0\tGRCh38\t\n')
    user = tmp_path / 'user.yaml'
    user.write_text(yaml.safe_dump({'reference': {'fasta': 'ref.fa', 'build': 'GRCh38'},
        'samples': {'sheet': 'samples.tsv'}, 'runtime': {'threads': 1}, 'workflow': {'memory_mb': 2048},
        'cohort': {'enabled': True, 'cohort_id': 'C1', 'input_manifest': 'cohort.tsv',
                   'small_variants': {'enabled': True, 'threads': 1, 'memory_gb': 1}}}))
    workflow = tmp_path / 'workflow'
    workflow.mkdir()
    (workflow / 'Snakefile').write_text('rule all:\n    input: []\n')
    monkeypatch.setattr(launcher.platform, 'system', lambda: 'Linux')
    monkeypatch.setattr(launcher, 'installed_workflow_root', lambda: workflow)
    return load_config(DEFAULT, user_config=user)


def test_prepare_launch_persists_scope_but_does_not_transform_ready_gvcf(ready_config, tmp_path):
    path = tmp_path / 'source.g.vcf.gz'
    before = path.read_bytes()
    effective, _, _, inputs = launcher.prepare_launch(ready_config, tmp_path / 'run')
    assert path in inputs
    assert Path(effective['cohort']['input_manifest']) == tmp_path / 'run/cohort.inputs.tsv'
    assert path.read_bytes() == before and not (tmp_path / 'run').exists()

    class SchedulerBoundary:
        def run(self, command, **kwargs):
            out = '9.26.1\n' if '--version' in command else ''
            for key in ('stdout_path', 'stderr_path'):
                if kwargs.get(key):
                    Path(kwargs[key]).write_text(out if key == 'stdout_path' else '')
            return SimpleNamespace(stdout=out, returncode=0, duration_seconds=0.001)

    assert launcher.WorkflowLauncher(SchedulerBoundary()).run(
        ready_config, run_dir=tmp_path / 'run', dry_run=True) == 'DRY_RUN_PASS'
    audit = json.loads(next((tmp_path / 'run/attempts').glob('*/downstream_input_audit.json')).read_text())
    evidence = audit['cohort'][0]['gvcf_structure_audit']
    assert evidence['records_examined'] == 1 and evidence['first_reference_block']['alt'] == '<*>'
    assert evidence['full_record_validation'] == 'NOT_RUN'
    assert path.read_bytes() == before


@pytest.mark.parametrize('updates', [{'sample': 'OTHER'}, {'contig': '1'}, {'length': 999}, {'build': 'GRCh37'}])
def test_exact_sample_contig_length_and_build_guards_are_not_bypassed(ready_config, tmp_path, updates):
    path = source(tmp_path, text(**updates))
    before = path.read_bytes()
    with pytest.raises((InputValidationError, ReferenceError)):
        launcher.prepare_launch(ready_config, tmp_path / 'run')
    assert not (tmp_path / 'run').exists() and path.read_bytes() == before


def test_missing_conventional_index_still_rejected(ready_config, tmp_path):
    (tmp_path / 'source.g.vcf.gz.tbi').unlink()
    with pytest.raises((InputValidationError, ReferenceError)):
        launcher.prepare_launch(ready_config, tmp_path / 'run')
    assert not (tmp_path / 'run').exists()
