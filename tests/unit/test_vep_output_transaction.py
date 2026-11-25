"""VEP failure preservation and actual output/log alias boundaries."""
from dataclasses import replace
from pathlib import Path
import os

import pytest

from hifivar.annotation import AnnotationInput, VariantCategory
from hifivar.command import CommandResult
from hifivar.exceptions import CommandExecutionError, OutputValidationError
from hifivar.reference import ReferenceGenome
from hifivar.vep import VepRequest, VepWrapper


VEP_COLUMNS = (
    'Uploaded_variation', 'Location', 'Allele', 'Gene', 'Feature',
    'Feature_type', 'Consequence', 'cDNA_position', 'CDS_position',
    'Protein_position', 'Amino_acids', 'Codons', 'Existing_variation',
)
VEP_HEADER = ('#' + '\t'.join(VEP_COLUMNS) + '\n').encode()
VEP_ROW = b'v1\tchr1:10\tT\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\n'
VALID_VEP = VEP_HEADER + VEP_ROW


@pytest.fixture
def vep_request(tmp_path):
    fasta = tmp_path / 'reference.fa'
    fasta.write_text('>chr1\n' + 'A' * 100 + '\n', encoding='utf-8')
    Path(str(fasta) + '.fai').write_text('chr1\t100\t6\t100\t101\n', encoding='utf-8')
    source = tmp_path / 'source.vcf'
    source.write_text('##fileformat=VCFv4.3\n##contig=<ID=chr1,length=100>\n'
                     '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n'
                     'chr1\t10\tv1\tA\tT\t30\tPASS\t.\tGT\t0/1\n', encoding='utf-8')
    cache = tmp_path / 'cache'
    cache.mkdir()
    (cache / 'resource').write_bytes(b'IMMUTABLE_CACHE')
    item = AnnotationInput('S1', VariantCategory.SMALL, source, 'DeepVariant',
                           ReferenceGenome.from_fasta(fasta, build='GRCh38'))
    return VepRequest(item, cache, '115', 'homo_sapiens', 'GRCh38', tmp_path / 'result.tsv')


class Runner:
    def __init__(self, mode='success', callback=None, payload=None):
        self.mode = mode
        self.callback = callback
        self.payload = VALID_VEP if payload is None else payload
        self.calls = []

    def require_executable(self, executable):
        self.calls.append(('require', executable))
        return Path(executable)

    def run(self, command, **kwargs):
        args = tuple(map(str, command))
        self.calls.append(('run', args, kwargs))
        if kwargs.get('dry_run'):
            return CommandResult(args, 0, None, None, 0, None, False)
        if '--help' in args:
            if self.mode == 'version_failure':
                raise CommandExecutionError('version probe failed')
            if self.callback:
                self.callback('version')
            return CommandResult(args, 0, 'ensembl-vep : 115', '', .01, None, True)
        target = Path(args[args.index('--output_file') + 1])
        if self.mode != 'missing_output':
            target.write_bytes(b'' if self.mode == 'empty_output' else self.payload)
        # VEP's ordinary native stats suffix belongs to the actual output stem.
        Path(str(target) + '_summary.html').write_text('NATIVE_STATS', encoding='utf-8')
        if self.callback:
            self.callback('analysis')
        if self.mode == 'annotation_failure':
            raise CommandExecutionError('annotation failed')
        return CommandResult(args, 0, '', '', .02, None, True)


@pytest.mark.parametrize('overwrite', [False, True])
@pytest.mark.parametrize('payload', [
    b'<html>tool startup failed</html>\n',
    b'##VEP version=116.3\n',
    b'#Uploaded_variation\tConsequence\nv1\tfeature\n',
    VEP_ROW + VEP_HEADER,
    VEP_HEADER + VEP_ROW.rstrip(b'\n') + b'\textra\n',
    VEP_HEADER + b'v1\tchr1:10\n',
    VEP_HEADER + VEP_ROW.replace(b'\tT\t', b'\t\t'),
    VEP_HEADER + VEP_ROW.replace(b'v1', b'v\x00'),
    VEP_HEADER + VEP_ROW.replace(b'v1', b'v\xff'),
    VEP_HEADER + VEP_ROW + VEP_HEADER,
    VEP_HEADER + b'##unexpected metadata after the column header\n',
    VEP_HEADER.replace(b'Location', b'Uploaded_variation') + VEP_ROW,
    VEP_HEADER.replace(b'Location', b'') + VEP_ROW,
], ids=['html', 'metadata_only', 'missing_columns', 'row_before_header',
        'extra_cell', 'short_row', 'empty_cell', 'nul_byte', 'invalid_utf8',
        'repeated_header', 'late_metadata', 'duplicate_column', 'empty_column'])
def test_invalid_native_tsv_cannot_publish_or_replace(vep_request, overwrite, payload):
    before = {p: p.read_bytes() for p in (
        vep_request.input.source_vcf, vep_request.input.reference.fasta,
        vep_request.input.reference.fai, vep_request.cache_directory / 'resource',
    )}
    if overwrite:
        vep_request.output_tsv.write_bytes(b'PREVIOUS_VALIDATED_ANNOTATION')
    with pytest.raises(OutputValidationError, match='VEP_TSV_INVALID'):
        VepWrapper(runner=Runner(payload=payload)).run(
            replace(vep_request, overwrite=overwrite)
        )
    if overwrite:
        assert vep_request.output_tsv.read_bytes() == b'PREVIOUS_VALIDATED_ANNOTATION'
    else:
        assert not vep_request.output_tsv.exists()
    assert all(path.read_bytes() == data for path, data in before.items())
    staged = list(vep_request.output_tsv.parent.glob('.vep-*/result.tsv'))
    assert len(staged) == 1 and staged[0].read_bytes() == payload


@pytest.mark.parametrize('payload', [
    VALID_VEP,
    b'##VEP version=116.3\n## Column descriptions:\n' + VALID_VEP,
    VEP_HEADER,
    VALID_VEP.replace(b'\n', b'\r\n'),
    VALID_VEP.rstrip(b'\n'),
    VEP_HEADER.rstrip(b'\n') + b'\tIMPACT\tCUSTOM\n'
    + VEP_ROW.rstrip(b'\n') + b'\tMODIFIER\t-\n',
    VEP_HEADER + VEP_ROW + VEP_ROW.replace(b'intergenic_variant', b'downstream_gene_variant'),
], ids=['ordinary', 'metadata', 'header_only', 'crlf', 'no_final_newline',
        'additional_columns', 'multiple_features'])
def test_legal_native_tsv_preserved_exactly(vep_request, payload):
    result = VepWrapper(runner=Runner(payload=payload)).run(vep_request)
    assert result.status.value == 'completed'
    assert vep_request.output_tsv.read_bytes() == payload


def test_invalid_overwrite_does_not_reach_atomic_commit(vep_request, monkeypatch):
    from hifivar import vep

    def unsafe(*_):
        pytest.fail('Invalid native output must not reach publication')

    monkeypatch.setattr(vep.os, 'replace', unsafe)
    vep_request.output_tsv.write_bytes(b'PREVIOUS')
    with pytest.raises(OutputValidationError, match='VEP_TSV_INVALID'):
        VepWrapper(runner=Runner(payload=b'nonempty garbage')).run(
            replace(vep_request, overwrite=True)
        )
    assert vep_request.output_tsv.read_bytes() == b'PREVIOUS'


def test_default_columns_can_be_reordered_without_rewriting(vep_request):
    header = VEP_HEADER.decode().rstrip('\n')[1:].split('\t')
    row = VEP_ROW.decode().rstrip('\n').split('\t')
    payload = ('#' + '\t'.join(reversed(header)) + '\n'
               + '\t'.join(reversed(row)) + '\n').encode()
    VepWrapper(runner=Runner(payload=payload)).run(vep_request)
    assert vep_request.output_tsv.read_bytes() == payload


def test_tsv_validation_does_not_read_the_whole_file(vep_request, monkeypatch):
    from hifivar.vep import _validate_vep_tsv

    vep_request.output_tsv.write_bytes(VALID_VEP)

    def whole_file_read(*args, **kwargs):
        pytest.fail('Native annotation validation must stream text')

    with monkeypatch.context() as patch:
        patch.setattr(Path, 'read_text', whole_file_read)
        patch.setattr(Path, 'read_bytes', whole_file_read)
        _validate_vep_tsv(vep_request.output_tsv)


@pytest.mark.parametrize('mode', ['version_failure', 'annotation_failure', 'missing_output', 'empty_output'])
def test_failed_overwrite_keeps_previous_tsv(vep_request, mode):
    vep_request.output_tsv.write_bytes(b'PREVIOUS_ANNOTATION')
    runner = Runner(mode)
    with pytest.raises((CommandExecutionError, OutputValidationError)):
        VepWrapper(runner=runner).run(replace(vep_request, overwrite=True))
    assert vep_request.output_tsv.read_bytes() == b'PREVIOUS_ANNOTATION'


def test_implicit_stdout_alias_is_rejected_before_probe(vep_request, tmp_path):
    log = tmp_path / 'job.log'
    companion = tmp_path / 'job.stdout.log'
    os.link(vep_request.input.source_vcf, companion)
    before = vep_request.input.source_vcf.read_bytes()
    runner = Runner()
    with pytest.raises(OutputValidationError, match='ANNOTATION_OUTPUT_PROTECTION'):
        VepWrapper(runner=runner).run(vep_request, stderr_path=log)
    assert runner.calls == []
    assert vep_request.input.source_vcf.read_bytes() == before


def test_native_sidecars_do_not_follow_final_prefix_alias(vep_request):
    cache = vep_request.cache_directory / 'resource'
    os.link(cache, Path(str(vep_request.output_tsv) + '_summary.html'))
    result = VepWrapper(runner=Runner()).run(vep_request)
    assert result.status.value == 'completed'
    assert cache.read_bytes() == b'IMMUTABLE_CACHE'


def test_successful_overwrite_keeps_old_until_native_output_validated(vep_request):
    vep_request.output_tsv.write_bytes(b'PREVIOUS')

    def see_previous(_):
        assert vep_request.output_tsv.read_bytes() == b'PREVIOUS'

    runner = Runner(callback=see_previous)
    result = VepWrapper(runner=runner).run(replace(vep_request, overwrite=True))
    assert vep_request.output_tsv.read_bytes().startswith(b'#Uploaded_variation')
    actual = Path(result.command[result.command.index('--output_file') + 1])
    assert actual != vep_request.output_tsv
    assert actual.parent.parent == vep_request.output_tsv.parent
    assert actual.parent.name.startswith('.vep-')
    assert Path(str(actual) + '_summary.html').read_text() == 'NATIVE_STATS'
    assert result.artifact.output_path == vep_request.output_tsv


def test_existing_output_without_overwrite_fails_before_probe(vep_request):
    vep_request.output_tsv.write_bytes(b'EXISTING')
    runner = Runner()
    with pytest.raises(OutputValidationError, match='already exists'):
        VepWrapper(runner=runner).run(vep_request)
    assert runner.calls == []
    assert not list(vep_request.output_tsv.parent.glob('.vep-*'))


def test_version_failure_creates_no_staging_directory(vep_request):
    with pytest.raises(CommandExecutionError):
        VepWrapper(runner=Runner('version_failure')).run(vep_request)
    assert not vep_request.output_tsv.exists()
    assert not list(vep_request.output_tsv.parent.glob('.vep-*'))


def test_failed_native_output_is_retained_beside_previous_result(vep_request):
    vep_request.output_tsv.write_bytes(b'OLD')
    with pytest.raises(CommandExecutionError):
        VepWrapper(runner=Runner('annotation_failure')).run(replace(vep_request, overwrite=True))
    staged = list(vep_request.output_tsv.parent.glob('.vep-*/result.tsv'))
    assert len(staged) == 1 and staged[0].read_bytes().startswith(b'#Uploaded_variation')
    assert vep_request.output_tsv.read_bytes() == b'OLD'


@pytest.mark.parametrize('phase', ['version', 'analysis'])
def test_output_created_during_execution_is_not_replaced(vep_request, phase):
    def create_late(current):
        if phase == current:
            vep_request.output_tsv.write_bytes(b'LATE_OTHER_RESULT')

    with pytest.raises(OutputValidationError, match='already exists'):
        VepWrapper(runner=Runner(callback=create_late)).run(vep_request)
    assert vep_request.output_tsv.read_bytes() == b'LATE_OTHER_RESULT'


@pytest.mark.parametrize('phase', ['version', 'analysis'])
def test_destination_alias_created_later_fails_closed(vep_request, phase):
    protected = vep_request.cache_directory / 'resource'

    def create_alias(current):
        if phase == current:
            os.link(protected, vep_request.output_tsv)

    with pytest.raises(OutputValidationError, match='ANNOTATION_OUTPUT_PROTECTION'):
        VepWrapper(runner=Runner(callback=create_alias)).run(replace(vep_request, overwrite=True))
    assert protected.read_bytes() == b'IMMUTABLE_CACHE'


def test_implicit_log_alias_created_by_probe_is_rechecked(vep_request, tmp_path):
    companion = tmp_path / 'job.stdout.log'

    def create_alias(phase):
        if phase == 'version':
            os.link(vep_request.input.source_vcf, companion)

    runner = Runner(callback=create_alias)
    with pytest.raises(OutputValidationError, match='ANNOTATION_OUTPUT_PROTECTION'):
        VepWrapper(runner=runner).run(vep_request, stderr_path=tmp_path / 'job.log')
    assert not any(c[0] == 'run' and '--output_file' in c[1] for c in runner.calls)


def test_staged_output_alias_is_rejected_before_publication(vep_request):
    cache = vep_request.cache_directory / 'resource'
    vep_request.output_tsv.write_bytes(b'OLD')
    runner = Runner()

    def alias_stage(phase):
        if phase == 'analysis':
            args = runner.calls[-1][1]
            path = Path(args[args.index('--output_file') + 1])
            path.unlink()  # Only this fake tool's just-created test output.
            os.link(cache, path)

    runner.callback = alias_stage
    with pytest.raises(OutputValidationError, match='ANNOTATION_OUTPUT_PROTECTION'):
        VepWrapper(runner=runner).run(replace(vep_request, overwrite=True))
    assert cache.read_bytes() == b'IMMUTABLE_CACHE'
    assert vep_request.output_tsv.read_bytes() == b'OLD'


def test_exclusive_publication_rejects_commit_time_race(vep_request, monkeypatch):
    from hifivar import vep

    def raced_link(source, destination):
        Path(destination).write_bytes(b'CONCURRENT_RESULT')
        raise FileExistsError('simulated race after final precheck')

    monkeypatch.setattr(vep.os, 'link', raced_link)
    with pytest.raises(OutputValidationError, match='atomically published'):
        VepWrapper(runner=Runner()).run(vep_request)
    assert vep_request.output_tsv.read_bytes() == b'CONCURRENT_RESULT'


def test_replace_failure_preserves_old_and_valid_staged_output(vep_request, monkeypatch):
    from hifivar import vep

    def denied(*_):
        raise PermissionError('simulated replacement denial')

    monkeypatch.setattr(vep.os, 'replace', denied)
    vep_request.output_tsv.write_bytes(b'PREVIOUS')
    with pytest.raises(OutputValidationError, match='atomically published'):
        VepWrapper(runner=Runner()).run(replace(vep_request, overwrite=True))
    assert vep_request.output_tsv.read_bytes() == b'PREVIOUS'
    assert len(list(vep_request.output_tsv.parent.glob('.vep-*/result.tsv'))) == 1


def test_no_hardlink_support_has_no_unsafe_fallback(vep_request, monkeypatch):
    from hifivar import vep

    def denied(*_):
        raise OSError('hardlinks unsupported')

    def unsafe(*_):
        pytest.fail('Must not fall back to replacing an unowned output')

    monkeypatch.setattr(vep.os, 'link', denied)
    monkeypatch.setattr(vep.os, 'replace', unsafe)
    with pytest.raises(OutputValidationError, match='atomically published'):
        VepWrapper(runner=Runner()).run(vep_request)
    assert not vep_request.output_tsv.exists()


def test_dry_run_uses_requested_final_path_without_side_effects(vep_request):
    runner = Runner()
    result = VepWrapper(runner=runner).run(vep_request, dry_run=True)
    assert result.command[result.command.index('--output_file') + 1] == str(vep_request.output_tsv.absolute())
    assert len(runner.calls) == 1 and runner.calls[0][2]['dry_run'] is True
    assert not vep_request.output_tsv.exists()
    assert not list(vep_request.output_tsv.parent.glob('.vep-*'))


def test_repeated_requests_have_unique_owned_staging(vep_request):
    paths = []
    for name in ['a.tsv', 'b.tsv']:
        request = replace(vep_request, output_tsv=vep_request.output_tsv.parent / name)
        result = VepWrapper(runner=Runner()).run(request)
        paths.append(Path(result.command[result.command.index('--output_file') + 1]))
    assert len({path.parent for path in paths}) == 2
    assert all(path.exists() for path in paths)
