"""Category-neutral, fail-closed identity checks for annotation source VCFs."""
from __future__ import annotations

import gzip
import hashlib
from pathlib import Path

import pytest

from hifivar.annotation import AnnotationInput, VariantCategory, audit_annotation_vcf_identity
from hifivar.annovar import AnnovarRequest, AnnovarWrapper
from hifivar.exceptions import InputValidationError, ReferenceError
from hifivar.reference import ReferenceGenome
from hifivar.vep import VepRequest, VepWrapper


@pytest.fixture
def reference(tmp_path: Path) -> ReferenceGenome:
    fasta = tmp_path / 'reference.fa'
    fasta.write_text('>chr1\n' + 'A' * 1000 + '\n', encoding='utf-8')
    Path(f'{fasta}.fai').write_text('chr1\t1000\t6\t1000\t1001\n', encoding='utf-8')
    return ReferenceGenome.from_fasta(fasta, build='GRCh38')


def write_vcf(
    path: Path, *, samples: tuple[str, ...] = ('S1',),
    contigs: tuple[str, ...] = ('##contig=<ID=chr1,length=1000>',),
    build: str | None = 'GRCh38', format_column: bool = True,
) -> Path:
    header = ['##fileformat=VCFv4.3', *contigs]
    if build is not None:
        header.append('##reference=' + build)
    columns = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO'
    if format_column:
        columns += '\tFORMAT' + ''.join('\t' + sample for sample in samples)
    text = '\n'.join(header) + '\n' + columns + '\n'
    text += 'chr1\t10\tv1\tA\tT\t.\tPASS\t.\tGT\t0/1\n'
    if path.suffix == '.gz':
        with gzip.open(path, 'wt', encoding='utf-8', newline='') as handle:
            handle.write(text)
        Path(f'{path}.tbi').write_bytes(b'fixture-index-presence-only')
    else:
        path.write_text(text, encoding='utf-8')
    return path


@pytest.mark.parametrize('category', tuple(VariantCategory))
def test_all_categories_accept_sampled_vcf_without_svtype(tmp_path, reference, category):
    vcf = write_vcf(tmp_path / f'{category.value}.vcf')
    before = vcf.read_bytes()
    item = AnnotationInput('S1', category, vcf, 'caller', reference)
    assert item.source_identity['sample_check'] == 'PASS_EXACT_SINGLE_SAMPLE'
    assert item.source_identity['contig_length_check'] == 'PASS'
    assert item.source_identity['reference_build_check'] == 'PASS'
    assert item.source_identity['record_validation'] == 'NOT_RUN_HEADER_ONLY'
    assert item.source_identity['reference_allele_sequence_check'] == 'NOT_RUN'
    assert vcf.read_bytes() == before


@pytest.mark.parametrize(('samples', 'format_column', 'message'), [
    (('WRONG',), True, 'SAMPLE_MISMATCH'),
    ((), True, 'missing or empty'),
    (('',), True, 'missing or empty'),
    (('S1', 'EXTRA'), True, 'SAMPLE_MISMATCH'),
    (('S1', 'S1'), True, 'duplicate sample'),
    (('S1', ''), True, 'missing or empty'),
    ((' S1',), True, 'missing or empty'),
    ((), False, 'genotyped VCF'),
])
def test_header_sample_identity_is_exact(tmp_path, reference, samples, format_column, message):
    vcf = write_vcf(tmp_path / 'bad.vcf', samples=samples, format_column=format_column)
    with pytest.raises(InputValidationError, match=message):
        AnnotationInput('S1', VariantCategory.SMALL, vcf, 'caller', reference)


@pytest.mark.parametrize('sample', ['', 'S 1', 'unknown', 'UnknownSampleName'])
def test_placeholder_or_invalid_requested_sample_is_rejected(tmp_path, reference, sample):
    vcf = write_vcf(tmp_path / 'input.vcf')
    with pytest.raises(InputValidationError):
        audit_annotation_vcf_identity(vcf, sample, reference)


@pytest.mark.parametrize(('contigs', 'error', 'message'), [
    ((), InputValidationError, 'lacks contig'),
    (('##contig=<ID=chr1>',), InputValidationError, 'lacks a length'),
    (('##contig=<ID=chr1,length=0>',), InputValidationError, 'Invalid'),
    (('##contig=<ID=chr1,length=-1>',), InputValidationError, 'Invalid'),
    (('##contig=<ID=chr1,length=1.0>',), InputValidationError, 'Invalid'),
    (('##contig=<ID=chr1,length=１０００>',), InputValidationError, 'Invalid'),
    (('##contig=<ID=chr1,length=999>',), ReferenceError, 'CONTIG_LENGTH_MISMATCH'),
    (('##contig=<ID=1,length=1000>',), ReferenceError, 'REFERENCE_CONTIG_MISMATCH'),
    (('##contig=<ID=chr2,length=1000>',), ReferenceError, 'REFERENCE_CONTIG_MISMATCH'),
    (('##contig=<ID=chr1,length=1000>', '##contig=<ID=chr1,length=1000>'), InputValidationError, 'Duplicate'),
    (('##contig=<ID=chr1,length=1000,length=1000>',), InputValidationError, 'Malformed'),
    (('##contig=<ID=chr1,length=1000',), InputValidationError, 'Malformed'),
])
def test_contig_identity_fails_closed(tmp_path, reference, contigs, error, message):
    vcf = write_vcf(tmp_path / 'bad.vcf', contigs=contigs)
    with pytest.raises(error, match=message):
        audit_annotation_vcf_identity(vcf, 'S1', reference)


@pytest.mark.parametrize('build', ['GRCh37', 'hg19', '/references/GRCh37/human.fa', 'GRCh37_and_GRCh38'])
def test_explicit_wrong_build_is_rejected(tmp_path, reference, build):
    vcf = write_vcf(tmp_path / 'bad.vcf', build=build)
    with pytest.raises(ReferenceError, match='REFERENCE_BUILD_MISMATCH'):
        audit_annotation_vcf_identity(vcf, 'S1', reference)


@pytest.mark.parametrize('build', [None, '/references/genome.fa'])
def test_missing_or_unknown_reference_header_remains_unverified(tmp_path, reference, build):
    vcf = write_vcf(tmp_path / 'input.vcf', build=build)
    result = audit_annotation_vcf_identity(vcf, 'S1', reference)
    assert result['reference_build_check'] == 'UNVERIFIED_MISSING_OR_UNRECOGNIZED_METADATA'


def test_build_alias_and_quoted_contig_attributes_are_supported(tmp_path, reference):
    vcf = write_vcf(tmp_path / 'input.vcf', build='hg38',
                    contigs=('##contig=<ID=chr1,length=1000,assembly="GRCh38,primary">',))
    assert audit_annotation_vcf_identity(vcf, 'S1', reference)['reference_build_check'] == 'PASS'


def test_explicit_contig_assembly_conflict_is_not_hidden_by_reference_header(tmp_path, reference):
    vcf = write_vcf(tmp_path / 'bad.vcf', build='GRCh38',
                    contigs=('##contig=<ID=chr1,length=1000,assembly=GRCh37>',))
    with pytest.raises(ReferenceError, match='REFERENCE_BUILD_MISMATCH'):
        audit_annotation_vcf_identity(vcf, 'S1', reference)


def test_compressed_vcf_preserves_input_and_requires_existing_index(tmp_path, reference):
    vcf = write_vcf(tmp_path / 'input.vcf.gz')
    index = Path(f'{vcf}.tbi')
    before = [hashlib.sha256(path.read_bytes()).hexdigest() for path in (vcf, index)]
    assert audit_annotation_vcf_identity(vcf, 'S1', reference)['input_modified'] is False
    assert before == [hashlib.sha256(path.read_bytes()).hexdigest() for path in (vcf, index)]
    index.unlink()
    with pytest.raises(InputValidationError, match='index'):
        audit_annotation_vcf_identity(vcf, 'S1', reference)


def test_id_provenance_does_not_claim_backend_record_filtering(tmp_path, reference):
    vcf = write_vcf(tmp_path / 'input.vcf')
    item = AnnotationInput('S1', VariantCategory.SMALL, vcf, 'caller', reference, ('selected-id',))
    cache = tmp_path / 'cache'
    database = tmp_path / 'database'
    cache.mkdir()
    database.mkdir()
    annovar = AnnovarRequest(item, database, 'db1', ('refGene',), ('g',), tmp_path / 'ann', '2025-03-02')
    vep = VepRequest(item, cache, '116', 'homo_sapiens', 'GRCh38', tmp_path / 'vep.tsv')
    ann_command = AnnovarWrapper().build_command(annovar)
    vep_command = VepWrapper().build_command(vep)
    assert str(vcf.absolute()) in ann_command and str(vcf.absolute()) in vep_command
    assert 'selected-id' not in ann_command and 'selected-id' not in vep_command
    payload = item.to_dict()
    assert payload['annovar_vep_input_scope'] == 'whole_source_vcf'
    assert 'not_annovar_vep_filter' in payload['source_variant_ids_role']


@pytest.mark.parametrize('ids', [('',), ('v1', 'v1')])
def test_invalid_source_ids_remain_rejected(tmp_path, reference, ids):
    vcf = write_vcf(tmp_path / 'input.vcf')
    with pytest.raises(InputValidationError, match='IDs'):
        AnnotationInput('S1', VariantCategory.SV, vcf, 'caller', reference, ids)
