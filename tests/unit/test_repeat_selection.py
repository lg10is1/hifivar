"""TRGT's real missing-ID convention and anchor-excluding overlap contracts."""
from dataclasses import replace
import gzip
import hashlib
import json
from pathlib import Path

import pytest

from hifivar.annotation import (AnnotationInput, RegionCategory, RegionDatabase,
    VariantCategory, annotate_region_overlaps, read_selected_variant_loci)
from hifivar.exceptions import InputValidationError, OutputValidationError
from hifivar.reference import ReferenceGenome
from hifivar.repeat_selection import read_selected_repeat_records


@pytest.fixture
def case(tmp_path):
    fasta = tmp_path / 'reference.fa'
    fasta.write_text('>chr1\n' + 'A' * 1000 + '\n')
    Path(str(fasta) + '.fai').write_text('chr1\t1000\t6\t1000\t1001\n')
    reference = ReferenceGenome.from_fasta(fasta, build='GRCh38')
    path = tmp_path / 'TRGT.vcf'
    path.write_text('##fileformat=VCFv4.3\n##contig=<ID=chr1,length=1000>\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n'
        'chr1\t100\t.\tAAAAA\tAAA\t.\t.\tTRID=L1;END=104;MOTIFS=A;STRUC=<TR>\tGT:AL:SD\t1/1:2,2:8,9\n'
        'chr1\t200\tother\tAAAAA\t.\t.\t.\tTRID=L2;END=204;MOTIFS=A;STRUC=<TR>\tGT:AL:SD\t0/0:4,4:9,9\n')
    item = AnnotationInput('S1', VariantCategory.TR, path, 'trgt', reference, ('TRID:L1',))
    return item, tmp_path


def test_missing_vcf_id_selects_explicit_trid_without_rewriting(case):
    item, _ = case
    before = item.source_vcf.read_bytes()
    loci = read_selected_variant_loci(item)
    assert len(loci) == 1
    locus = loci[0]
    assert (locus.variant_id, locus.start, locus.end) == ('TRID:L1', 101, 104)
    meta = locus.repeat_record.to_dict()
    assert meta['raw_vcf_id'] == '.' and meta['trid'] == 'L1'
    assert meta['record_ordinal'] == 1 and meta['source_sha256'] == hashlib.sha256(before).hexdigest()
    assert item.source_vcf.read_bytes() == before


def test_selection_order_and_nonmissing_raw_id_preserved(case):
    item, _ = case
    selected = read_selected_repeat_records(item.source_vcf, ('TRID:L2', 'TRID:L1'),
        reference=item.reference, sample_id='S1')
    assert [x.raw_vcf_id for x in selected] == ['other', '.']
    assert [x.ordinal for x in selected] == [2, 1]


def test_revalidation_hashes_shared_source_once_per_boundary(case, monkeypatch):
    from hifivar import repeat_selection
    item, _ = case
    records = read_selected_repeat_records(item.source_vcf, ('TRID:L1', 'TRID:L2'),
        reference=item.reference, sample_id='S1')
    original, calls = repeat_selection._sha, []
    def observed(path):
        calls.append(path)
        return original(path)
    monkeypatch.setattr(repeat_selection, '_sha', observed)
    repeat_selection.revalidate_repeat_records(records)
    assert calls == [item.source_vcf]


def test_gzip_stream_and_no_call_not_promoted(case):
    item, root = case
    text = item.source_vcf.read_text().replace('1/1:2,2:8,9', './.:.,.:0,0')
    compressed = root / 'input.vcf.gz'
    with gzip.open(compressed, 'wt') as f:
        f.write(text)
    # Synthetic parser fixture only; native index readability is a Linux gate.
    Path(str(compressed) + '.tbi').write_bytes(b'SYNTHETIC_INDEX_PRESENCE_NOT_TABIX')
    selected = read_selected_repeat_records(compressed, ('TRID:L1',), reference=item.reference, sample_id='S1')
    assert selected[0].trid == 'L1'  # interval inspection is not a CALLED/truth assertion


def test_compressed_source_without_index_is_rejected(case):
    item, root = case
    compressed = root / 'missing_index.vcf.gz'
    with gzip.open(compressed, 'wt') as f:
        f.write(item.source_vcf.read_text())
    with pytest.raises(InputValidationError, match='index missing'):
        read_selected_repeat_records(compressed, ('TRID:L1',), reference=item.reference, sample_id='S1')


def test_legacy_tr_ordinary_id_route_remains_unchanged(case):
    item, _ = case
    ordinary = replace(item, source_variant_ids=('other',))
    locus = read_selected_variant_loci(ordinary)[0]
    assert (locus.variant_id, locus.start, locus.end) == ('other', 200, 204)
    assert locus.repeat_record is None


@pytest.mark.parametrize('selectors', [(), ('L1',), ('TRID:',), ('TRID:.',), ('TRID:L1', 'TRID:L1'),
    ('TRID:L1', 'other'), ('TRID:L 1',), ('TRID:L,1',), ('TRID:L1', []), 'TRID:L1', None])
def test_invalid_or_mixed_selectors_fail(case, selectors):
    item, _ = case
    with pytest.raises(InputValidationError):
        read_selected_repeat_records(item.source_vcf, selectors, reference=item.reference, sample_id='S1')


@pytest.mark.parametrize('old,new', [('END=104', 'END=100'), ('END=104', 'END=1001'),
    ('END=104', 'END=x'), ('END=104', 'MISSING_END=104'),
    ('MOTIFS=A;', ''), ('STRUC=<TR>', 'OTHER=x'), ('AAAAA\tAAA', 'AAAA\tAAA'),
    ('TRID=L1;', 'TRID=L1;TRID=L2;'), ('END=104;', 'END=104;END=105;')])
def test_bad_selected_repeat_contract_fails(case, old, new):
    item, _ = case
    item.source_vcf.write_text(item.source_vcf.read_text().replace(old, new))
    with pytest.raises(InputValidationError):
        read_selected_repeat_records(item.source_vcf, ('TRID:L1',), reference=item.reference, sample_id='S1')


def test_duplicate_selected_trid_is_ambiguous(case):
    item, _ = case
    text = item.source_vcf.read_text()
    item.source_vcf.write_text(text + text.splitlines(keepends=True)[3])
    with pytest.raises(InputValidationError, match='duplicate'):
        read_selected_variant_loci(item)


def test_no_ordinary_id_fallback(case):
    item, _ = case
    item.source_vcf.write_text(item.source_vcf.read_text().replace('TRID=L1;', 'TRID=OTHER;').replace('\t.\tAAAAA', '\tL1\tAAAAA', 1))
    with pytest.raises(InputValidationError, match='absent'):
        read_selected_variant_loci(item)


def test_header_sample_mismatch_fails(case):
    item, _ = case
    item.source_vcf.write_text(item.source_vcf.read_text().replace('\tS1\n', '\tOTHER\n'))
    with pytest.raises(InputValidationError):
        read_selected_variant_loci(item)


def test_explicit_namespace_requires_trgt_source(case):
    item, _ = case
    with pytest.raises(InputValidationError, match='TRGT'):
        read_selected_variant_loci(replace(item, source_tool='other-caller'))


def test_repeat_locus_cannot_be_relabelled_or_shifted(case):
    item, _ = case
    locus = read_selected_variant_loci(item)[0]
    with pytest.raises(InputValidationError, match='differs'):
        replace(locus, start=100)
    with pytest.raises(InputValidationError, match='differs'):
        replace(locus, variant_category=VariantCategory.SV)


def test_anchor_and_end_boundary_excluded_and_traceable(case):
    item, root = case
    bed = root / 'features.bed'
    bed.write_text('chr1\t99\t100\tanchor_only\nchr1\t100\t104\trepeat_core\nchr1\t104\t105\tafter_repeat\n')
    result = annotate_region_overlaps(read_selected_variant_loci(item),
        [RegionDatabase(RegionCategory.REPEAT, bed, 'fixed-1', 'GRCh38')],
        reference=item.reference, output_tsv=root / 'overlap.tsv')
    assert [x.feature_id for x in result.overlaps] == ['repeat_core']
    row = result.overlaps[0].to_dict()
    assert json.loads(row['repeat_json'])['raw_vcf_id'] == '.'
    assert row['coordinate_semantics'] == 'TRGT_REPEAT_1_BASED_INCLUSIVE_EXCLUDING_LEFT_ANCHOR'
    assert result.to_dict()['repeat_annotations'][0]['overlap_count'] == 1


def test_zero_overlap_keeps_selected_repeat_provenance(case):
    item, root = case
    bed = root / 'empty_hit.bed'
    bed.write_text('chr1\t99\t100\tanchor_only\n')
    result = annotate_region_overlaps(read_selected_variant_loci(item),
        [RegionDatabase(RegionCategory.GENE, bed, 'fixed-1', 'GRCh38')],
        reference=item.reference, output_tsv=root / 'zero.tsv')
    assert not result.overlaps and result.repeat_annotations[0]['overlap_count'] == 0
    assert 'repeat_json' in (root / 'zero.tsv').read_text()


def test_changed_source_fails_before_output(case):
    item, root = case
    loci = read_selected_variant_loci(item)
    item.source_vcf.write_text(item.source_vcf.read_text() + '\n')
    with pytest.raises(InputValidationError, match='changed'):
        annotate_region_overlaps(loci, (), reference=item.reference, output_tsv=root / 'new.tsv')
    assert not (root / 'new.tsv').exists()


def test_output_cannot_replace_selected_source(case):
    item, _ = case
    before = item.source_vcf.read_bytes()
    with pytest.raises(OutputValidationError):
        annotate_region_overlaps(read_selected_variant_loci(item), (), reference=item.reference,
            output_tsv=item.source_vcf, overwrite=True)
    assert item.source_vcf.read_bytes() == before


def test_sv_id_with_prefix_still_uses_ordinary_id_contract(case):
    item, _ = case
    item.source_vcf.write_text(item.source_vcf.read_text().replace('\t.\tAAAAA', '\tTRID:SV\tAAAAA', 1))
    sv = AnnotationInput('S1', VariantCategory.SV, item.source_vcf, 'sv-caller', item.reference, ('TRID:SV',))
    locus = read_selected_variant_loci(sv)[0]
    assert locus.start == 100 and locus.repeat_record is None
