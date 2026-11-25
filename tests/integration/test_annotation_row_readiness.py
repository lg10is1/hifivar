"""Ready-row annotation must not validate unrelated future caller outputs."""
from __future__ import annotations

from pathlib import Path
import runpy
from types import SimpleNamespace

import pytest

from hifivar.annotation import (
    AnnotationArtifact, AnnotationResult, AnnotationRunStatus, AnnotationSource,
    VariantCategory, read_annotation_input_rows, read_annotation_inputs,
)
from hifivar.context import AnalysisContext
from hifivar.exceptions import InputValidationError
from hifivar.reference import ReferenceGenome
import hifivar.vep


@pytest.fixture
def prepared(tmp_path):
    fasta = tmp_path / 'reference.fa'
    fasta.write_text('>chr1\n' + 'A' * 1000 + '\n', encoding='utf-8')
    Path(f'{fasta}.fai').write_text('chr1\t1000\t6\t1000\t1001\n', encoding='utf-8')
    reference = ReferenceGenome.from_fasta(fasta, build='GRCh38')
    source = tmp_path / 'S1.small.vcf'
    source.write_text(
        '##fileformat=VCFv4.3\n##contig=<ID=chr1,length=1000>\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n'
        'chr1\t10\tv1\tA\tT\t.\tPASS\t.\tGT\t0/1\n', encoding='utf-8',
    )
    manifest = tmp_path / 'inputs.tsv'
    manifest.write_text(
        'sample\tvariant_category\tsource_vcf\tsource_tool\tsource_variant_ids\n'
        'S1\tsmall\tS1.small.vcf\tdeepvariant\tv1\n'
        'S2\tsv\tfuture.S2.sv.vcf\tjasmine\tv2\n', encoding='utf-8',
    )
    return reference, source, manifest


def test_selection_occurs_before_future_source_validation(prepared):
    reference, source, manifest = prepared
    assert len(read_annotation_input_rows(manifest)) == 2
    selected = read_annotation_inputs(manifest, reference=reference, sample_id='S1', variant_category='small')
    assert len(selected) == 1 and selected[0].source_vcf == source
    assert selected[0].source_variant_ids == ('v1',)
    with pytest.raises(InputValidationError, match='future.S2.sv.vcf'):
        read_annotation_inputs(manifest, reference=reference)
    with pytest.raises(InputValidationError, match='future.S2.sv.vcf'):
        read_annotation_inputs(manifest, reference=reference, sample_id='S2', variant_category=VariantCategory.SV)


@pytest.mark.parametrize('filters', [
    {'sample_id': 'S1'}, {'variant_category': 'small'},
    {'sample_id': 'S1', 'variant_category': VariantCategory.SMALL},
])
def test_each_optional_filter_preserves_all_row_default(prepared, filters):
    reference, source, manifest = prepared
    assert read_annotation_inputs(manifest, reference=reference, **filters)[0].source_vcf == source


@pytest.mark.parametrize('filters', [
    {'sample_id': 'MISSING'}, {'variant_category': 'tr'},
    {'sample_id': 'S1', 'variant_category': 'sv'},
    {'sample_id': ''}, {'sample_id': 42}, {'variant_category': 'UNKNOWN'},
])
def test_empty_or_invalid_selection_is_not_success(prepared, filters):
    reference, _, manifest = prepared
    with pytest.raises(InputValidationError):
        read_annotation_inputs(manifest, reference=reference, **filters)


def test_empty_plan_is_rejected_at_metadata_stage(tmp_path):
    manifest = tmp_path / 'empty.tsv'
    manifest.write_text('sample\tvariant_category\tsource_vcf\tsource_tool\tsource_variant_ids\n', encoding='utf-8')
    with pytest.raises(InputValidationError, match='at least one row'):
        read_annotation_input_rows(manifest)


def test_duplicate_unselected_rows_are_still_invalid_plan(prepared):
    reference, _, manifest = prepared
    with manifest.open('a', encoding='utf-8') as handle:
        handle.write('S2\tsv\tanother.future.vcf\tjasmine\tv3\n')
    with pytest.raises(InputValidationError, match='duplicates'):
        read_annotation_inputs(manifest, reference=reference, sample_id='S1', variant_category='small')


@pytest.mark.parametrize('row', [
    'S2\tsv\tfuture.vcf\tjasmine\n',
    'S2\tsv\tfuture.vcf\tjasmine\tv2\textra\n',
])
def test_malformed_unselected_metadata_is_not_silently_ignored(prepared, row):
    reference, _, manifest = prepared
    with manifest.open('a', encoding='utf-8') as handle:
        handle.write(row)
    with pytest.raises(InputValidationError, match='column count'):
        read_annotation_inputs(manifest, reference=reference, sample_id='S1', variant_category='small')


def test_default_all_inputs_interface_remains_available(prepared):
    reference, _, manifest = prepared
    ready_only = manifest.with_name('ready.tsv')
    ready_only.write_text('\n'.join(manifest.read_text(encoding='utf-8').splitlines()[:2]) + '\n', encoding='utf-8')
    assert len(read_annotation_inputs(ready_only, reference=reference)) == 1


def test_workflow_bridge_consumes_only_current_ready_row(tmp_path, prepared, monkeypatch):
    reference, source, manifest = prepared
    original = source.read_bytes()
    cache = tmp_path / 'cache'
    cache.mkdir()
    seen = []

    class FakeVepWrapper:
        def __init__(self, **kwargs):
            pass

        def run(self, request, **kwargs):
            seen.append(request.input)
            request.output_tsv.parent.mkdir(parents=True, exist_ok=True)
            request.output_tsv.write_text('#Uploaded_variation\tConsequence\nv1\tfeature\n', encoding='utf-8')
            artifact = AnnotationArtifact(request.input, AnnotationSource.VEP,
                                          request.output_tsv, 'tsv', '116.3', ())
            return AnnotationResult(request.input, AnnotationSource.VEP,
                                    AnnotationRunStatus.COMPLETED, ('fake-vep',),
                                    '116.3', artifact=artifact)

    monkeypatch.setattr(hifivar.vep, 'VepWrapper', FakeVepWrapper)
    monkeypatch.setattr(AnalysisContext, 'from_config', lambda config: SimpleNamespace(reference=reference))
    output = tmp_path / 'annotation-result'
    snakemake = SimpleNamespace(
        config={'annotation': {
            'vep_enabled': True, 'vep_executable': 'fake-vep',
            'annovar_executable': 'unused-annovar', 'vep_cache_directory': str(cache),
            'vep_cache_version': '116', 'vep_species': 'homo_sapiens', 'vep_assembly': 'GRCh38',
        }},
        input=SimpleNamespace(manifest=str(manifest)),
        output=SimpleNamespace(result=str(output)),
        wildcards=SimpleNamespace(sample='S1', category='small'),
    )
    script = Path(__file__).resolve().parents[2] / 'workflow/scripts/run_annotation.py'
    runpy.run_path(str(script), init_globals={'snakemake': snakemake})
    assert len(seen) == 1 and seen[0].sample_id == 'S1'
    assert (output / 'phase11.provenance.json').is_file()
    assert source.read_bytes() == original
    assert not (manifest.parent / 'future.S2.sv.vcf').exists()
