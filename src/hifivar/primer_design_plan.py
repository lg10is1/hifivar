"""Strict local candidate-primer plans bound to selected VCF alleles.

No calling, online sequence upload or experimental validation is performed.
FASTA access uses FAI byte offsets, never a whole-genome in-memory string.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import re
from pathlib import Path

import yaml

from hifivar.exceptions import InputValidationError, OutputValidationError
from hifivar.primer_design import PrimerDesignSettings, VariantTemplate
from hifivar.primer3 import Primer3Wrapper
from hifivar.review_cli import _output_path, _reserve, _write_text
from hifivar.serialization import write_json_atomic

_DNA = re.compile(r'[ACGT]+\Z')
_ID = re.compile(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,99}\Z')


class _Loader(yaml.SafeLoader):
    pass


def _mapping(loader, node, deep=False):
    pairs = loader.construct_pairs(node, deep=deep)
    result = {}
    for key, value in pairs:
        if not isinstance(key, str) or key in result:
            raise InputValidationError('Primer YAML keys must be unique strings.')
        result[key] = value
    return result


_Loader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, _mapping)


def _sha(path: Path) -> str:
    value = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            value.update(chunk)
    return value.hexdigest()


def _stat(path: Path) -> dict:
    value = path.stat()
    if not path.is_file():
        raise InputValidationError(f'Not a regular input file: {path}')
    return {'bytes': value.st_size, 'mtime_ns': value.st_mtime_ns}


def _file(value, base: Path) -> Path:
    if not isinstance(value, str) or not value.strip() or any(ord(x) < 32 for x in value):
        raise InputValidationError('Primer input paths must be nonempty single-line strings.')
    path = Path(value).expanduser()
    path = (base / path if not path.is_absolute() else path).resolve()
    if not path.is_file():
        raise InputValidationError(f'Primer input file does not exist: {path}')
    return path


def read_fai(path: Path) -> dict[str, tuple[int, int, int, int]]:
    """Read bounded conventional uncompressed FASTA offsets."""
    if path.stat().st_size > 8 * 1024 * 1024:
        raise InputValidationError('FASTA index exceeds the 8 MiB bounded plan limit.')
    result = {}
    for line in path.read_text(encoding='utf-8').splitlines():
        fields = line.split('\t')
        if len(fields) < 5 or fields[0] in result or not fields[0]:
            raise InputValidationError('Malformed or duplicate FASTA index entry.')
        try:
            length, offset, line_bases, line_width = map(int, fields[1:5])
        except ValueError as error:
            raise InputValidationError('Noninteger FASTA index value.') from error
        if length < 1 or offset < 0 or not 1 <= line_bases <= 100000000 or line_width-line_bases not in {0, 1, 2}:
            raise InputValidationError('Invalid FASTA index dimensions.')
        result[fields[0]] = (length, offset, line_bases, line_width)
    if not result:
        raise InputValidationError('Empty FASTA index.')
    return result


def fetch_sequence(fasta: Path, entry: tuple[int, int, int, int], start: int, end: int) -> str:
    """Fetch a bounded 1-based closed interval from an uncompressed indexed FASTA."""
    length, offset, line_bases, line_width = entry
    if type(start) is not int or type(end) is not int or not 1 <= start <= end <= length or end-start+1 > 100000:
        raise InputValidationError('Primer sequence interval must be valid and at most 100000 bp.')
    first = start - 1
    last = end - 1
    begin = offset + (first // line_bases) * line_width + first % line_bases
    finish = offset + (last // line_bases) * line_width + last % line_bases + 1
    if not 0 <= begin < finish <= fasta.stat().st_size or finish-begin > 300002:
        raise InputValidationError('FASTA index points outside the file or exceeds bounded read size.')
    with fasta.open('rb') as handle:
        handle.seek(begin)
        raw = handle.read(finish - begin)
    try:
        sequence = raw.replace(b'\r', b'').replace(b'\n', b'').decode('ascii').upper()
    except UnicodeError as error:
        raise InputValidationError('FASTA index does not address ASCII sequence.') from error
    if len(sequence) != end-start+1 or re.fullmatch('[ACGTN]+', sequence) is None:
        raise InputValidationError('FASTA sequence/index mismatch or unsupported sequence alphabet.')
    return sequence


def _selected_records(path: Path, selections: list[dict], lengths: dict) -> list[dict]:
    """Stream an explicit biallelic selection; retain exact raw-record identity."""
    wanted = {(x['sample'], x['contig'], x['position'], x['ref'], x['alt']): x for x in selections}
    if len(wanted) != len(selections):
        raise InputValidationError('Duplicate selected sample/allele in one VCF.')
    found = {}
    samples = None
    declarations = {}
    opener = gzip.open if str(path).lower().endswith('.gz') else open
    try:
        with opener(path, 'rt', encoding='utf-8', newline='') as handle:
            ordinal = 0
            for line in iter(lambda: handle.readline(2 * 1024 * 1024 + 1), ''):
                if len(line) > 2 * 1024 * 1024:
                    raise InputValidationError('VCF record exceeds the 2 MiB candidate limit.')
                if line.startswith('##contig='):
                    name = re.search(r'(?:<|,)ID=([^,>]+)', line)
                    size = re.search(r'(?:<|,)length=(\d+)', line, re.I)
                    if name and size:
                        if name[1] in declarations:
                            raise InputValidationError('Duplicate VCF contig declaration.')
                        declarations[name[1]] = int(size[1])
                    continue
                if line.startswith('#CHROM\t'):
                    if samples is not None:
                        raise InputValidationError('Duplicate VCF #CHROM header.')
                    samples = line.rstrip('\r\n').split('\t')[9:]
                    if not samples or len(samples) != len(set(samples)) or any(not x for x in samples):
                        raise InputValidationError('VCF requires unique nonempty sample columns.')
                    continue
                if line.startswith('#'):
                    continue
                if samples is None:
                    raise InputValidationError('VCF data precedes sample header.')
                ordinal += 1
                fields = line.rstrip('\r\n').split('\t')
                if len(fields) != 9 + len(samples):
                    raise InputValidationError('VCF column count differs from sample header.')
                try:
                    pos = int(fields[1])
                except ValueError as error:
                    raise InputValidationError('Invalid VCF position.') from error
                for sample in {x['sample'] for x in selections}:
                    key = (sample, fields[0], pos, fields[3], fields[4])
                    if key not in wanted:
                        continue
                    if key in found or sample not in samples:
                        raise InputValidationError('Ambiguous selected VCF record or missing sample.')
                    if declarations.get(fields[0]) != lengths.get(fields[0], (None,))[0]:
                        raise InputValidationError('Selected VCF contig length does not match reference FAI.')
                    fmt = fields[8].split(':')
                    if len(fmt) != len(set(fmt)) or 'GT' not in fmt:
                        raise InputValidationError('Selected VCF record requires an unambiguous GT field.')
                    data = fields[9 + samples.index(sample)].split(':')
                    gt_index = fmt.index('GT')
                    gt = data[gt_index] if gt_index < len(data) else '.'
                    if re.fullmatch(r'[01](?:[|/][01])*', gt) is None or '1' not in re.split(r'[|/]', gt):
                        raise InputValidationError('Selected sample must have a called ALT genotype; ambiguous/missing GT is deferred.')
                    found[key] = {'record_ordinal': ordinal, 'record_id': fields[2], 'genotype': gt,
                                  'record_sha256': hashlib.sha256(line.encode('utf-8')).hexdigest()}
    except (OSError, UnicodeError, EOFError) as error:
        raise InputValidationError(f'Unable to stream selected VCF: {error}') from error
    if len(found) != len(wanted):
        raise InputValidationError('One or more exact selected VCF alleles/sample identities were not found.')
    return [found[(x['sample'], x['contig'], x['position'], x['ref'], x['alt'])] for x in selections]


def run_primer_design_plan(plan_path: str | Path, output_directory: str | Path, *, dry_run: bool = False) -> dict:
    """Design candidate pairs locally; no off-target or wet-lab success is inferred."""
    plan = Path(plan_path).expanduser().resolve()
    output = _output_path(output_directory)
    if not plan.is_file() or not 0 < plan.stat().st_size <= 1024 * 1024:
        raise InputValidationError('Primer plan must be a nonempty UTF-8 YAML file <=1 MiB.')
    raw = plan.read_bytes()
    try:
        doc = yaml.load(raw.decode('utf-8-sig'), Loader=_Loader)
    except (yaml.YAMLError, UnicodeError, RecursionError) as error:
        raise InputValidationError('Invalid bounded primer YAML.') from error
    if not isinstance(doc, dict) or set(doc) - {'schema_version', 'reference', 'primer3_executable', 'settings', 'targets', 'flank_bp'}:
        raise InputValidationError('Unknown primer plan fields.')
    if type(doc.get('schema_version')) is not int or doc['schema_version'] != 1:
        raise InputValidationError('Primer plan schema_version must be integer 1.')
    reference = doc.get('reference')
    if not isinstance(reference, dict) or set(reference) != {'fasta', 'build'} or not isinstance(reference['build'], str) or not reference['build'].strip():
        raise InputValidationError('Reference requires explicit fasta and build.')
    fasta = _file(reference['fasta'], plan.parent)
    if fasta.suffix.lower() in {'.gz', '.bgz'}:
        raise InputValidationError('Primer first version requires an uncompressed indexed FASTA.')
    fai = _file(str(fasta) + '.fai', plan.parent)
    lengths = read_fai(fai)
    flank = doc.get('flank_bp', 1000)
    if type(flank) is not int or not 100 <= flank <= 20000:
        raise InputValidationError('flank_bp must be an integer in 100..20000.')
    values = doc.get('settings', {})
    if not isinstance(values, dict):
        raise InputValidationError('Primer settings must be a mapping.')
    try:
        settings = PrimerDesignSettings(**values)
    except TypeError as error:
        raise InputValidationError(f'Unknown primer settings: {error}') from error
    targets = doc.get('targets')
    if not isinstance(targets, list) or not 1 <= len(targets) <= 100:
        raise InputValidationError('Primer plan requires 1..100 explicit targets.')
    groups = {}
    ordered = []
    for item in targets:
        required = {'target_id', 'vcf', 'sample', 'contig', 'position', 'ref', 'alt'}
        if not isinstance(item, dict) or set(item) != required:
            raise InputValidationError('Primer target requires exactly target_id/vcf/sample/contig/position/ref/alt.')
        if not isinstance(item['target_id'], str) or _ID.fullmatch(item['target_id']) is None or item['target_id'] in {x['target_id'] for x in ordered}:
            raise InputValidationError('Primer target_id must be unique and ASCII-safe, at most100 characters.')
        if not isinstance(item['sample'], str) or not item['sample'] or any(x.isspace() for x in item['sample']):
            raise InputValidationError('Primer sample must be a nonempty single field.')
        if not isinstance(item['contig'], str) or item['contig'] not in lengths or type(item['position']) is not int or item['position'] < 1:
            raise InputValidationError('Invalid primer target contig/position.')
        if any(not isinstance(item[k], str) or _DNA.fullmatch(item[k]) is None for k in ('ref', 'alt')) or item['ref'] == item['alt']:
            raise InputValidationError('Primer first version requires explicit uppercase biallelic DNA; symbolic/complex alleles are deferred.')
        selected = dict(item)
        selected['vcf'] = _file(item['vcf'], plan.parent)
        ordered.append(selected)
        groups.setdefault(selected['vcf'], []).append(selected)
    inputs = [plan, fasta, fai, *groups]
    for resource in inputs:
        if output.resolve() == resource.resolve() or output.resolve() in resource.resolve().parents:
            raise OutputValidationError('Primer output overlaps a protected input.')
    before = {str(p): _stat(p) for p in inputs}
    identities = {str(p): _sha(p) for p in inputs if p != fasta}
    records = {}
    for vcf, items in groups.items():
        for item, record in zip(items, _selected_records(vcf, items, lengths), strict=True):
            records[item['target_id']] = record
    templates = []
    for item in ordered:
        length = lengths[item['contig']][0]
        start = max(1, item['position'] - flank)
        end = min(length, item['position'] + len(item['ref']) - 1 + flank)
        sequence = fetch_sequence(fasta, lengths[item['contig']], start, end)
        record = records[item['target_id']]
        templates.append(VariantTemplate(sample=item['sample'], contig=item['contig'], position=item['position'],
            ref=item['ref'], alt=item['alt'], sequence_start=start, reference_sequence=sequence,
            source_variant_id=record['record_id'], source_record_sha256=record['record_sha256']))
    executable = doc.get('primer3_executable', 'primer3_core')
    if not isinstance(executable, str) or not executable.strip():
        raise InputValidationError('primer3_executable must be a nonempty string.')
    if ('/' in executable or '\\' in executable) and not Path(executable).is_absolute():
        executable = str((plan.parent / executable).resolve())
    wrapper = Primer3Wrapper(executable=executable)
    base = {'status': 'PLANNED' if dry_run else 'CANDIDATE_DESIGN_COMPLETE', 'target_count': len(ordered),
            'reference_build': reference['build'], 'reference_fasta': str(fasta),
            'reference_hash_policy': 'FAI_SHA256_AND_LOCAL_TEMPLATE_HASH_WITH_FASTA_STAT_NOT_WHOLE_FASTA_HASH',
            'source_sha256': identities, 'specificity_status': 'NOT_ASSESSED', 'experimental_status': 'NOT_RUN',
            'output_directory': str(output)}
    if dry_run:
        base['candidate_pair_count'] = None
        return base
    _reserve(output)
    results = []
    for item, template in zip(ordered, templates, strict=True):
        result = wrapper.design(template, output / item['target_id'], settings)
        result['target_id'] = item['target_id']
        result['source_vcf'] = str(item['vcf'])
        result['source_record'] = records[item['target_id']]
        results.append(result)
    if before != {str(p): _stat(p) for p in inputs} or identities != {str(p): _sha(p) for p in inputs if p != fasta}:
        raise InputValidationError('Primer inputs changed during candidate design; no completion receipt written.')
    rows = []
    products = []
    band_rows = []
    for result in results:
        for index, pair in enumerate(result['candidates'], 1):
            row = {'target_id': result['target_id'], 'candidate': index, 'assay': result['assay'],
                   **{k: v for k, v in pair.items() if not k.endswith('_product_sequence')},
                   'specificity_status': 'NOT_ASSESSED', 'experimental_status': 'NOT_RUN'}
            rows.append({k: json.dumps(v, ensure_ascii=False, sort_keys=True) if isinstance(v, (dict, list, tuple)) else v for k,v in row.items()})
            gt = result['source_record']['genotype']
            alleles = set(re.split(r'[|/]', gt))
            expected = {name: pair[name+'_product_bp'] for code,name in (('0','reference'),('1','alternate')) if code in alleles}
            unique_sizes = sorted(set(expected.values()))
            band_rows.append({'target_id':result['target_id'], 'candidate':index, 'source_genotype':gt,
                'ploidy_from_gt':len(re.split(r'[|/]',gt)), 'predicted_allele_products_bp':';'.join(f'{k}:{v}' for k,v in expected.items()),
                'distinct_product_lengths_bp':';'.join(map(str,unique_sizes)),
                'predicted_length_pattern':'ONE_LENGTH_CLASS' if len(unique_sizes)==1 else 'TWO_LENGTH_CLASSES',
                'gel_discrimination':'NO_FOR_SNV' if result['assay']=='PCR_SANGER' else 'NOT_ASSESSED',
                'actual_bands':'NOT_MEASURED', 'interpretation':'COMPUTATIONAL_PRODUCTS_NOT_OBSERVED_GEL;NO_INTENSITY_OR_ALLELE_DROPOUT_PREDICTION'})
            for allele in ('reference', 'alternate'):
                products.append(f">{result['target_id']}.pair{index}.{allele}\n{pair[allele+'_product_sequence']}\n")
    fields = list(dict.fromkeys(key for row in rows for key in row)) or ['target_id', 'candidate', 'assay']
    stream = io.StringIO(newline='')
    writer = csv.DictWriter(stream, fieldnames=fields, delimiter='\t', lineterminator='\n')
    writer.writeheader(); writer.writerows(rows)
    _write_text(output / 'primers.tsv', stream.getvalue())
    sequence_fields = ('target_id', 'candidate', 'assay', 'forward_sequence_5to3', 'reverse_sequence_5to3',
                       'reference_product_bp', 'alternate_product_bp', 'specificity_status', 'experimental_status')
    sequence_path = output / 'primer_sequences.tsv'
    sequence_stream = io.StringIO(newline='')
    sequence_writer = csv.DictWriter(sequence_stream, fieldnames=sequence_fields, delimiter='\t',
                                     lineterminator='\n', extrasaction='ignore')
    sequence_writer.writeheader()
    sequence_writer.writerows(rows)
    _write_text(sequence_path, sequence_stream.getvalue())
    _write_text(output / 'predicted_amplicons.fasta', ''.join(products))
    band_stream = io.StringIO(newline='')
    band_fields = list(band_rows[0]) if band_rows else ['target_id','candidate','source_genotype','actual_bands']
    band_writer = csv.DictWriter(band_stream,fieldnames=band_fields,delimiter='\t',lineterminator='\n')
    band_writer.writeheader(); band_writer.writerows(band_rows)
    _write_text(output / 'predicted_bands.tsv', band_stream.getvalue())
    base.update(candidate_pair_count=len(rows), targets=results,
                targets_with_requested_count=sum(x['requested_count_satisfied'] for x in results),
                manual_review_status='NOT_REVIEWED', primer_sequences_file=str(sequence_path))
    write_json_atomic(base, output / 'primer_design_summary.json')
    return base
