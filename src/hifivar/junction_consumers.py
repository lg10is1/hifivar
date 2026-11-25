"""Explicit selected-record BND consumers; anchors are never event spans.

Reads native VCFs or an explicitly bound, reversible BND bundle. Only selected
records and their declared mate IDs are materialized. Missing links remain QC,
not invented partners or biological validity decisions.
"""
from __future__ import annotations

from dataclasses import dataclass
import gzip
import json
from pathlib import Path
from typing import Mapping, Sequence

from hifivar.breakends import BndPolicy, BreakendJunction, parse_breakend
from hifivar.exceptions import InputValidationError
from hifivar.reference import ReferenceGenome
from hifivar.validation import compute_sha256


def _identity(path: Path) -> dict:
    if not path.is_file() or path.is_symlink():
        raise InputValidationError('Junction consumers require regular non-symlink source files.')
    before = path.stat()
    digest = compute_sha256(path)
    after = path.stat()
    if (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns) != (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns):
        raise InputValidationError('Junction input changed while hashing.')
    return {'path': str(path.resolve()), 'sha256': digest, 'bytes': after.st_size}


def validate_bnd_consumer_bundle(manifest_path: Path, source_vcf: Path,
                                 reference: ReferenceGenome, sample_id: str) -> dict:
    """Bind every mandatory bundle file, exact reference and full record round trip.

    This proves internal bytes/roles, not authenticity against malicious edits or
    biological correctness. It does not rerun tools, index or change any input.
    """
    from hifivar.bnd_records import verify_record_round_trip
    from hifivar.annotation import audit_annotation_vcf_identity
    path = Path(manifest_path).absolute()
    manifest_identity = _identity(path)
    if manifest_identity['bytes'] > 1024 * 1024:
        raise InputValidationError('BND consumer manifest exceeds 1 MiB.')
    try:
        data = json.loads(path.read_text(encoding='utf-8'))
        if type(data.get('schema_version')) is not int or data['schema_version'] != 1 or data['status'] != 'completed':
            raise InputValidationError('BND consumer requires a completed supported bundle.')
        policy = BndPolicy.from_config(data['policy'])
        if policy is None or data.get('native_output_role_unchanged') is not True:
            raise InputValidationError('BND bundle policy/native role is invalid.')
        artifacts, copies = data['artifact_identities'], data['native_binary_copies']
        if not isinstance(artifacts, dict) or set(copies) != {'vcf', 'index', 'raw', 'evidence'}:
            raise InputValidationError('BND bundle mandatory native roles are missing.')
        root = path.parent.resolve()
        identities = {}
        for name, expected in artifacts.items():
            part = Path(name)
            if part.is_absolute() or '..' in part.parts or not part.parts:
                raise InputValidationError('BND artifact path escapes its bundle.')
            artifact = path.parent / part
            if not artifact.resolve().is_relative_to(root):
                raise InputValidationError('BND artifact leaves its bundle.')
            actual = _identity(artifact)
            if (actual['sha256'], actual['bytes']) != (expected['sha256'], expected['bytes']):
                raise InputValidationError('BND bundle artifact bytes changed.')
            identities[name] = actual
        derived_name = f'{sample_id}.bndrepr-v1.sv.vcf.gz'
        required = {derived_name, derived_name + '.tbi', 'records.jsonl', 'junctions.tsv', 'native_result.json', *copies.values()}
        if not required.issubset(identities) or Path(identities[derived_name]['path']) != Path(source_vcf).resolve():
            raise InputValidationError('Selected derived VCF is not the exact bundle output role.')
        if data['native_inputs_before'] != data['native_inputs_after']:
            raise InputValidationError('BND bundle upstream immutability was not established.')
        for ref in (reference.fasta, reference.fai):
            actual = _identity(ref)
            expected = data['native_inputs_before'].get(actual['path'])
            if not expected or (actual['sha256'], actual['bytes']) != (expected['sha256'], expected['bytes']):
                raise InputValidationError('BND consumer reference differs from bundle reference.')
        native = path.parent / copies['vcf']
        audit_annotation_vcf_identity(native, sample_id, reference)
        round_trip = verify_record_round_trip(native, source_vcf, path.parent / 'records.jsonl')
        if _identity(path) != manifest_identity or any(_identity(Path(info['path'])) != info for info in identities.values()):
            raise InputValidationError('BND bundle changed during consumer validation.')
        return {'role': 'EXPLICIT_VERIFIED_BND_DERIVED_WITH_NATIVE_LEDGER',
                'manifest': manifest_identity, 'artifacts': identities,
                'native_vcf': identities[copies['vcf']], 'record_round_trip': round_trip,
                'policy': policy.to_dict(), 'biological_validity': 'NOT_ASSESSED'}
    except (OSError, UnicodeError, ValueError, KeyError, TypeError, AttributeError) as error:
        raise InputValidationError(f'Cannot validate BND consumer bundle: {error}') from error


def _info(text: str) -> dict[str, str | None]:
    result = {}
    for token in (() if text == '.' else text.split(';')):
        key, separator, value = token.partition('=')
        if not key or key in result or (separator and not value):
            raise InputValidationError('Selected junction record has malformed or duplicate INFO.')
        result[key] = value if separator else None
    return result


def _mate_ids(info: Mapping[str, str | None]) -> tuple[str, ...]:
    if 'MATEID' not in info or info['MATEID'] == '.':
        return ()
    value = info['MATEID']
    if not isinstance(value, str) or any(not part for part in value.split(',')):
        raise InputValidationError('Selected BND has an invalid MATEID declaration.')
    return tuple(value.split(','))


def _record_junction(fields: list[str], info: Mapping[str, str | None],
                     lengths: Mapping[str, int]) -> BreakendJunction | None:
    alt = fields[4]
    is_bnd = (info.get('SVTYPE') == 'BND' or '[' in alt or ']' in alt
              or (alt != '.' and (alt.startswith('.') or alt.endswith('.'))))
    if not is_bnd:
        return None  # Ordinary multi-ALT SV/TR records retain their existing locus contract.
    return parse_breakend(fields[0], int(fields[1]), fields[3], alt, lengths, svtype=info.get('SVTYPE'))


def _lines(path: Path):
    opener = gzip.open if path.name.lower().endswith('.gz') else open
    with opener(path, 'rb') as stream:
        yield from stream


@dataclass(frozen=True, slots=True)
class SelectedJunctionRecord:
    source_vcf: Path
    source_sha256: str
    ordinal: int
    record_id: str
    record_sha256: str
    contig: str
    position: int
    end: int
    ref: str
    alt: str
    junction: BreakendJunction | None
    mate_qc: dict
    source_role: str = 'EXPLICIT_NATIVE_VCF'
    bundle_binding: dict | None = None

    def to_dict(self) -> dict:
        return {'source_vcf': str(self.source_vcf), 'source_sha256': self.source_sha256,
                'record_ordinal': self.ordinal, 'record_id': self.record_id,
                'source_record_key': f'sha256:{self.source_sha256}:record:{self.ordinal}',
                'selected_record_sha256': self.record_sha256, 'source_role': self.source_role,
                'contig': self.contig, 'position': self.position, 'ref': self.ref, 'alt': self.alt,
                'junction': self.junction.to_dict() if self.junction else None,
                'mate_qc': self.mate_qc, 'bundle_binding': self.bundle_binding,
                'coordinate_semantics': 'VCF_ANCHORS_NOT_EVENT_SPAN' if self.junction else 'VCF_DECLARED_LOCUS',
                'annotation_is_confidence': False, 'phase_accuracy': 'NOT_ASSESSED'}


def read_selected_junction_records(source_vcf: Path, identifiers: Sequence[str], *,
                                   reference: ReferenceGenome, sample_id: str,
                                   bundle_manifest: Path | None = None) -> tuple[SelectedJunctionRecord, ...]:
    """Stream exact selected IDs and bounded explicit mate links, without guessing.

    Repeated selected IDs are ambiguous and rejected. Paired ALT coordinates are
    retained even when MATEID is absent/unresolved; single breakends have no mate.
    Virtual 0/N+1 anchors are metadata, never clamped into reference intervals.
    """
    from hifivar.annotation import audit_annotation_vcf_identity
    import hashlib
    selected = tuple(identifiers)
    if not selected or any(not isinstance(x, str) or not x or x == '.' for x in selected) or len(set(selected)) != len(selected):
        raise InputValidationError('Junction selection requires unique, nonmissing explicit record IDs.')
    source = Path(source_vcf).absolute()
    audit_annotation_vcf_identity(source, sample_id, reference)
    before = _identity(source)
    lengths = {c.name: c.length for c in reference.contigs}
    binding = validate_bnd_consumer_bundle(bundle_manifest, source, reference, sample_id) if bundle_manifest is not None else None
    found, derived_header, ordinal = {}, False, 0
    try:
        for raw in _lines(source):
            if raw.startswith(b'#'):
                derived_header |= raw.startswith(b'##HIFIVAR_BND_')
                continue
            if not raw.strip():
                continue
            ordinal += 1
            fields = raw.decode('utf-8').rstrip('\r\n').split('\t')
            if len(fields) < 8 or fields[2] not in selected:
                continue
            identifier = fields[2]
            if identifier in found:
                raise InputValidationError('Selected source record ID is ambiguous: duplicate records.')
            info = _info(fields[7])
            if {'HIFIVAR_BND_REPR', 'HIFIVAR_BND_LEDGER'}.intersection(info):
                derived_header = True
            position = int(fields[1])
            junction = _record_junction(fields, info, lengths)
            end = position if junction else int(info.get('END', str(position)))
            if not junction and (position < 1 or end < position or end > lengths[fields[0]]):
                raise InputValidationError('Selected non-BND locus is outside ordinary reference bounds.')
            found[identifier] = (ordinal, fields, info, junction, end, hashlib.sha256(raw).hexdigest())
        if set(found) != set(selected):
            raise InputValidationError('Selected source variant IDs are absent from VCF.')
        if derived_header != (binding is not None):
            raise InputValidationError('Derived BND input requires its explicit verified bundle; native inputs must not attach a derived bundle.')
        mate_names = {name for _, _, info, junction, _, _ in found.values() if junction for name in _mate_ids(info)}
        mates, counts = {}, {name: 0 for name in mate_names}
        for raw in _lines(source) if mate_names else ():
            if raw.startswith(b'#') or not raw.strip():
                continue
            fields = raw.decode('utf-8').rstrip('\r\n').split('\t')
            if len(fields) >= 8 and fields[2] in mate_names:
                name = fields[2]; counts[name] += 1
                if counts[name] == 1:
                    mates[name] = fields
        result = []
        for identifier in selected:
            number, fields, info, junction, end, digest = found[identifier]
            links = []
            for name in _mate_ids(info) if junction else ():
                count = counts[name]
                link = {'mate_id': name, 'matching_record_count': count,
                        'id_resolution': 'UNRESOLVED' if count == 0 else 'AMBIGUOUS_DUPLICATE_ID' if count > 1 else 'RESOLVED_ID_ONLY'}
                if count == 1:
                    target = mates[name]
                    other = _record_junction(target, _info(target[7]), lengths)
                    link.update(target_is_breakend=other is not None,
                                exact_alt_mate_coordinate=(junction.mate_contig, junction.mate_position) == (target[0], int(target[1])),
                                reciprocal_id_link=identifier in _mate_ids(_info(target[7])),
                                self_link=name == identifier)
                    if other:
                        link['orientation_consistent'] = junction.local_side == other.mate_side and junction.mate_side == other.local_side
                links.append(link)
            qc = {'id_resolution': 'DECLARED_LINK_QC' if links else 'NO_MATEID', 'links': links,
                  'coordinate_inequality_is_invalid': False, 'biological_validity': 'NOT_ASSESSED'}
            result.append(SelectedJunctionRecord(source, before['sha256'], number, identifier, digest,
                fields[0], int(fields[1]), end, fields[3], fields[4], junction, qc,
                binding['role'] if binding else 'EXPLICIT_NATIVE_VCF', binding))
        if _identity(source) != before:
            raise InputValidationError('Selected junction source changed during reading.')
        return tuple(result)
    except (OSError, EOFError, UnicodeError, ValueError, KeyError, TypeError) as error:
        raise InputValidationError(f'Cannot read selected junction source: {error}') from error


def revalidate_junction_record(record: SelectedJunctionRecord) -> None:
    if _identity(record.source_vcf)['sha256'] != record.source_sha256:
        raise InputValidationError('Selected BND source bytes changed after selection.')
    if record.bundle_binding:
        for info in (record.bundle_binding['manifest'], *record.bundle_binding['artifacts'].values()):
            if _identity(Path(info['path'])) != info:
                raise InputValidationError('Selected BND bundle changed after selection.')


__all__ = ['SelectedJunctionRecord', 'read_selected_junction_records', 'validate_bnd_consumer_bundle', 'revalidate_junction_record']
