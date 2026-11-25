"""Bounded IGV session identity checks, independent of human interpretation."""
from __future__ import annotations

import hashlib
from pathlib import Path
import re
from urllib.parse import unquote, urlparse
import xml.etree.ElementTree as ET

from hifivar.exceptions import OutputValidationError


def _local_path(value: str | None, parent: Path) -> Path:
    if not isinstance(value, str) or not value.strip() or any(ord(c) < 32 for c in value):
        raise OutputValidationError('IGV session has a missing or invalid local resource path.')
    if value.startswith('file:'):
        parsed = urlparse(value)
        if parsed.netloc not in {'', 'localhost'}:
            raise OutputValidationError('IGV session remote file hosts are unsupported.')
        value = unquote(parsed.path)
        if re.match(r'^/[A-Za-z]:/', value):
            value = value[1:]
    elif re.match(r'^[a-zA-Z][a-zA-Z0-9+.-]*://', value):
        raise OutputValidationError('IGV session resources must be local files.')
    path = Path(value).expanduser()
    return (parent / path if not path.is_absolute() else path).resolve()


def validate_igv_session(path: str | Path, *, reference: Path, alignment: Path,
                         source_vcf: Path, contig: str, start: int, end: int,
                         tolerance_bp: int = 2) -> dict[str, object]:
    """Require the saved native session to bind reference, tracks and locus.

    IGV's text locus can round the requested viewport by a few bases. The
    explicit tolerance is recorded. Track presence does not establish data
    coverage, genotype, rendering quality or independent variant validity.
    """
    session = Path(path)
    if (not isinstance(contig, str) or not contig or type(start) is not int
            or type(end) is not int or not 1 <= start <= end
            or type(tolerance_bp) is not int or not 0 <= tolerance_bp <= 2):
        raise OutputValidationError('Invalid expected IGV session locus or tolerance.')
    try:
        if not session.is_file() or not 0 < session.stat().st_size <= 2 * 1024 * 1024:
            raise OutputValidationError('IGV session is missing/empty or larger than 2 MiB.')
        with session.open('rb') as handle:
            raw = handle.read(2 * 1024 * 1024 + 1)
        if len(raw) > 2 * 1024 * 1024:
            raise OutputValidationError('IGV session grew beyond the 2 MiB bounded limit.')
        text = raw.decode('utf-8-sig')
        if '\x00' in text or re.search(r'<!\s*(?:DOCTYPE|ENTITY)\b', text, re.IGNORECASE):
            raise OutputValidationError('IGV session DTD/entities or non-UTF8 XML are forbidden.')
        root = ET.fromstring(text)
    except (OSError, UnicodeError, ET.ParseError) as error:
        raise OutputValidationError(f'Unable to parse bounded IGV session: {error}') from error
    if root.tag != 'Session':
        raise OutputValidationError('IGV session root must be Session.')
    expected_reference = Path(reference).resolve()
    actual_reference = _local_path(root.get('genome'), session.parent)
    if actual_reference != expected_reference:
        raise OutputValidationError('IGV session reference differs from the selected reference.')
    resources = {_local_path(item.get('path'), session.parent)
                 for item in root.findall('./Resources/Resource')}
    expected_resources = {Path(alignment).resolve(), Path(source_vcf).resolve()}
    if not expected_resources <= resources:
        raise OutputValidationError('IGV session lacks the selected BAM/CRAM or VCF resource.')
    locus_text = root.get('locus', '')
    locus = re.fullmatch(r'([^:\s]+):([0-9,]+)-([0-9,]+)', locus_text)
    if locus is None:
        raise OutputValidationError('IGV session locus is absent, whole-genome or malformed.')
    observed_start = int(locus[2].replace(',', ''))
    observed_end = int(locus[3].replace(',', ''))
    if (locus[1] != contig or not 1 <= observed_start <= observed_end
            or abs(observed_start - start) > tolerance_bp
            or abs(observed_end - end) > tolerance_bp):
        raise OutputValidationError('IGV session locus differs from the expected screenshot window.')
    return {'session_path': str(session), 'sha256': hashlib.sha256(raw).hexdigest(),
            'reference_path_match': True, 'selected_resources_present': True,
            'expected_locus': f'{contig}:{start}-{end}', 'observed_locus': locus_text,
            'coordinate_tolerance_bp': tolerance_bp,
            'status': 'REFERENCE_RESOURCES_LOCUS_MATCH', 'human_review': 'NOT_REVIEWED',
            'data_coverage': 'NOT_ASSESSED', 'variant_accuracy': 'NOT_ASSESSED'}


__all__ = ['validate_igv_session']
