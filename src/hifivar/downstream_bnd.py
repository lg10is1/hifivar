"""Bind an explicitly selected derived role to a same-attempt BND producer.

This supplements, never replaces, the native harmonized VCF output role.
Receipt readback and bundle round-trip checks are not biological validation.
"""
from __future__ import annotations

from pathlib import Path
from typing import Mapping

from hifivar.breakends import BndPolicy
from hifivar.exceptions import InputValidationError
from hifivar.junction_consumers import validate_bnd_consumer_bundle
from hifivar.reference import ReferenceGenome


def generated_bnd_roles(receipt: Mapping, reference: ReferenceGenome, sample: str) -> dict:
    """Verify exact producer outputs, fixed policy, copied native bytes and ledger.

    The caller must already have read a receipt bound to the selected attempt.
    Every bundle member must be an explicitly sealed output, not merely a file
    next to a producer result. No files are created, normalized or indexed.
    """
    from hifivar.attempt_receipts import artifact_identity
    if receipt.get('track') != 'sv.harmonization' or receipt.get('samples') != [sample]:
        raise InputValidationError('Derived BND source needs the selected sample harmonization receipt.')
    outputs = receipt.get('outputs', {})
    required = {'vcf', 'index', 'evidence', 'bnd_manifest', 'bnd_vcf', 'bnd_index', 'bnd_ledger',
                'bnd_junctions', 'bnd_native_result', 'bnd_native_vcf', 'bnd_native_index',
                'bnd_native_raw', 'bnd_native_evidence'}
    if not isinstance(outputs, dict) or not required.issubset(outputs):
        raise InputValidationError('Generated BND receipt lacks its complete explicitly sealed bundle roles.')
    policy = BndPolicy.from_config(receipt.get('metadata', {}).get('bnd_representation', {}))
    if policy is None:
        raise InputValidationError('Generated BND receipt lacks an enabled fixed representation policy.')
    manifest = Path(outputs['bnd_manifest']['path'])
    source = Path(outputs['bnd_vcf']['path'])
    binding = validate_bnd_consumer_bundle(manifest, source, reference, sample)
    if binding['policy'] != policy.to_dict():
        raise InputValidationError('Generated BND bundle policy differs from its producer receipt.')
    roles = {
        'bnd_manifest': manifest, 'bnd_vcf': source, 'bnd_index': Path(str(source)+'.tbi'),
        'bnd_ledger': manifest.parent/'records.jsonl', 'bnd_junctions': manifest.parent/'junctions.tsv',
        'bnd_native_result': manifest.parent/'native_result.json',
        'bnd_native_vcf': manifest.parent/'native'/f'{sample}.harmonized.sv.vcf.gz',
        'bnd_native_index': manifest.parent/'native'/f'{sample}.harmonized.sv.vcf.gz.tbi',
        'bnd_native_raw': manifest.parent/'native'/f'{sample}.jasmine.raw.vcf',
        'bnd_native_evidence': manifest.parent/'native'/f'{sample}.sv.evidence.tsv',
    }
    for role, path in roles.items():
        if outputs[role] != artifact_identity(path):
            raise InputValidationError('Generated BND bundle member differs from its sealed output role.')
    for native_role, copied_role in (('vcf', 'bnd_native_vcf'), ('index', 'bnd_native_index'),
                                      ('evidence', 'bnd_native_evidence')):
        actual = artifact_identity(outputs[native_role]['path'])
        if outputs[native_role] != actual or any(actual[k] != outputs[copied_role][k] for k in ('bytes', 'sha256')):
            raise InputValidationError('Derived BND native copy differs from its canonical producer output.')
    return {'vcf': outputs['bnd_vcf'], 'index': outputs['bnd_index'],
            'bnd_bundle_manifest': outputs['bnd_manifest'], 'bundle_binding': binding}


__all__ = ['generated_bnd_roles']
