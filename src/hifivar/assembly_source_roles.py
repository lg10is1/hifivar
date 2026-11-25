"""Source assembly-input roles bound to files and the current assembly receipt.

These are workflow input declarations, never per-variant haplotype support or
biological phasing accuracy. Unbound external files retain unknown roles.
"""
from __future__ import annotations
from collections.abc import Sequence
from pathlib import Path
from hifivar.assembly import HaplotypeAssemblyArtifact
from hifivar.ready_assemblies import ReadyAssemblyArtifact
from hifivar.exceptions import InputValidationError

_ROLES = ('haplotype1', 'haplotype2')
_ASSEMBLY = {'pav', 'svim_asm'}
_READ = {'sawfish', 'sniffles2', 'pbsv', 'cutesv'}
_SCOPE = 'EXPLICIT_ASSEMBLY_REQUEST_ROLES_BOUND_TO_CURRENT_ASSEMBLY_OUTPUTS'


def _owner(context_path: Path) -> tuple[Path, str]:
    path = Path(context_path).absolute()
    if path.name != 'execution_context.json' or path.parent.parent.name != 'attempts':
        raise InputValidationError('Assembly role binding requires an owned execution context.')
    return path.parent.parent.parent, path.parent.name


def _snapshot(sample_id: str, roles: Sequence[tuple[str, Path]], context_path: Path) -> dict:
    from hifivar.attempt_receipts import artifact_identity, read_step_receipt, receipt_name
    root, attempt = _owner(context_path)
    from hifivar.attempt_receipts import _context
    context, identity, _ = _context(context_path, live=False)
    names = [role for role, _ in roles]
    if (not 1 <= len(names) <= 2 or any(not isinstance(x, str) or x not in _ROLES for x in names)
            or len(set(names)) != len(names)):
        raise InputValidationError('Assembly roles must explicitly select unique haplotype1/haplotype2 inputs.')
    if identity.get('config', {}).get('assembly_sv', {}).get('input_manifest'):
        from hifivar.assembly_reuse import read_ready_handoff
        _, items = read_ready_handoff(root, context_path.parent, context.get('ready_assembly_handoff'),
                                     verify_inputs=True, sample_id=sample_id)
        selected = {item.role.value: item for item in items}
        if set(names) != set(selected) or any(Path(path).resolve() != selected[role].path for role, path in roles):
            raise InputValidationError('Actual ready assembly roles/paths differ from the frozen handoff.')
        return {'schema_version': 1, 'sample_id': sample_id,
                'scope': 'EXPLICIT_READY_ASSEMBLY_INPUTS_BOUND_TO_THIS_ATTEMPT_VALIDATION_ONLY',
                'phase_accuracy': 'NOT_ASSESSED',
                'assemblies': [{'role': role, 'fasta': artifact_identity(Path(path).absolute())} for role, path in roles],
                'input_validation': context['ready_assembly_handoff'], 'new_hifiasm_execution': False}
    producer = read_step_receipt(root, attempt, 'assembly', (sample_id,), verify_outputs=False)
    if producer is None:
        raise InputValidationError('Assembly role inputs lack a current-attempt assembly receipt.')
    entries, seen = [], set()
    for role, path in roles:
        info = artifact_identity(Path(path).absolute())
        if producer['outputs'].get(role + '_fasta') != info:
            raise InputValidationError('Declared assembly role file differs from its current producer output.')
        stat = Path(info['path']).stat()
        physical = (stat.st_dev, stat.st_ino)
        if physical in seen:
            raise InputValidationError('Different assembly roles cannot alias one physical FASTA.')
        seen.add(physical)
        entries.append({'role': role, 'fasta': info})
    receipt = root / 'attempts' / attempt / 'receipts' / receipt_name('assembly', (sample_id,))
    return {'schema_version': 1, 'sample_id': sample_id, 'scope': _SCOPE,
            'phase_accuracy': 'NOT_ASSESSED', 'assemblies': entries,
            'assembly_receipt': artifact_identity(receipt)}


def assembly_role_binding(sample_id: str, assemblies: Sequence[HaplotypeAssemblyArtifact | ReadyAssemblyArtifact],
                          context_path: Path | None) -> dict | None:
    """Snapshot actual request inputs; standalone wrappers do not invent receipts."""
    if context_path is None:
        return None
    if (not isinstance(assemblies, Sequence) or isinstance(assemblies, (str, bytes))
            or any(not isinstance(item, (HaplotypeAssemblyArtifact, ReadyAssemblyArtifact)) for item in assemblies)):
        raise InputValidationError('Assembly role binding requires explicit assembly artifact objects.')
    if any(item.sample_id != sample_id for item in assemblies):
        raise InputValidationError('Assembly role sample differs from the request sample.')
    from hifivar.attempt_receipts import _context
    _, identity, _ = _context(context_path, live=False)
    ready_mode = bool(identity.get('config', {}).get('assembly_sv', {}).get('input_manifest'))
    if any(isinstance(item, ReadyAssemblyArtifact) != ready_mode for item in assemblies):
        raise InputValidationError('Assembly artifact origin differs from the explicit frozen input mode.')
    for item in assemblies:
        if isinstance(item, ReadyAssemblyArtifact):
            item.verify()
    return _snapshot(sample_id, tuple((item.role.value, item.path) for item in assemblies), context_path)


def _verify_binding(value: object, sample_id: str, context_path: Path) -> dict:
    if not isinstance(value, dict) or type(value.get('schema_version')) is not int or value.get('schema_version') != 1 or value.get('sample_id') != sample_id:
        raise InputValidationError('Assembly-SV receipt lacks an explicit assembly role binding.')
    entries = value.get('assemblies')
    if not isinstance(entries, list) or any(not isinstance(item, dict) or set(item) != {'role', 'fasta'}
            or not isinstance(item['fasta'], dict) or not isinstance(item['fasta'].get('path'), str) for item in entries):
        raise InputValidationError('Assembly role binding contains malformed input identities.')
    current = _snapshot(sample_id, tuple((item['role'], Path(item['fasta']['path'])) for item in entries), context_path)
    if current != value:
        raise InputValidationError('Assembly input roles, bytes or upstream receipt changed.')
    return current


def source_haplotype_binding(sample_id: str, caller: str, *, external: bool,
                             context_path: Path | None) -> dict:
    """Resolve declared source input roles without inspecting filenames or GT."""
    if caller in _READ:
        return {'roles': [], 'status': 'NOT_APPLICABLE_READ_SOURCE', 'binding': None}
    if caller not in _ASSEMBLY:
        raise InputValidationError('Unsupported caller for assembly source roles.')
    if external or context_path is None:
        return {'roles': [], 'status': 'UNVERIFIED_EXTERNAL_ASSEMBLY_INPUT_ROLES' if external else
                'UNVERIFIED_ASSEMBLY_INPUT_ROLES_NO_ATTEMPT_BINDING', 'binding': None}
    from hifivar.attempt_receipts import read_step_receipt
    root, attempt = _owner(context_path)
    producer = read_step_receipt(root, attempt, 'assembly_sv.' + caller, (sample_id,), verify_outputs=False)
    if producer is None:
        raise InputValidationError('Assembly source has no current-attempt assembly-SV receipt.')
    value = _verify_binding(producer.get('metadata', {}).get('assembly_role_binding'), sample_id, context_path)
    return {'roles': [item['role'] for item in value['assemblies']],
            'status': ('CURRENT_ATTEMPT_BOUND_READY_SOURCE_INPUT_ROLES_NOT_VARIANT_SUPPORT'
                       if value.get('new_hifiasm_execution') is False else
                       'CURRENT_ATTEMPT_BOUND_SOURCE_ASSEMBLY_INPUT_ROLES_NOT_VARIANT_SUPPORT'),
            'binding': value}
