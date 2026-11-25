"""Deterministic per-sample Phase 9 source selection, not cross-sample clustering."""
from __future__ import annotations
from pathlib import Path
from hifivar.exceptions import ConfigurationError

READ_CALLERS = ('sawfish', 'sniffles2', 'pbsv', 'cutesv')
ASSEMBLY_CALLERS = ('pav', 'svim_asm')
CALLER_ORDER = READ_CALLERS + ASSEMBLY_CALLERS


def generated_callers(config: dict) -> tuple[str, ...]:
    """Select only explicitly enabled producers in their established order."""
    return tuple(name for group, names in (('sv', READ_CALLERS), ('assembly_sv', ASSEMBLY_CALLERS))
                 for name in names if config.get(group, {}).get('enabled') is True
                 and config[group].get(name, {}).get('enabled') is True)


def harmonization_sources(config: dict, sample: str, run_root: Path) -> tuple[dict, ...]:
    """Declare exact source roles without executing tools or guessing paths.

    Generated callers precede external callers. Both use CALLER_ORDER; YAML
    mapping order cannot change SUPP_VEC interpretation. External source files
    retain the existing exact sample/caller filename contract.
    """
    section = config.get('sv', {}).get('harmonization', {})
    external_all = section.get('input_vcfs', {})
    if not isinstance(external_all, dict):
        raise ConfigurationError('sv.harmonization.input_vcfs must be a mapping.')
    external = external_all.get(sample, {})
    if not isinstance(external, dict) or set(external) - set(CALLER_ORDER):
        raise ConfigurationError('External harmonization source has an unsupported caller or mapping.')
    generated = generated_callers(config)
    if set(generated) & set(external):
        raise ConfigurationError('A harmonization caller cannot be both generated and external.')
    result = []
    for caller in (*generated, *(name for name in CALLER_ORDER if name in external)):
        suffix = '.assembly.sv.vcf.gz' if caller in ASSEMBLY_CALLERS else '.sv.vcf.gz'
        expected = f'{sample}.{caller}{suffix}'
        origin = 'EXTERNAL' if caller in external else 'GENERATED'
        if origin == 'EXTERNAL':
            value = external[caller]
            if not isinstance(value, str) or not value.strip():
                raise ConfigurationError('External harmonization VCF must be a nonempty explicit path.')
            path = Path(value)
            if path.name != expected:
                raise ConfigurationError(f'External harmonization VCF must be named {expected!r}.')
        else:
            directory = run_root / 'results' / ('assembly_sv' if caller in ASSEMBLY_CALLERS else 'sv')
            path = directory / sample / expected if caller in ASSEMBLY_CALLERS else directory / expected
        if not path.is_absolute() or '..' in path.parts:
            raise ConfigurationError('Harmonization sources require normalized absolute paths.')
        result.append({'caller': caller, 'origin': origin, 'vcf': path,
                       'index': Path(str(path) + '.tbi'),
                       'track': ('assembly_sv.' if caller in ASSEMBLY_CALLERS else 'sv.') + caller})
    if not result:
        raise ConfigurationError('Harmonization requires generated or external evidence for every sample.')
    return tuple(result)
