"""YAML configuration loading and validation for HiFiVar."""

from __future__ import annotations

from collections.abc import Iterator, Mapping
from copy import deepcopy
from dataclasses import dataclass, field
from pathlib import Path

import yaml

from hifivar.exceptions import ConfigurationError
from hifivar.logging_utils import get_logger, parse_log_level


ConfigDict = dict[str, object]
ConfigPath = str | Path

_LOGGER = get_logger(__name__)
_ALLOWED_SCHEMA = {
    "project": frozenset({"name"}),
    "reference": frozenset({"fasta", "build"}),
    "samples": frozenset({"sheet"}),
    "runtime": frozenset({"threads", "tmpdir"}),
    "paths": frozenset({"workdir", "outdir"}),
    "logging": frozenset({"level", "file"}),
    "workflow": frozenset({"preset"}),
    "alignment": frozenset(
        {
            "tool",
            "output_format",
            "threads",
            "memory_mb",
            "runtime_minutes",
            "overwrite",
            "pbmm2_preset",
            "pbmm2_log_level",
            "index_threads",
            "bam_index_format",
        }
    ),
    "small": frozenset(
        {
            "enabled",
            "execution_mode",
            "deepvariant_executable",
            "deepvariant_image",
            "model_type",
            "threads",
            "memory_mb",
            "runtime_minutes",
            "overwrite",
        }
    ),
    "sv": frozenset(
        {
            "enabled",
            "overwrite",
            "sawfish",
            "sniffles2",
            "pbsv",
            "cutesv",
            "finalization",
            "harmonization",
        }
    ),
    "assembly_sv": frozenset({"enabled", "overwrite", "pav", "svim_asm"}),
    "assembly": frozenset(
        {
            "enabled",
            "backend",
            "executable",
            "threads",
            "memory_mb",
            "runtime_minutes",
            "overwrite",
        }
    ),
    "phasing": frozenset(
        {
            "enabled",
            "backend",
            "executable",
            "tabix_executable",
            "threads",
            "memory_mb",
            "runtime_minutes",
            "overwrite",
        }
    ),
    "tr": frozenset(
        {
            "enabled",
            "catalog",
            "catalog_reference_build",
            "executable",
            "bcftools_executable",
            "samtools_executable",
            "threads",
            "memory_mb",
            "runtime_minutes",
            "preset",
            "karyotype",
            "overwrite",
        }
    ),
    "review": frozenset(
        {
            "enabled",
            "selection_file",
            "igv_executable",
            "flank_bp",
            "threads",
            "memory_mb",
            "runtime_minutes",
            "overwrite",
        }
    ),
    "annotation": frozenset(
        {
            "enabled", "input_manifest", "overwrite", "threads",
            "memory_mb", "runtime_minutes", "annovar_enabled",
            "annovar_executable", "annovar_version", "annovar_database_root",
            "annovar_database_version", "annovar_protocols", "annovar_operations",
            "vep_enabled", "vep_executable", "vep_cache_directory",
            "vep_cache_version", "vep_species", "vep_assembly",
            "overlap_enabled", "gene_bed", "gene_version", "exon_bed",
            "exon_version", "regulatory_bed", "regulatory_version",
            "repeat_bed", "repeat_version", "segdup_bed", "segdup_version",
            "functional_enabled", "functional_selection_file",
            "alphagenome_model_version", "alphagenome_modalities",
        }
    ),
    "cohort": frozenset(
        {
            "enabled", "cohort_id", "input_manifest", "overwrite",
            "small_variants", "sv", "tr",
        }
    ),
    "benchmark": frozenset(
        {"enabled", "benchmark_id", "sample_id", "overwrite", "small_variants", "sv", "assembly_sv", "tr"}
    ),
}
_REQUIRED_FIELDS = {
    "project": frozenset({"name"}),
    "reference": frozenset({"fasta", "build"}),
    "samples": frozenset({"sheet"}),
    "runtime": frozenset({"threads", "tmpdir"}),
    "paths": frozenset({"workdir", "outdir"}),
    "logging": frozenset({"level", "file"}),
    "workflow": frozenset({"preset"}),
    "alignment": frozenset(
        {
            "tool",
            "output_format",
            "threads",
            "memory_mb",
            "runtime_minutes",
            "overwrite",
            "pbmm2_preset",
            "pbmm2_log_level",
            "index_threads",
            "bam_index_format",
        }
    ),
}
_PATH_FIELDS = (
    ("reference", "fasta"),
    ("samples", "sheet"),
    ("runtime", "tmpdir"),
    ("paths", "workdir"),
    ("paths", "outdir"),
    ("logging", "file"),
    ("tr", "catalog"),
    ("review", "selection_file"),
    ("annotation", "input_manifest"),
    ("annotation", "annovar_database_root"),
    ("annotation", "vep_cache_directory"),
    ("annotation", "gene_bed"),
    ("annotation", "exon_bed"),
    ("annotation", "regulatory_bed"),
    ("annotation", "repeat_bed"),
    ("annotation", "segdup_bed"),
    ("annotation", "functional_selection_file"),
    ("cohort", "input_manifest"),
)
_WORKFLOW_PRESETS = frozenset(
    {"fast", "standard", "comprehensive", "cohort", "trio"}
)
_ALIGNMENT_TOOLS = frozenset({"pbmm2", "minimap2"})
_ALIGNMENT_FORMATS = frozenset({"bam", "cram"})
_PBMM2_PRESETS = frozenset({"CCS", "HIFI"})
_PBMM2_LOG_LEVELS = frozenset({"TRACE", "DEBUG", "INFO", "WARN", "FATAL"})
_BAM_INDEX_FORMATS = frozenset({"auto", "bai", "csi"})
_DEEPVARIANT_EXECUTION_MODES = frozenset({"native", "docker", "apptainer"})
_DEEPVARIANT_MODEL_TYPES = frozenset({"PACBIO"})
_TRGT_PRESETS = frozenset({"wgs", "targeted"})
_TRGT_KARYOTYPES = frozenset({"auto", "XX", "XY"})
_SV_CALLER_KEYS = {
    "sawfish": frozenset(
        {"enabled", "executable", "threads", "memory_mb", "runtime_minutes", "disable_cnv"}
    ),
    "sniffles2": frozenset(
        {"enabled", "executable", "threads", "memory_mb", "runtime_minutes", "minimum_support", "minimum_sv_length"}
    ),
    "pbsv": frozenset(
        {"enabled", "executable", "threads", "memory_mb", "runtime_minutes"}
    ),
    "cutesv": frozenset(
        {
            "enabled", "executable", "threads", "memory_mb", "runtime_minutes",
            "minimum_support", "minimum_sv_size", "max_cluster_bias_ins",
            "diff_ratio_merging_ins", "max_cluster_bias_del", "diff_ratio_merging_del",
            "genotype",
        }
    ),
    "finalization": frozenset({"bgzip_executable", "tabix_executable"}),
    "harmonization": frozenset(
        {
            "enabled", "backend", "jasmine_executable", "truvari_executable",
            "bgzip_executable", "tabix_executable", "threads", "memory_mb",
            "runtime_minutes", "max_dist", "distance_type", "overwrite", "input_vcfs",
        }
    ),
}

_ASSEMBLY_SV_KEYS = {
    "pav": frozenset(
        {"enabled", "executable", "snakefile", "version", "threads", "memory_mb", "runtime_minutes"}
    ),
    "svim_asm": frozenset(
        {
            "enabled", "executable", "minimap2_executable", "samtools_executable",
            "bgzip_executable", "tabix_executable", "threads", "memory_mb", "runtime_minutes",
        }
    ),
}

@dataclass
class HiFiVarConfig(Mapping[str, object]):
    """Merged HiFiVar configuration with lightweight source provenance.

    The object exposes mapping-style read access. Nested values remain mutable
    for callers that need an in-memory adjustment, while every load owns an
    independent deep copy and cannot contaminate a later load.
    """

    _data: ConfigDict = field(repr=False)
    sources: dict[str, Path]

    def __post_init__(self) -> None:
        """Detach stored values and provenance from caller-owned mappings."""
        self._data = deepcopy(self._data)
        self.sources = dict(self.sources)

    def __getitem__(self, key: str) -> object:
        """Return a top-level configuration value."""
        return self._data[key]

    def __iter__(self) -> Iterator[str]:
        """Iterate over top-level configuration keys."""
        return iter(self._data)

    def __len__(self) -> int:
        """Return the number of top-level configuration sections."""
        return len(self._data)

    def __repr__(self) -> str:
        """Describe structure and provenance without dumping config values."""
        return (
            f"HiFiVarConfig(sections={tuple(self._data)!r}, "
            f"sources={tuple(self.sources)!r})"
        )

    def to_dict(self) -> ConfigDict:
        """Return an independent dictionary containing the effective config."""
        return deepcopy(self._data)


def load_yaml(path: ConfigPath) -> ConfigDict:
    """Load one UTF-8 YAML configuration file with ``yaml.safe_load``.

    An empty YAML document is represented by an empty dictionary. Any non-empty
    document must have a mapping at its root.
    """
    config_path = Path(path).expanduser()
    _LOGGER.debug("Loading YAML configuration from %s", config_path)

    try:
        text = config_path.read_text(encoding="utf-8")
    except OSError as error:
        raise ConfigurationError(
            f"Unable to read configuration file '{config_path}': {error}"
        ) from error

    try:
        loaded = yaml.safe_load(text)
    except yaml.YAMLError as error:
        raise ConfigurationError(
            f"Invalid YAML in configuration file '{config_path}': {error}"
        ) from error

    if loaded is None:
        return {}
    if not isinstance(loaded, dict):
        raise ConfigurationError(
            f"Configuration file '{config_path}' must contain a mapping at "
            f"the YAML root, not {type(loaded).__name__}."
        )

    return loaded


def deep_merge(
    base: Mapping[str, object],
    override: Mapping[str, object],
) -> ConfigDict:
    """Recursively merge mappings without mutating either input.

    Nested mappings are merged. Every other override value, including a list,
    replaces the earlier value as an independent deep copy.
    """
    merged: ConfigDict = deepcopy(dict(base))
    for key, override_value in override.items():
        base_value = merged.get(key)
        if isinstance(base_value, Mapping) and isinstance(
            override_value, Mapping
        ):
            merged[key] = deep_merge(base_value, override_value)
        else:
            merged[key] = deepcopy(override_value)
    return merged


def merge_configs(*configs: Mapping[str, object]) -> ConfigDict:
    """Merge configuration layers from lowest to highest priority."""
    merged: ConfigDict = {}
    for config in configs:
        merged = deep_merge(merged, config)
    return merged


def validate_config(
    config: Mapping[str, object],
    *,
    require_complete: bool = False,
) -> None:
    """Validate current HiFiVar schema keys and basic value types.

    Partial mappings such as presets are accepted by default. ``load_config``
    requests complete validation after all layers have been merged.
    """
    if not isinstance(config, Mapping):
        raise ConfigurationError("Configuration root must be a mapping.")

    _validate_known_keys(config)
    if require_complete:
        _validate_required_fields(config)

    project = _optional_section(config, "project")
    if project is not None and "name" in project:
        name = project["name"]
        if not isinstance(name, str) or not name.strip():
            raise ConfigurationError("project.name must be a non-empty string.")

    reference = _optional_section(config, "reference")
    if reference is not None:
        if "fasta" in reference:
            _validate_optional_path(reference["fasta"], "reference.fasta")
        if "build" in reference:
            build = reference["build"]
            if build is not None and (
                not isinstance(build, str) or not build.strip()
            ):
                raise ConfigurationError(
                    "reference.build must be a non-empty string or null."
                )

    samples = _optional_section(config, "samples")
    if samples is not None and "sheet" in samples:
        _validate_optional_path(samples["sheet"], "samples.sheet")

    runtime = _optional_section(config, "runtime")
    if runtime is not None:
        if "threads" in runtime:
            threads = runtime["threads"]
            if (
                not isinstance(threads, int)
                or isinstance(threads, bool)
                or threads <= 0
            ):
                raise ConfigurationError(
                    "runtime.threads must be a positive integer."
                )
        if "tmpdir" in runtime:
            _validate_optional_path(runtime["tmpdir"], "runtime.tmpdir")

    paths = _optional_section(config, "paths")
    if paths is not None:
        for key in ("workdir", "outdir"):
            if key in paths:
                _validate_optional_path(paths[key], f"paths.{key}")

    logging_config = _optional_section(config, "logging")
    if logging_config is not None:
        if "level" in logging_config:
            try:
                parse_log_level(logging_config["level"])
            except ConfigurationError as error:
                raise ConfigurationError(
                    f"Invalid logging.level: {error}"
                ) from error
        if "file" in logging_config:
            _validate_optional_path(logging_config["file"], "logging.file")

    workflow = _optional_section(config, "workflow")
    if workflow is not None and "preset" in workflow:
        preset = workflow["preset"]
        if not isinstance(preset, str) or preset not in _WORKFLOW_PRESETS:
            allowed = ", ".join(sorted(_WORKFLOW_PRESETS))
            raise ConfigurationError(
                f"workflow.preset must be one of: {allowed}."
            )

    alignment = _optional_section(config, "alignment")
    if alignment is not None:
        if "tool" in alignment:
            _validate_choice(
                alignment["tool"],
                "alignment.tool",
                _ALIGNMENT_TOOLS,
            )
        if "output_format" in alignment:
            _validate_choice(
                alignment["output_format"],
                "alignment.output_format",
                _ALIGNMENT_FORMATS,
            )
        for key in ("threads", "memory_mb", "runtime_minutes", "index_threads"):
            if key in alignment:
                _validate_positive_integer(
                    alignment[key],
                    f"alignment.{key}",
                )
        if "overwrite" in alignment and not isinstance(
            alignment["overwrite"],
            bool,
        ):
            raise ConfigurationError("alignment.overwrite must be a boolean.")
        if "pbmm2_preset" in alignment:
            _validate_choice(
                alignment["pbmm2_preset"],
                "alignment.pbmm2_preset",
                _PBMM2_PRESETS,
                case_sensitive=False,
            )
        if "pbmm2_log_level" in alignment:
            _validate_choice(
                alignment["pbmm2_log_level"],
                "alignment.pbmm2_log_level",
                _PBMM2_LOG_LEVELS,
                case_sensitive=False,
            )
        if "bam_index_format" in alignment:
            _validate_choice(
                alignment["bam_index_format"],
                "alignment.bam_index_format",
                _BAM_INDEX_FORMATS,
                case_sensitive=False,
            )

    small = _optional_section(config, "small")
    if small is not None:
        if "enabled" in small and not isinstance(small["enabled"], bool):
            raise ConfigurationError("small.enabled must be a boolean.")
        if "execution_mode" in small:
            _validate_choice(
                small["execution_mode"],
                "small.execution_mode",
                _DEEPVARIANT_EXECUTION_MODES,
                case_sensitive=False,
            )
        if "deepvariant_executable" in small:
            executable = small["deepvariant_executable"]
            if not isinstance(executable, str) or not executable.strip():
                raise ConfigurationError(
                    "small.deepvariant_executable must be a non-empty string."
                )
        if "deepvariant_image" in small:
            image = small["deepvariant_image"]
            if image is not None and (
                not isinstance(image, str) or not image.strip()
            ):
                raise ConfigurationError(
                    "small.deepvariant_image must be a non-empty string or null."
                )
        if "model_type" in small:
            _validate_choice(
                small["model_type"],
                "small.model_type",
                _DEEPVARIANT_MODEL_TYPES,
                case_sensitive=False,
            )
        for key in ("threads", "memory_mb", "runtime_minutes"):
            if key in small:
                _validate_positive_integer(small[key], f"small.{key}")
        if "overwrite" in small and not isinstance(small["overwrite"], bool):
            raise ConfigurationError("small.overwrite must be a boolean.")
        mode = small.get("execution_mode")
        image = small.get("deepvariant_image")
        if isinstance(mode, str):
            normalized_mode = mode.casefold()
            if normalized_mode in {"docker", "apptainer"} and (
                not isinstance(image, str) or not image.strip()
            ):
                raise ConfigurationError(
                    f"small.deepvariant_image is required for {normalized_mode} mode."
                )
            if normalized_mode == "native" and image is not None:
                raise ConfigurationError(
                    "small.deepvariant_image must be null for native mode."
                )

    sv = _optional_section(config, "sv")
    if sv is not None:
        for key in ("enabled", "overwrite"):
            if key in sv and not isinstance(sv[key], bool):
                raise ConfigurationError(f"sv.{key} must be a boolean.")
        for section_name, allowed_keys in _SV_CALLER_KEYS.items():
            subsection = sv.get(section_name)
            if subsection is None:
                continue
            if not isinstance(subsection, Mapping):
                raise ConfigurationError(f"sv.{section_name} must be a mapping.")
            for key in subsection:
                if not isinstance(key, str) or key not in allowed_keys:
                    raise ConfigurationError(f"Unknown configuration key: sv.{section_name}.{key}.")
            if section_name == "harmonization":
                for key in ("enabled", "overwrite"):
                    if key in subsection and not isinstance(subsection[key], bool):
                        raise ConfigurationError(f"sv.harmonization.{key} must be a boolean.")
                if "backend" in subsection:
                    _validate_choice(subsection["backend"], "sv.harmonization.backend", frozenset({"jasmine"}))
                for key in ("jasmine_executable", "truvari_executable", "bgzip_executable", "tabix_executable", "distance_type"):
                    if key in subsection and (not isinstance(subsection[key], str) or not subsection[key].strip()):
                        raise ConfigurationError(f"sv.harmonization.{key} must be a non-empty string.")
                for key in ("threads", "memory_mb", "runtime_minutes", "max_dist"):
                    if key in subsection:
                        _validate_positive_integer(subsection[key], f"sv.harmonization.{key}")
                if "input_vcfs" in subsection:
                    inputs = subsection["input_vcfs"]
                    if not isinstance(inputs, Mapping):
                        raise ConfigurationError("sv.harmonization.input_vcfs must be a mapping.")
                    allowed_callers = {"sawfish", "sniffles2", "pbsv", "cutesv", "pav", "svim_asm"}
                    for sample_id, callers in inputs.items():
                        if not isinstance(sample_id, str) or not sample_id.strip():
                            raise ConfigurationError("sv.harmonization.input_vcfs sample IDs must be non-empty strings.")
                        if not isinstance(callers, Mapping):
                            raise ConfigurationError(f"sv.harmonization.input_vcfs.{sample_id} must be a mapping.")
                        for caller, path in callers.items():
                            if caller not in allowed_callers:
                                raise ConfigurationError(
                                    f"Unknown Phase 9 external caller for sample {sample_id}: {caller}."
                                )
                            if not isinstance(path, str) or not path.strip():
                                raise ConfigurationError(
                                    f"sv.harmonization.input_vcfs.{sample_id}.{caller} must be a non-empty path."
                                )
                continue
            for key, value in subsection.items():
                field = f"sv.{section_name}.{key}"
                if key in {"enabled", "disable_cnv", "genotype"}:
                    if not isinstance(value, bool):
                        raise ConfigurationError(f"{field} must be a boolean.")
                elif key.endswith("executable") or key == "executable":
                    if not isinstance(value, str) or not value.strip():
                        raise ConfigurationError(f"{field} must be a non-empty string.")
                elif key in {"diff_ratio_merging_ins", "diff_ratio_merging_del"}:
                    if not isinstance(value, (int, float)) or isinstance(value, bool) or not 0 <= float(value) <= 1:
                        raise ConfigurationError(f"{field} must be between 0 and 1.")
                elif key == "minimum_support" and section_name == "sniffles2" and value is None:
                    continue
                else:
                    _validate_positive_integer(value, field)

    assembly_sv = _optional_section(config, "assembly_sv")
    if assembly_sv is not None:
        for key in ("enabled", "overwrite"):
            if key in assembly_sv and not isinstance(assembly_sv[key], bool):
                raise ConfigurationError(f"assembly_sv.{key} must be a boolean.")
        for section_name, allowed_keys in _ASSEMBLY_SV_KEYS.items():
            subsection = assembly_sv.get(section_name)
            if subsection is None:
                continue
            if not isinstance(subsection, Mapping):
                raise ConfigurationError(f"assembly_sv.{section_name} must be a mapping.")
            for key in subsection:
                if key not in allowed_keys:
                    raise ConfigurationError(f"Unknown configuration key: assembly_sv.{section_name}.{key}.")
            for key, value in subsection.items():
                field = f"assembly_sv.{section_name}.{key}"
                if key == "enabled":
                    if not isinstance(value, bool):
                        raise ConfigurationError(f"{field} must be a boolean.")
                elif key in {"threads", "memory_mb", "runtime_minutes"}:
                    _validate_positive_integer(value, field)
                elif not isinstance(value, str) or not value.strip():
                    raise ConfigurationError(f"{field} must be a non-empty string.")

    tr = _optional_section(config, "tr")
    if tr is not None:
        for key in ("enabled", "overwrite"):
            if key in tr and not isinstance(tr[key], bool):
                raise ConfigurationError(f"tr.{key} must be a boolean.")
        if "catalog" in tr:
            _validate_optional_path(tr["catalog"], "tr.catalog")
        if "catalog_reference_build" in tr:
            build = tr["catalog_reference_build"]
            if build is not None and (not isinstance(build, str) or not build.strip()):
                raise ConfigurationError("tr.catalog_reference_build must be a non-empty string or null.")
        for key in ("executable", "bcftools_executable", "samtools_executable"):
            if key in tr and (not isinstance(tr[key], str) or not str(tr[key]).strip()):
                raise ConfigurationError(f"tr.{key} must be a non-empty string.")
        for key in ("threads", "memory_mb", "runtime_minutes"):
            if key in tr:
                _validate_positive_integer(tr[key], f"tr.{key}")
        if "preset" in tr:
            _validate_choice(tr["preset"], "tr.preset", _TRGT_PRESETS, case_sensitive=False)
        if "karyotype" in tr:
            _validate_choice(tr["karyotype"], "tr.karyotype", _TRGT_KARYOTYPES)

    phasing = _optional_section(config, "phasing")
    if phasing is not None:
        for key in ("enabled", "overwrite"):
            if key in phasing and not isinstance(phasing[key], bool):
                raise ConfigurationError(f"phasing.{key} must be a boolean.")
        if "backend" in phasing:
            _validate_choice(phasing["backend"], "phasing.backend", frozenset({"hiphase"}))
        for key in ("executable", "tabix_executable"):
            if key in phasing and (
                not isinstance(phasing[key], str) or not phasing[key].strip()
            ):
                raise ConfigurationError(f"phasing.{key} must be a non-empty string.")
        for key in ("threads", "memory_mb", "runtime_minutes"):
            if key in phasing:
                _validate_positive_integer(phasing[key], f"phasing.{key}")
    review = _optional_section(config, "review")
    if review is not None:
        for key in ("enabled", "overwrite"):
            if key in review and not isinstance(review[key], bool):
                raise ConfigurationError(f"review.{key} must be a boolean.")
        if "selection_file" in review:
            _validate_optional_path(review["selection_file"], "review.selection_file")
        if review.get("enabled") and not review.get("selection_file"):
            raise ConfigurationError("review.enabled requires review.selection_file.")
        if "igv_executable" in review and (
            not isinstance(review["igv_executable"], str)
            or not review["igv_executable"].strip()
        ):
            raise ConfigurationError("review.igv_executable must be a non-empty string.")
        if "flank_bp" in review:
            flank = review["flank_bp"]
            if not isinstance(flank, int) or isinstance(flank, bool) or flank < 0:
                raise ConfigurationError("review.flank_bp must be a non-negative integer.")
        for key in ("threads", "memory_mb", "runtime_minutes"):
            if key in review:
                _validate_positive_integer(review[key], f"review.{key}")
    annotation = _optional_section(config, "annotation")
    if annotation is not None:
        boolean_keys = (
            "enabled", "overwrite", "annovar_enabled", "vep_enabled",
            "overlap_enabled", "functional_enabled",
        )
        for key in boolean_keys:
            if key in annotation and not isinstance(annotation[key], bool):
                raise ConfigurationError(f"annotation.{key} must be a boolean.")
        path_keys = (
            "input_manifest", "annovar_database_root", "vep_cache_directory",
            "gene_bed", "exon_bed", "regulatory_bed", "repeat_bed",
            "segdup_bed", "functional_selection_file",
        )
        for key in path_keys:
            if key in annotation:
                _validate_optional_path(annotation[key], f"annotation.{key}")
        string_keys = (
            "annovar_executable", "annovar_version", "annovar_database_version",
            "vep_executable", "vep_cache_version", "vep_species", "vep_assembly",
            "gene_version", "exon_version", "regulatory_version",
            "repeat_version", "segdup_version", "alphagenome_model_version",
        )
        for key in string_keys:
            value = annotation.get(key)
            if value is not None and (not isinstance(value, str) or not value.strip()):
                raise ConfigurationError(f"annotation.{key} must be a non-empty string or null.")
        for key in ("threads", "memory_mb", "runtime_minutes"):
            if key in annotation:
                _validate_positive_integer(annotation[key], f"annotation.{key}")
        for key in ("annovar_protocols", "annovar_operations", "alphagenome_modalities"):
            if key in annotation:
                value = annotation[key]
                if not isinstance(value, list) or any(
                    not isinstance(item, str) or not item.strip() for item in value
                ):
                    raise ConfigurationError(f"annotation.{key} must be a list of non-empty strings.")
        if "annovar_protocols" in annotation and "annovar_operations" in annotation and (
            len(annotation["annovar_protocols"]) != len(annotation["annovar_operations"])
        ):
            raise ConfigurationError("ANNOVAR protocol/operation lists must have equal length.")
        if annotation.get("enabled"):
            if not annotation.get("input_manifest"):
                raise ConfigurationError("annotation.enabled requires annotation.input_manifest.")
            if not any(annotation.get(key) is True for key in ("annovar_enabled", "vep_enabled", "overlap_enabled")):
                raise ConfigurationError("annotation.enabled requires at least one annotation backend.")
        if annotation.get("annovar_enabled") and require_complete:
            required = (
                "annovar_executable", "annovar_version", "annovar_database_root",
                "annovar_database_version", "annovar_protocols", "annovar_operations",
            )
            missing = [key for key in required if not annotation.get(key)]
            if missing:
                raise ConfigurationError(f"ANNOVAR configuration is incomplete: {missing!r}.")
        if annotation.get("vep_enabled") and require_complete:
            required = (
                "vep_executable", "vep_cache_directory", "vep_cache_version",
                "vep_species", "vep_assembly",
            )
            missing = [key for key in required if not annotation.get(key)]
            if missing:
                raise ConfigurationError(f"VEP configuration is incomplete: {missing!r}.")
        if annotation.get("overlap_enabled"):
            pairs = (
                ("gene_bed", "gene_version"), ("exon_bed", "exon_version"),
                ("regulatory_bed", "regulatory_version"),
                ("repeat_bed", "repeat_version"), ("segdup_bed", "segdup_version"),
            )
            if not any(annotation.get(path_key) for path_key, _ in pairs):
                raise ConfigurationError("Region overlap requires at least one BED database.")
            for path_key, version_key in pairs:
                if bool(annotation.get(path_key)) != bool(annotation.get(version_key)):
                    raise ConfigurationError(
                        f"annotation.{path_key} and annotation.{version_key} must be configured together."
                    )
        if annotation.get("functional_enabled"):
            required = ("functional_selection_file", "alphagenome_model_version", "alphagenome_modalities")
            missing = [key for key in required if not annotation.get(key)]
            if missing:
                raise ConfigurationError(f"Functional prioritization configuration is incomplete: {missing!r}.")
    cohort = _optional_section(config, "cohort")
    if cohort is not None:
        for key in ("enabled", "overwrite"):
            if key in cohort and not isinstance(cohort[key], bool):
                raise ConfigurationError(f"cohort.{key} must be a boolean.")
        if "cohort_id" in cohort:
            value = cohort["cohort_id"]
            if value is not None and (not isinstance(value, str) or not value.strip()):
                raise ConfigurationError("cohort.cohort_id must be a non-empty string or null.")
        if "input_manifest" in cohort:
            _validate_optional_path(cohort["input_manifest"], "cohort.input_manifest")
        if cohort.get("enabled") and (not cohort.get("cohort_id") or not cohort.get("input_manifest")):
            raise ConfigurationError("cohort.enabled requires cohort.cohort_id and cohort.input_manifest.")
        track_keys = {
            "small_variants": frozenset({"enabled", "glnexus_executable", "bcftools_executable", "preset", "threads", "memory_gb", "runtime_minutes"}),
            "sv": frozenset({"enabled", "threads", "memory_mb", "runtime_minutes"}),
            "tr": frozenset({"enabled", "threads", "memory_mb", "runtime_minutes"}),
        }
        for track_name, allowed in track_keys.items():
            track = cohort.get(track_name)
            if track is None:
                continue
            if not isinstance(track, Mapping):
                raise ConfigurationError(f"cohort.{track_name} must be a mapping.")
            unknown = set(track).difference(allowed)
            if unknown:
                raise ConfigurationError(f"Unknown configuration key: cohort.{track_name}.{sorted(unknown)[0]}.")
            if "enabled" in track and not isinstance(track["enabled"], bool):
                raise ConfigurationError(f"cohort.{track_name}.enabled must be a boolean.")
            for key in ("threads", "memory_mb", "memory_gb", "runtime_minutes"):
                if key in track:
                    _validate_positive_integer(track[key], f"cohort.{track_name}.{key}")
            for key in ("glnexus_executable", "bcftools_executable", "preset"):
                if key in track and (not isinstance(track[key], str) or not track[key].strip()):
                    raise ConfigurationError(f"cohort.{track_name}.{key} must be a non-empty string.")
    benchmark = _optional_section(config, "benchmark")
    if benchmark is not None:
        for key in ("enabled", "overwrite"):
            if key in benchmark and not isinstance(benchmark[key], bool):
                raise ConfigurationError(f"benchmark.{key} must be a boolean.")
        for key in ("benchmark_id", "sample_id"):
            value = benchmark.get(key)
            if value is not None and (not isinstance(value, str) or not value.strip()):
                raise ConfigurationError(f"benchmark.{key} must be a non-empty string or null.")
        track_keys = {
            "small_variants": frozenset({"enabled", "query_vcf", "truth_vcf", "truth_version", "truth_source", "confident_bed", "confident_bed_version", "happy_executable", "happy_version", "metrics_compression", "engine", "summary_filter", "threads", "memory_mb", "runtime_minutes", "stratifications"}),
            "sv": frozenset({"enabled", "query_vcf", "truth_vcf", "truth_version", "truth_source", "confident_bed", "confident_bed_version", "region_class", "truvari_executable", "threads", "memory_mb", "runtime_minutes", "refdist", "pctseq", "pctsize", "pctovl", "sizemin", "sizemax", "bnddist", "pass_only", "size_bins"}),
            "assembly_sv": frozenset({"enabled", "query_vcf", "truth_vcf", "truth_version", "truth_source", "confident_bed", "confident_bed_version", "region_class", "truvari_executable", "threads", "memory_mb", "runtime_minutes", "refdist", "pctseq", "pctsize", "pctovl", "sizemin", "sizemax", "bnddist", "pass_only", "size_bins"}),
            "tr": frozenset({"enabled", "query_vcf", "truth_vcf", "truth_version", "truth_source", "catalog_id", "threads", "memory_mb", "runtime_minutes"}),
        }
        any_enabled = False
        for track_name, allowed in track_keys.items():
            track = benchmark.get(track_name)
            if track is None:
                continue
            if not isinstance(track, Mapping):
                raise ConfigurationError(f"benchmark.{track_name} must be a mapping.")
            unknown = set(track).difference(allowed)
            if unknown:
                raise ConfigurationError(f"Unknown configuration key: benchmark.{track_name}.{sorted(unknown)[0]}.")
            if "enabled" in track and not isinstance(track["enabled"], bool):
                raise ConfigurationError(f"benchmark.{track_name}.enabled must be a boolean.")
            any_enabled = any_enabled or track.get("enabled") is True
            for key in ("threads", "memory_mb", "runtime_minutes"):
                if key in track:
                    _validate_positive_integer(track[key], f"benchmark.{track_name}.{key}")
            for key in allowed.intersection({"query_vcf", "truth_vcf", "confident_bed"}):
                if key in track: _validate_optional_path(track[key], f"benchmark.{track_name}.{key}")
            for key in allowed.intersection({"truth_version", "truth_source", "confident_bed_version", "region_class", "happy_executable", "happy_version", "engine", "summary_filter", "truvari_executable", "catalog_id"}):
                value = track.get(key)
                if value is not None and (not isinstance(value, str) or not value.strip()):
                    raise ConfigurationError(f"benchmark.{track_name}.{key} must be a non-empty string or null.")
            if "metrics_compression" in track:
                _validate_choice(
                    track["metrics_compression"],
                    f"benchmark.{track_name}.metrics_compression",
                    frozenset({"gzip", "plain"}),
                )
            for key in ("pctseq", "pctsize", "pctovl"):
                if key in track and track[key] is not None and (isinstance(track[key], bool) or not isinstance(track[key], (int, float)) or not 0 <= float(track[key]) <= 1):
                    raise ConfigurationError(f"benchmark.{track_name}.{key} must be between 0 and 1 or null.")
            for key in ("refdist", "sizemin", "sizemax", "bnddist"):
                if key in track and track[key] is not None and (isinstance(track[key], bool) or not isinstance(track[key], int) or track[key] < 0):
                    raise ConfigurationError(f"benchmark.{track_name}.{key} must be a non-negative integer or null.")
            if "pass_only" in track and not isinstance(track["pass_only"], bool):
                raise ConfigurationError(f"benchmark.{track_name}.pass_only must be a boolean.")
            if "size_bins" in track:
                bins=track["size_bins"]
                if not isinstance(bins,list) or any(isinstance(value,bool) or not isinstance(value,int) or value <= 0 for value in bins) or sorted(set(bins)) != bins:
                    raise ConfigurationError(f"benchmark.{track_name}.size_bins must be unique ascending positive integers.")
            if "stratifications" in track:
                regions = track["stratifications"]
                if not isinstance(regions, list):
                    raise ConfigurationError(f"benchmark.{track_name}.stratifications must be a list.")
                for region in regions:
                    if not isinstance(region, Mapping) or set(region) != {"name", "path", "version", "region_class"} or any(not isinstance(region[key], str) or not region[key].strip() for key in region):
                        raise ConfigurationError("Each benchmark stratification requires non-empty name/path/version/region_class strings.")
            if track.get("enabled"):
                required = {"query_vcf", "truth_vcf", "truth_version", "truth_source"}
                if track_name == "small_variants": required |= {"confident_bed", "confident_bed_version"}
                if track_name == "tr": required |= {"catalog_id"}
                missing = sorted(key for key in required if not track.get(key))
                if missing: raise ConfigurationError(f"benchmark.{track_name} configuration is incomplete: {missing!r}.")
        if benchmark.get("enabled") and (not benchmark.get("benchmark_id") or not benchmark.get("sample_id") or not any_enabled):
            raise ConfigurationError("benchmark.enabled requires benchmark_id, sample_id, and at least one enabled track.")
        if benchmark.get("enabled") and (reference is None or not reference.get("build")):
            raise ConfigurationError("benchmark.enabled requires an explicit reference.build.")
    assembly = _optional_section(config, "assembly")
    if assembly is not None:
        for key in ("enabled", "overwrite"):
            if key in assembly and not isinstance(assembly[key], bool):
                raise ConfigurationError(f"assembly.{key} must be a boolean.")
        if "backend" in assembly:
            _validate_choice(assembly["backend"], "assembly.backend", frozenset({"hifiasm"}))
        if "executable" in assembly and (
            not isinstance(assembly["executable"], str) or not assembly["executable"].strip()
        ):
            raise ConfigurationError("assembly.executable must be a non-empty string.")
        for key in ("threads", "memory_mb", "runtime_minutes"):
            if key in assembly:
                _validate_positive_integer(assembly[key], f"assembly.{key}")



def load_config(
    default_config: ConfigPath,
    preset: ConfigPath | None = None,
    user_config: ConfigPath | None = None,
) -> HiFiVarConfig:
    """Load, validate, and merge default, preset, and user configuration.

    Precedence is deterministic: user config overrides preset, and preset
    overrides defaults. Each layer is checked before merging so an unknown key
    cannot be hidden by a later layer.
    """
    layers: list[ConfigDict] = []
    sources: dict[str, Path] = {}
    requested_sources = (
        ("default", default_config),
        ("preset", preset),
        ("user", user_config),
    )

    for source_name, source_path in requested_sources:
        if source_path is None:
            continue
        normalized_path = Path(source_path).expanduser()
        _LOGGER.debug(
            "Loading %s configuration from %s",
            source_name,
            normalized_path,
        )
        layer = load_yaml(normalized_path)
        validate_config(layer)
        layers.append(layer)
        sources[source_name] = normalized_path

    effective_config = merge_configs(*layers)
    validate_config(effective_config, require_complete=True)
    effective_config = _expand_user_paths(effective_config)

    return HiFiVarConfig(effective_config, sources)


def write_effective_config(
    config: HiFiVarConfig,
    output_path: ConfigPath,
) -> None:
    """Write only the effective configuration as reloadable UTF-8 YAML."""
    effective_config = config.to_dict()
    validate_config(effective_config, require_complete=True)
    destination = Path(output_path).expanduser()
    _LOGGER.debug("Writing effective configuration to %s", destination)

    try:
        destination.parent.mkdir(parents=True, exist_ok=True)
        serialized = yaml.safe_dump(
            effective_config,
            allow_unicode=True,
            default_flow_style=False,
            sort_keys=False,
        )
        destination.write_text(serialized, encoding="utf-8")
    except (OSError, yaml.YAMLError) as error:
        raise ConfigurationError(
            f"Unable to write effective configuration '{destination}': {error}"
        ) from error


def _validate_known_keys(config: Mapping[str, object]) -> None:
    """Reject unknown or non-string keys at current schema levels."""
    for section_name in config:
        if not isinstance(section_name, str):
            raise ConfigurationError(
                "Top-level configuration keys must be strings."
            )
        if section_name not in _ALLOWED_SCHEMA:
            raise ConfigurationError(
                f"Unknown top-level configuration key: {section_name}."
            )

        section = config[section_name]
        if not isinstance(section, Mapping):
            raise ConfigurationError(
                f"Configuration section {section_name} must be a mapping."
            )
        for key in section:
            if not isinstance(key, str):
                raise ConfigurationError(
                    f"Keys in configuration section {section_name} must be "
                    "strings."
                )
            if key not in _ALLOWED_SCHEMA[section_name]:
                raise ConfigurationError(
                    f"Unknown configuration key: {section_name}.{key}."
                )


def _validate_required_fields(config: Mapping[str, object]) -> None:
    """Ensure a merged effective configuration has all foundation fields."""
    for section_name, required_keys in _REQUIRED_FIELDS.items():
        section = config.get(section_name)
        if not isinstance(section, Mapping):
            raise ConfigurationError(
                f"Missing required configuration section: {section_name}."
            )
        missing_keys = required_keys.difference(section)
        if missing_keys:
            missing = ", ".join(
                f"{section_name}.{key}" for key in sorted(missing_keys)
            )
            raise ConfigurationError(
                f"Missing required configuration key(s): {missing}."
            )


def _optional_section(
    config: Mapping[str, object],
    name: str,
) -> Mapping[str, object] | None:
    """Return a schema-validated optional section."""
    section = config.get(name)
    if section is None:
        return None
    if not isinstance(section, Mapping):
        raise ConfigurationError(f"Configuration section {name} must be a mapping.")
    return section


def _validate_optional_path(value: object, field_name: str) -> None:
    """Require path settings to be strings or null without checking existence."""
    if value is not None and not isinstance(value, str):
        raise ConfigurationError(f"{field_name} must be a string path or null.")


def _validate_positive_integer(value: object, field_name: str) -> None:
    """Require a non-boolean positive integer configuration value."""
    if not isinstance(value, int) or isinstance(value, bool) or value <= 0:
        raise ConfigurationError(f"{field_name} must be a positive integer.")


def _validate_choice(
    value: object,
    field_name: str,
    choices: frozenset[str],
    *,
    case_sensitive: bool = True,
) -> None:
    """Require one string from a small explicit vocabulary."""
    if not isinstance(value, str):
        allowed = ", ".join(sorted(choices))
        raise ConfigurationError(f"{field_name} must be one of: {allowed}.")
    selected = value if case_sensitive else value.upper()
    available = choices if case_sensitive else frozenset(
        choice.upper() for choice in choices
    )
    if selected not in available:
        allowed = ", ".join(sorted(choices))
        raise ConfigurationError(f"{field_name} must be one of: {allowed}.")


def _expand_user_paths(config: Mapping[str, object]) -> ConfigDict:
    """Expand leading ``~`` without resolving or checking configured paths."""
    expanded = deepcopy(dict(config))
    for section_name, key in _PATH_FIELDS:
        section = expanded.get(section_name)
        if not isinstance(section, dict):
            continue
        value = section.get(key)
        if isinstance(value, str) and value.startswith("~"):
            try:
                section[key] = str(Path(value).expanduser())
            except RuntimeError as error:
                raise ConfigurationError(
                    f"Unable to expand user path for {section_name}.{key}: "
                    f"{value!r}."
                ) from error
    return expanded


__all__ = [
    "HiFiVarConfig",
    "deep_merge",
    "load_config",
    "load_yaml",
    "merge_configs",
    "validate_config",
    "write_effective_config",
]
