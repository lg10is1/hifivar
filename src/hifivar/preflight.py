"""Opt-in, non-mutating checks for a fresh, bounded calling-launcher run.

This is not a universal workflow validator, benchmark, installer or OS sandbox.
Only explicit tool probes execute existing wrapper version/read-only methods.
"""
from __future__ import annotations

import importlib.util
import math
import os
import platform
import re
import shutil
import sys
from collections import Counter
from collections.abc import Callable
from dataclasses import asdict, dataclass, field
from pathlib import Path

from hifivar import __version__
from hifivar.command import CommandResult, CommandRunner
from hifivar.config import HiFiVarConfig
from hifivar.exceptions import ConfigurationError, HiFiVarError, InputValidationError
from hifivar.launch_inputs import IndexedBamPreflight, parse_header_records, unique_bam_index
from hifivar.launcher import prepare_launch
from hifivar.package_resources import installed_workflow_root
from hifivar.sample import InputType
from hifivar.serialization import utc_now_iso8601


@dataclass(frozen=True)
class PreflightCheck:
    """One scoped observation, distinct from accuracy or complete readiness."""

    key: str
    status: str
    message: str
    required: bool = True
    details: dict[str, object] = field(default_factory=dict)


class ProbeRunner:
    """Keep optional probes shell-free, bounded and without file destinations."""

    def __init__(self, runner: CommandRunner, timeout: float) -> None:
        self.runner, self.timeout, self.steps = runner, timeout, []

    def require_executable(self, executable: str | Path) -> Path:
        return self.runner.require_executable(executable)

    def run(self, command: list[str], **kwargs) -> CommandResult:
        if kwargs.get("stdout_path") or kwargs.get("stderr_path"):
            raise InputValidationError("Preflight probes cannot create tool log files.")
        kwargs["timeout"] = min(kwargs.get("timeout") or self.timeout, self.timeout)
        # Do not copy credentials/Python injection into version/header probes.
        keys = {"PATH", "HOME", "LANG", "LC_ALL", "TZ", "CONDA_PREFIX", "LD_LIBRARY_PATH",
                "DYLD_LIBRARY_PATH", "SYSTEMROOT", "SystemRoot", "WINDIR", "PATHEXT",
                "TMPDIR", "TMP", "TEMP"}
        kwargs.update(env={k: v for k, v in os.environ.items() if k in keys}, inherit_env=False)
        step = {"command": list(map(str, command)), "timeout_seconds": kwargs["timeout"]}
        self.steps.append(step)
        try:
            result = self.runner.run(command, **kwargs)
        except (HiFiVarError, OSError) as error:
            step.update(status="FAILED", error_type=type(error).__name__)
            raise
        step.update(status="COMPLETED", returncode=result.returncode,
                    duration_seconds=result.duration_seconds)
        return result


def _tool_groups(config: dict, runner: ProbeRunner, *, check_primary_alignments: bool = True) -> list[tuple[str, dict, Callable]]:
    """Use existing wrappers as owners of all bioinformatics version commands."""
    from hifivar.cutesv import CuteSvWrapper
    from hifivar.deepvariant import DeepVariantRuntime, DeepVariantWrapper
    from hifivar.hiphase import HiPhaseWrapper
    from hifivar.hifiasm import HifiasmWrapper
    from hifivar.glnexus import GLnexusWrapper
    from hifivar.jasmine import JasmineWrapper
    from hifivar.minimap2 import Minimap2Wrapper
    from hifivar.pbmm2 import Pbmm2Wrapper
    from hifivar.pbsv import PbsvWrapper
    from hifivar.samtools import SamtoolsWrapper
    from hifivar.sawfish import SawfishWrapper
    from hifivar.sniffles2 import Sniffles2Wrapper
    from hifivar.sv import BgzipTabixWrapper
    from hifivar.svim_asm import SvimAsmWrapper
    from hifivar.trgt import TrgtWrapper
    from hifivar.truvari import TruvariWrapper
    from hifivar.vep import VepWrapper

    samtools = config["workflow"]["samtools_executable"]
    groups = ([("samtools", {"samtools": samtools}, SamtoolsWrapper(runner=runner, executable=samtools).detect_version)]
              if check_primary_alignments else [])
    if config["alignment"]["enabled"]:
        a = config["alignment"]
        backend = Pbmm2Wrapper if a["tool"] == "pbmm2" else Minimap2Wrapper
        executable = a[f"{a['tool']}_executable"]
        groups.append((a["tool"], {a["tool"]: executable}, backend(runner=runner, executable=executable).detect_version))
    if config["small"]["enabled"]:
        runtime = DeepVariantRuntime.from_config(config)
        groups.append(("deepvariant", {"deepvariant": runtime.launcher},
                       DeepVariantWrapper(runner=runner, runtime=runtime).detect_version))
    sv = config["sv"]
    if sv["enabled"]:
        for name, wrapper in (("sawfish", SawfishWrapper), ("sniffles2", Sniffles2Wrapper),
                              ("pbsv", PbsvWrapper), ("cutesv", CuteSvWrapper)):
            if sv[name]["enabled"]:
                executable = sv[name]["executable"]
                groups.append((name, {name: executable}, wrapper(runner=runner, executable=executable).detect_version))
        if sv["pbsv"]["enabled"] or sv["cutesv"]["enabled"]:
            f = sv["finalization"]
            groups.append(("sv_finalization", {"bgzip": f["bgzip_executable"], "tabix": f["tabix_executable"]},
                BgzipTabixWrapper(runner=runner, **f).detect_versions))
    h = sv["harmonization"]
    if h["enabled"]:
        groups.append(("jasmine", {"jasmine": h["jasmine_executable"],
            "bgzip": h["bgzip_executable"], "tabix": h["tabix_executable"]},
            JasmineWrapper(runner=runner, executable=h["jasmine_executable"],
                bgzip_executable=h["bgzip_executable"], tabix_executable=h["tabix_executable"]).detect_versions))
        groups.append(("truvari", {"truvari": h["truvari_executable"]},
            TruvariWrapper(runner=runner, executable=h["truvari_executable"]).detect_version))
    if config["tr"]["enabled"]:
        tr = config["tr"]
        groups.append(("trgt", {"trgt": tr["executable"], "bcftools": tr["bcftools_executable"],
            "samtools": tr["samtools_executable"]}, TrgtWrapper(runner=runner, executable=tr["executable"],
                bcftools_executable=tr["bcftools_executable"], samtools_executable=tr["samtools_executable"]).detect_versions))
    if config["phasing"]["enabled"]:
        p = config["phasing"]
        # HiPhase checks tabix availability, not its version, in the existing wrapper.
        groups.append(("hiphase", {"hiphase": p["executable"], "tabix": p["tabix_executable"]},
            HiPhaseWrapper(runner=runner, executable=p["executable"], tabix_executable=p["tabix_executable"]).detect_version))
    if config["assembly"]["enabled"]:
        executable = config["assembly"]["executable"]
        groups.append(("assembly.hifiasm", {"hifiasm": executable},
            HifiasmWrapper(runner=runner, executable=executable).detect_version))
    assembly_sv = config["assembly_sv"]
    if assembly_sv["enabled"] and assembly_sv["svim_asm"]["enabled"]:
        a = assembly_sv["svim_asm"]
        dependencies = {name: a[name + "_executable"]
                        for name in ("minimap2", "samtools", "bgzip", "tabix")}
        dependencies["svim_asm"] = a["executable"]
        groups.append(("assembly_sv.svim_asm", dependencies,
            SvimAsmWrapper(runner=runner, **{key: a[key] for key in (
                "executable", "minimap2_executable", "samtools_executable",
                "bgzip_executable", "tabix_executable")}).detect_versions))
    annotation = config["annotation"]
    if annotation["enabled"] and annotation["vep_enabled"]:
        executable = annotation["vep_executable"]
        groups.append(("annotation.vep", {"vep": executable},
            VepWrapper(runner=runner, executable=executable).detect_version))
    cohort = config["cohort"]
    if cohort["enabled"] and cohort["small_variants"]["enabled"]:
        c = cohort["small_variants"]
        groups.append(("cohort.glnexus", {"glnexus": c["glnexus_executable"],
            "bcftools": c["bcftools_executable"]},
            GLnexusWrapper(runner=runner, executable=c["glnexus_executable"],
                bcftools_executable=c["bcftools_executable"]).detect_versions))
    return groups


def _host_resources(root: Path) -> dict:
    """Report available filesystem/OS metadata without promising quota/headroom."""
    parent = root.parent
    while not parent.exists() and parent != parent.parent:
        parent = parent.parent
    usage = shutil.disk_usage(parent)
    data = {"filesystem_checked_at": str(parent), "disk_free_bytes": usage.free,
            "cpu_count": os.cpu_count(), "parent_writable_advisory": os.access(parent, os.W_OK),
            "quota": "NOT_CHECKED", "required_disk_bytes": "NOT_ESTIMATED",
            "thread_headroom": "NOT_MEASURED", "available_memory": "NOT_MEASURED"}
    if hasattr(os, "statvfs"):
        stat = os.statvfs(parent)
        data["inodes_available"] = stat.f_favail
        data["inodes_total"] = stat.f_files
    try:
        import resource
    except ImportError:
        data["process_limit"] = "NOT_AVAILABLE_ON_PLATFORM"
    else:
        for name, flag in (("process_limit", resource.RLIMIT_NPROC), ("file_descriptor_limit", resource.RLIMIT_NOFILE)):
            soft, hard = resource.getrlimit(flag)
            data[name] = {"soft": soft, "hard": hard, "infinity_value": resource.RLIM_INFINITY}
    return data


def run_preflight(config: HiFiVarConfig, run_directory: str | Path, *, probe_tools: bool = False,
                  timeout_seconds: float = 30, runner: CommandRunner | None = None) -> dict:
    """Check only the fresh hifivar-run scope; do not create or repair anything."""
    if isinstance(timeout_seconds, bool) or not isinstance(timeout_seconds, (int, float)) or not math.isfinite(timeout_seconds) or timeout_seconds <= 0:
        raise ConfigurationError("Preflight timeout must be a finite positive number.")
    checks: list[PreflightCheck] = []
    probes = ProbeRunner(runner or CommandRunner(), timeout_seconds)
    root = Path(run_directory).expanduser().absolute()

    def add(key, status, message, *, required=True, **details):
        checks.append(PreflightCheck(key, status, message, required, details))

    def report():
        counts = dict(Counter(check.status for check in checks))
        state = "BLOCKED" if counts.get("FAIL") else (
            "INCOMPLETE" if any(c.required and c.status == "NOT_CHECKED" for c in checks) else "CHECKED_SCOPE_PASS")
        return {"schema_version": 1, "hifivar_version": __version__, "created_at": utc_now_iso8601(),
                "scope": "FRESH_CALLING_LAUNCHER_ONLY", "status": state,
                "exit_code": {"CHECKED_SCOPE_PASS": 0, "BLOCKED": 1, "INCOMPLETE": 2}[state],
                "run_directory": str(root), "probe_tools": probe_tools, "counts": counts,
                "checks": [asdict(c) for c in checks], "probe_steps": probes.steps,
                "analysis_executed": False, "ready_for_all_workflows": False,
                "scientific_accuracy_verified": False, "raw_content_immutability_verified": False}

    try:
        if root.resolve() != root or root == Path(root.anchor) or root == Path.home() or root.exists():
            raise InputValidationError("Use a new non-root, non-linked run directory; resume is not this preflight's scope.")
        add("destination", "PASS", "Fresh destination planned; no directory created.")
        effective, context, _, inputs = prepare_launch(config, root)
        add("launcher_contract", "PASS", "Enabled tracks, budgets and lightweight inputs satisfy the existing launcher.",
            samples=list(context.sample_ids), input_count=len(inputs))
    except (HiFiVarError, OSError) as error:
        add("launcher_contract", "FAIL", str(error), error_type=type(error).__name__)
        return report()
    add("platform", "PASS" if platform.system() == "Linux" else "FAIL", "Production calling requires Linux.",
        observed=platform.system())
    add("reference_build", "PASS" if context.reference.build else "FAIL", "Reference build must be explicitly declared.",
        declared=context.reference.build, contigs=len(context.reference.contigs),
        validation="FASTA_PREFIX_AND_FAI_NOT_FULL_SEQUENCE_IDENTITY")
    try:
        workflow = installed_workflow_root()
        add("workflow_resources", "PASS", "Installed workflow located.", path=str(workflow))
        resources = _host_resources(root)
        inode_ok = not resources.get("inodes_total") or resources["inodes_available"] > 0
        add("host_resources", "PASS" if resources["disk_free_bytes"] > 0 and resources["parent_writable_advisory"] and inode_ok else "FAIL",
            "Filesystem observations only; no write test or disk requirement estimate.", **resources)
        add("resource_headroom", "NOT_CHECKED", "Quota, actual per-UID thread headroom and WGS disk/RAM sufficiency are not certified.", required=False)
    except (HiFiVarError, OSError) as error:
        add("workflow_or_resources", "FAIL", str(error), error_type=type(error).__name__)
    if effective["tr"]["enabled"]:
        from hifivar.phase5 import resolve_trgt_karyotype
        from hifivar.tr import TandemRepeatCatalog
        try:
            if not effective["tr"]["catalog_reference_build"]:
                raise ConfigurationError("Preflight requires explicit tr.catalog_reference_build.")
            TandemRepeatCatalog(Path(effective["tr"]["catalog"]), effective["tr"]["catalog_reference_build"]).validate(context.reference)
            lengths = {contig.name: contig.length for contig in context.reference.contigs}
            with Path(effective["tr"]["catalog"]).open(encoding="utf-8") as stream:
                for line in stream:
                    line = line.strip()
                    if line and not line.startswith("#"):
                        fields = line.split("\t")
                        if len(fields) < 3 or fields[0] not in lengths:
                            raise InputValidationError("TRGT catalog changed or contains invalid contigs/fields.")
                        if int(fields[2]) > lengths[fields[0]]:
                            raise InputValidationError("TRGT catalog interval exceeds its declared reference contig.")
            for record in context.samples:
                resolve_trgt_karyotype(effective["tr"]["karyotype"], record.sex, record.sample.sample_id)
            add("tr_catalog", "PASS", "Catalog structure/exact names/bounds/declared build and karyotype checked; not pathogenicity.")
        except (HiFiVarError, OSError, ValueError) as error:
            add("tr_catalog", "FAIL", str(error), error_type=type(error).__name__)
    spec = importlib.util.find_spec("snakemake")
    add("snakemake_package", "PASS" if spec else "FAIL", "Snakemake must be installed in the active Python environment.")
    if probe_tools and spec:
        try:
            result = probes.run([sys.executable, "-I", "-B", "-m", "snakemake", "--version"])
            version = (result.stdout or "").strip()
            if not re.fullmatch(r"[89]\.\d+(?:\.\d+)?", version):
                raise ConfigurationError("Installed isolated Snakemake version must satisfy >=8,<10.")
            add("snakemake_version", "PASS", "Version resolved in isolated active interpreter.", version=version)
        except (HiFiVarError, OSError) as error:
            add("snakemake_version", "FAIL", "Isolated Snakemake version probe failed.", error_type=type(error).__name__)
    elif spec:
        add("snakemake_version", "NOT_CHECKED", "Opt in to tool probes to resolve the isolated interpreter version.")
    primary_alignment_checks = effective["alignment"]["enabled"] or any(
        item.sample.input.input_type is InputType.BAM for item in context.samples)
    for name, executables, callback in _tool_groups(effective, probes, check_primary_alignments=primary_alignment_checks):
        found = {label: probes.runner.find_executable(executable) for label, executable in executables.items()}
        missing = [label for label, path in found.items() if path is None]
        add(f"tool.{name}.presence", "FAIL" if missing else "PASS", "Configured enabled dependency lookup.",
            configured=executables, located={k: str(v) if v else None for k, v in found.items()}, missing=missing)
        if missing:
            continue
        container = name == "deepvariant" and effective["small"]["execution_mode"] != "native"
        if container or not probe_tools:
            add(f"tool.{name}.version", "NOT_CHECKED",
                "Container is not started or pulled by preflight." if container else "Version probe not requested.")
        else:
            try:
                version = callback()
                add(f"tool.{name}.version", "PASS", "Existing wrapper resolved the version, not tool accuracy.",
                    version=version, version_source="command",
                    validation_scope="WRAPPER_NATIVE_VERSION_ONLY_NOT_ANALYSIS_OR_DATABASE_CONTENT")
                if name == "jasmine":
                    from hifivar.jasmine import VCF44_MIN_TABIX
                    observed = version["tabix"]
                    compatible = tuple(int(item) for item in observed.split(".")) >= VCF44_MIN_TABIX
                    add("tool.jasmine.vcf44_indexer", "PASS" if compatible else "FAIL",
                        "Launcher harmonization deployment requires a VCFv4.4-compatible indexer; no header downgrade.",
                        observed=observed, minimum=".".join(map(str, VCF44_MIN_TABIX)))
                    if effective['sv']['harmonization'].get('input_identity', {}).get('enabled'):
                        from hifivar.sv_input_identity import SOURCE_IDENTITY_MIN_HTSLIB
                        observed_bgzip = version['bgzip']
                        bgzip_compatible = tuple(int(item) for item in observed_bgzip.split('.')) >= SOURCE_IDENTITY_MIN_HTSLIB
                        add('tool.jasmine.source_identity_compressor', 'PASS' if bgzip_compatible else 'FAIL',
                            'Explicit source identity uses the same bgzip compatibility floor as its wrapper.',
                            observed=observed_bgzip, minimum='.'.join(map(str, SOURCE_IDENTITY_MIN_HTSLIB)))
            except (HiFiVarError, OSError) as error:
                add(f"tool.{name}.version", "FAIL", "Wrapper version probe failed; inspect the enabled environment.", error_type=type(error).__name__)
    # Input-contract routing is broader than this bounded native probe registry.
    # Do not let a newly routable external backend inherit a false complete pass.
    unprobed = []
    # PAV's wrapper exposes an adapter version, not an observed PAV version;
    # ANNOVAR uses an explicit declared source version. Neither is upgraded to
    # native evidence by checking an adapter or echoing configuration.
    if effective['assembly_sv']['enabled'] and effective['assembly_sv']['pav']['enabled']:
        unprobed.append('assembly_sv.pav')
    if effective['annotation']['enabled'] and effective['annotation']['annovar_enabled']:
        unprobed.append('annotation.annovar')
    for name in unprobed:
        add('tool.' + name + '.probe_scope', 'NOT_CHECKED',
            'This enabled external backend is not covered by the bounded probe registry; actual wrapper acceptance is separate.')
    for record in context.samples:
        sample = record.sample
        if sample.input.input_type is InputType.FASTQ:
            add(f"input.{sample.sample_id}", "PASS", "FASTQ prefix validated; full molecule/record scan and future alignment are not performed.")
            continue
        if not probe_tools:
            add(f"input.{sample.sample_id}", "NOT_CHECKED", "BAM header/sample/reference/EOF/indexed access require explicit tool probes.")
            continue
        try:
            bam = sample.input.files[0].resolve()
            checker = IndexedBamPreflight(probes, effective["workflow"]["samtools_executable"])
            header = checker.read_header(bam)
            if parse_header_records(header)["@SQ"]:
                outcome = checker.run(bam, unique_bam_index(bam), sample.sample_id, context.reference)
                add(f"input.{sample.sample_id}", "PASS", "Existing wrapper verified header/EOF and one indexed region, not all records.", **outcome)
            else:
                from hifivar.ubam import validate_ubam_header
                if not effective["alignment"]["enabled"] or effective["alignment"]["tool"] != "pbmm2":
                    raise ConfigurationError("uBAM input requires explicitly enabled native pbmm2 routing.")
                validate_ubam_header(header, sample.sample_id)
                add(f"input.{sample.sample_id}", "NOT_CHECKED", "uBAM header compatible; complete flag/RG record scan deliberately not run.", kind="UBAM")
        except (HiFiVarError, OSError) as error:
            add(f"input.{sample.sample_id}", "FAIL", "BAM input probe failed; no input repair or sample/contig conversion occurred.", error_type=type(error).__name__)
    return report()


__all__ = ["PreflightCheck", "run_preflight"]
