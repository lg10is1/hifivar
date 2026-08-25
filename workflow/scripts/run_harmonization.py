"""Snakemake bridge for Phase 9 Jasmine harmonization and Truvari comparison."""

from pathlib import Path

from hifivar.assembly_sv import SVEvidenceSource
from hifivar.context import AnalysisContext
from hifivar.harmonization import EvidenceRunStatus, SVEvidenceSourceArtifact
from hifivar.jasmine import JasmineWrapper
from hifivar.phase9 import run_phase9
from hifivar.truvari import TruvariWrapper


sample = str(snakemake.wildcards.sample)  # type: ignore[name-defined]
config = snakemake.config  # type: ignore[name-defined]
context = AnalysisContext.from_config(config)
section = config["sv"]["harmonization"]
source_paths = [Path(str(item)) for item in snakemake.input.sources]  # type: ignore[name-defined]

if len(source_paths) % 2:
    raise ValueError("Phase 9 workflow inputs must be VCF/index pairs.")

sources = []
for vcf, index in zip(source_paths[::2], source_paths[1::2]):
    assembly_derived = ".assembly.sv.vcf.gz" in vcf.name
    suffix = ".assembly.sv.vcf.gz" if assembly_derived else ".sv.vcf.gz"
    caller = vcf.name.removeprefix(f"{sample}.").removesuffix(suffix)
    sources.append(
        SVEvidenceSourceArtifact(
            sample_id=sample,
            source=SVEvidenceSource.ASSEMBLY if assembly_derived else SVEvidenceSource.READ,
            caller=caller,
            vcf_path=vcf,
            index_path=index,
            status=EvidenceRunStatus.COMPLETED,
            haplotypes=("haplotype1", "haplotype2") if assembly_derived else (),
        )
    )

actual = {(item.source, item.caller) for item in sources}
sv_section = config.get("sv", {})
for caller in ("sawfish", "sniffles2", "pbsv", "cutesv"):
    key = (SVEvidenceSource.READ, caller)
    caller_config = sv_section.get(caller, {})
    if key not in actual and (
        sv_section.get("enabled") is not True
        or not isinstance(caller_config, dict)
        or caller_config.get("enabled", True) is not True
    ):
        sources.append(
            SVEvidenceSourceArtifact(
                sample, SVEvidenceSource.READ, caller, None, None,
                EvidenceRunStatus.DISABLED,
            )
        )

assembly_section = config.get("assembly_sv", {})
for caller in ("pav", "svim_asm"):
    key = (SVEvidenceSource.ASSEMBLY, caller)
    caller_config = assembly_section.get(caller, {})
    if key not in actual and (
        assembly_section.get("enabled") is not True
        or not isinstance(caller_config, dict)
        or caller_config.get("enabled", True) is not True
    ):
        sources.append(
            SVEvidenceSourceArtifact(
                sample, SVEvidenceSource.ASSEMBLY, caller, None, None,
                EvidenceRunStatus.DISABLED,
            )
        )

report = run_phase9(
    context,
    {sample: tuple(sources)},
    output_directory=Path(str(snakemake.output.vcf)).parent.parent,  # type: ignore[name-defined]
    work_directory=Path(str(snakemake.params.work_root)),  # type: ignore[name-defined]
    config=config,
    jasmine_wrapper=JasmineWrapper(
        executable=str(section["jasmine_executable"]),
        bgzip_executable=str(section["bgzip_executable"]),
        tabix_executable=str(section["tabix_executable"]),
    ),
    truvari_wrapper=TruvariWrapper(executable=str(section["truvari_executable"])),
)
report.write_json(
    Path(str(snakemake.output.provenance)),  # type: ignore[name-defined]
    overwrite=bool(section.get("overwrite", False)),
)
