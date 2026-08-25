"""Snakemake bridge delegating Phase 8 PAV execution to PavWrapper."""

from pathlib import Path

from hifivar.assembly import AssemblyRole, HaplotypeAssemblyArtifact
from hifivar.assembly_sv import AssemblySvCaller, AssemblySvRequest, AssemblySvResources
from hifivar.pav import PavWrapper
from hifivar.reference import ReferenceGenome

sample = str(snakemake.wildcards.sample)  # type: ignore[name-defined]
section = snakemake.config["assembly_sv"]  # type: ignore[name-defined]
caller = section["pav"]
reference = ReferenceGenome.from_fasta(
    Path(str(snakemake.input.reference)),  # type: ignore[name-defined]
    build=snakemake.config["reference"].get("build"),  # type: ignore[name-defined]
)
assemblies = tuple(
    HaplotypeAssemblyArtifact(
        sample, role, path, path, "phase7_workflow_handoff", path.stat().st_size
    )
    for role, path in (
        (AssemblyRole.HAPLOTYPE1, Path(str(snakemake.input.hap1))),  # type: ignore[name-defined]
        (AssemblyRole.HAPLOTYPE2, Path(str(snakemake.input.hap2))),  # type: ignore[name-defined]
    )
)
request = AssemblySvRequest(
    sample, AssemblySvCaller.PAV, reference, assemblies,
    Path(str(snakemake.params.work)), Path(str(snakemake.output.vcf)),  # type: ignore[name-defined]
    AssemblySvResources(int(snakemake.threads), int(snakemake.resources.mem_mb), int(snakemake.resources.runtime_min)),  # type: ignore[name-defined]
    bool(section.get("overwrite", False)),
)
PavWrapper(
    snakefile=str(caller["snakefile"]), executable=str(caller.get("executable", "snakemake")),
    pav_version=str(caller.get("version", "VERSION_PENDING_LINUX_VERIFICATION")),
).run(request, stderr_path=Path(str(snakemake.log[0])))  # type: ignore[name-defined]
