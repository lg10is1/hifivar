"""Snakemake bridge delegating Phase 8 SVIM-asm to SvimAsmWrapper."""

from pathlib import Path

from hifivar.assembly import AssemblyRole, HaplotypeAssemblyArtifact
from hifivar.assembly_sv import AssemblySvCaller, AssemblySvRequest, AssemblySvResources
from hifivar.reference import ReferenceGenome
from hifivar.svim_asm import SvimAsmWrapper

sample = str(snakemake.wildcards.sample)  # type: ignore[name-defined]
section = snakemake.config["assembly_sv"]  # type: ignore[name-defined]
caller = section["svim_asm"]
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
    sample, AssemblySvCaller.SVIM_ASM, reference, assemblies,
    Path(str(snakemake.params.work)), Path(str(snakemake.output.vcf)),  # type: ignore[name-defined]
    AssemblySvResources(int(snakemake.threads), int(snakemake.resources.mem_mb), int(snakemake.resources.runtime_min)),  # type: ignore[name-defined]
    bool(section.get("overwrite", False)),
)
SvimAsmWrapper(
    executable=str(caller.get("executable", "svim-asm")),
    minimap2_executable=str(caller.get("minimap2_executable", "minimap2")),
    samtools_executable=str(caller.get("samtools_executable", "samtools")),
    bgzip_executable=str(caller.get("bgzip_executable", "bgzip")),
    tabix_executable=str(caller.get("tabix_executable", "tabix")),
).run(request, stderr_path=Path(str(snakemake.log[0])))  # type: ignore[name-defined]
