"""Snakemake bridge delegating Phase 7 exclusively to HifiasmWrapper."""

from pathlib import Path

from hifivar.assembly import AssemblyRequest, AssemblyResources, AssemblyRole
from hifivar.attempt_receipts import attempt_step
from hifivar.hifiasm import HifiasmWrapper
from hifivar.sample import InputDataset, Sample


sample_id = str(snakemake.wildcards.sample)  # type: ignore[name-defined]
sample = Sample(
    sample_id,
    InputDataset.from_files(
        [Path(str(path)) for path in snakemake.input.reads]  # type: ignore[name-defined]
    ),
)
section = snakemake.config["assembly"]  # type: ignore[name-defined]
request = AssemblyRequest(
    sample,
    Path(str(snakemake.params.output_prefix)),  # type: ignore[name-defined]
    Path(str(snakemake.params.assembly_directory)),  # type: ignore[name-defined]
    AssemblyResources(
        threads=int(snakemake.threads),  # type: ignore[name-defined]
        memory_mb=int(snakemake.resources.mem_mb),  # type: ignore[name-defined]
        runtime_minutes=int(snakemake.resources.runtime_min),  # type: ignore[name-defined]
    ),
    overwrite=bool(section.get("overwrite", False)),
)
wrapper = HifiasmWrapper(executable=str(section.get("executable", "hifiasm")))
receipt_outputs = {
    "primary_fasta": request.fasta_paths[AssemblyRole.PRIMARY],
    "haplotype1_fasta": request.fasta_paths[AssemblyRole.HAPLOTYPE1],
    "haplotype2_fasta": request.fasta_paths[AssemblyRole.HAPLOTYPE2],
}
with attempt_step("assembly", (sample_id,), outputs=receipt_outputs) as receipt:
    result = wrapper.run(request, stderr_path=Path(str(snakemake.log[0])))  # type: ignore[name-defined]
    receipt.complete(
        outputs=receipt_outputs, tool_versions={"hifiasm": result.hifiasm_version},
        validation={"validated": True, "scope": "WRAPPER_REQUIRED_GFA_SEGMENTS_STREAM_CONVERSION_AND_DERIVED_FASTA_VALIDATION_NOT_ASSEMBLY_ACCURACY"},
        metadata={"reference_independent": True, "tool_version_sources": {"hifiasm": "COMMAND_PROBED"},
                  "assembly_artifacts": [item.to_dict() for item in result.artifact.assemblies]},
    )
