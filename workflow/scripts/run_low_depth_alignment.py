"""Align one Phase 15 raw FASTQ subset and create its explicit BAM index."""

from pathlib import Path

from hifivar.alignment import AlignmentOutputFormat, AlignmentRequest, AlignmentResources, AlignmentTool
from hifivar.alignment_postprocess import AlignmentArtifact, AlignmentIndexFormat, AlignmentIndexRequest
from hifivar.pbmm2 import Pbmm2Options, Pbmm2Wrapper
from hifivar.reference import ReferenceGenome
from hifivar.sample import InputDataset, Sample
from hifivar.samtools import SamtoolsWrapper
from hifivar.serialization import write_json_atomic

section = snakemake.config["low_depth"]  # type: ignore[name-defined]
alignment = snakemake.config["alignment"]  # type: ignore[name-defined]
sample = Sample(
    str(section["sample_id"]),
    InputDataset.from_files((Path(str(snakemake.input.fastq)),)),  # type: ignore[name-defined]
)
reference = ReferenceGenome.from_fasta(
    Path(str(snakemake.input.reference)),  # type: ignore[name-defined]
    build=snakemake.config["reference"]["build"],  # type: ignore[name-defined]
)
request = AlignmentRequest(
    sample=sample,
    reference=reference,
    output_path=Path(str(snakemake.output.bam)),  # type: ignore[name-defined]
    tool=AlignmentTool.PBMM2,
    output_format=AlignmentOutputFormat.BAM,
    resources=AlignmentResources(
        threads=int(snakemake.threads),  # type: ignore[name-defined]
        memory_mb=int(snakemake.resources.mem_mb),  # type: ignore[name-defined]
        runtime_minutes=int(snakemake.resources.runtime_min),  # type: ignore[name-defined]
    ),
    overwrite=bool(section.get("overwrite", False)),
)
alignment_result = Pbmm2Wrapper(
    executable=str(section["end_to_end"].get("pbmm2_executable", "pbmm2")),
    options=Pbmm2Options.from_config(snakemake.config),  # type: ignore[name-defined]
).run(request, stderr_path=Path(str(snakemake.log[0])))  # type: ignore[name-defined]
artifact = AlignmentArtifact.from_result(alignment_result)
index_request = AlignmentIndexRequest(
    artifact=artifact,
    index_format=AlignmentIndexFormat.BAI,
    output_path=Path(str(snakemake.output.index)),  # type: ignore[name-defined]
    threads=int(alignment.get("index_threads", 1)),
    overwrite=bool(section.get("overwrite", False)),
)
index_result = SamtoolsWrapper(
    executable=str(section["end_to_end"].get("samtools_executable", "samtools")),
).run_index(index_request, stderr_path=Path(str(snakemake.log[0])))  # type: ignore[name-defined]
write_json_atomic(
    {"alignment": alignment_result.to_dict(), "index": index_result.to_dict()},
    Path(str(snakemake.output.provenance)),  # type: ignore[name-defined]
    overwrite=bool(section.get("overwrite", False)),
    artifact_name="Low-depth alignment provenance",
)
