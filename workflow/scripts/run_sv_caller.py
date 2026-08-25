from pathlib import Path

from hifivar.alignment import AlignmentOutputFormat
from hifivar.alignment_postprocess import AlignmentArtifact, AlignmentSortOrder, AlignmentSource
from hifivar.cutesv import CuteSvRequest, CuteSvResources, CuteSvWrapper
from hifivar.pbsv import PbsvRequest, PbsvResources, PbsvWrapper
from hifivar.reference import ReferenceGenome
from hifivar.sawfish import SawfishRequest, SawfishResources, SawfishWrapper
from hifivar.sniffles2 import Sniffles2Request, Sniffles2Resources, Sniffles2Wrapper
from hifivar.sv import BgzipTabixWrapper, SvCaller, VcfFinalizeRequest, create_structural_variant_artifact


caller = str(snakemake.wildcards.caller)  # type: ignore[name-defined]
sample = str(snakemake.wildcards.sample)  # type: ignore[name-defined]
alignment_path = Path(str(snakemake.input.alignment))  # type: ignore[name-defined]
reference = ReferenceGenome.from_fasta(
    Path(str(snakemake.input.reference)),  # type: ignore[name-defined]
    build=snakemake.config["reference"].get("build"),  # type: ignore[name-defined]
)
artifact = AlignmentArtifact(
    sample_id=sample,
    path=alignment_path,
    output_format=(AlignmentOutputFormat.BAM if alignment_path.suffix.lower() == ".bam" else AlignmentOutputFormat.CRAM),
    reference=reference,
    source=AlignmentSource.EXISTING,
    sort_order=AlignmentSortOrder.UNKNOWN,
    index_path=Path(str(snakemake.input.alignment_index)),  # type: ignore[name-defined]
)
sv_config = snakemake.config["sv"]  # type: ignore[name-defined]
caller_config = sv_config[caller]
overwrite = bool(sv_config.get("overwrite", False))
output_vcf = Path(str(snakemake.output.vcf))  # type: ignore[name-defined]
work = Path(str(snakemake.params.workdir))  # type: ignore[name-defined]
log = Path(str(snakemake.log[0]))  # type: ignore[name-defined]
threads = int(snakemake.threads)  # type: ignore[name-defined]
memory_mb = int(snakemake.resources.mem_mb)  # type: ignore[name-defined]
runtime_minutes = int(snakemake.resources.runtime_min)  # type: ignore[name-defined]

if caller == "sawfish":
    request = SawfishRequest(
        artifact, output_vcf, work,
        SawfishResources(threads, memory_mb, runtime_minutes),
        overwrite, bool(caller_config.get("disable_cnv", False)),
    )
    result = SawfishWrapper(executable=str(caller_config["executable"])).run(request, stderr_path=log)
    create_structural_variant_artifact(caller=SvCaller.SAWFISH, sample_id=sample, reference=reference, vcf_path=output_vcf, caller_version=result.tool_version, commands=tuple(command.args for command in result.commands))
elif caller == "sniffles2":
    request = Sniffles2Request(
        artifact, output_vcf, Sniffles2Resources(threads, memory_mb, runtime_minutes),
        caller_config.get("minimum_support"), int(caller_config["minimum_sv_length"]), overwrite,
    )
    result = Sniffles2Wrapper(executable=str(caller_config["executable"])).run(request, stderr_path=log)
    create_structural_variant_artifact(caller=SvCaller.SNIFFLES2, sample_id=sample, reference=reference, vcf_path=output_vcf, caller_version=result.tool_version, commands=(result.command.args,))
elif caller == "pbsv":
    request = PbsvRequest(
        artifact, work / f"{sample}.pbsv.svsig.gz", work / f"{sample}.pbsv.raw.vcf",
        PbsvResources(threads, memory_mb, runtime_minutes), overwrite,
    )
    result = PbsvWrapper(executable=str(caller_config["executable"])).run(request, stderr_path=log)
    final_config = sv_config["finalization"]
    BgzipTabixWrapper(bgzip_executable=str(final_config["bgzip_executable"]), tabix_executable=str(final_config["tabix_executable"])).run(
        VcfFinalizeRequest(SvCaller.PBSV, sample, reference, request.raw_vcf, output_vcf, result.tool_version, tuple(command.args for command in result.commands), overwrite), stderr_path=log
    )
elif caller == "cutesv":
    request = CuteSvRequest(
        artifact=artifact,
        raw_vcf=work / f"{sample}.cutesv.raw.vcf",
        work_directory=work / "tmp",
        resources=CuteSvResources(threads, memory_mb, runtime_minutes),
        minimum_support=int(caller_config["minimum_support"]),
        minimum_sv_size=int(caller_config["minimum_sv_size"]),
        max_cluster_bias_ins=int(caller_config["max_cluster_bias_ins"]),
        diff_ratio_merging_ins=float(caller_config["diff_ratio_merging_ins"]),
        max_cluster_bias_del=int(caller_config["max_cluster_bias_del"]),
        diff_ratio_merging_del=float(caller_config["diff_ratio_merging_del"]),
        genotype=bool(caller_config["genotype"]),
        overwrite=overwrite,
    )
    result = CuteSvWrapper(executable=str(caller_config["executable"])).run(request, stderr_path=log)
    final_config = sv_config["finalization"]
    BgzipTabixWrapper(bgzip_executable=str(final_config["bgzip_executable"]), tabix_executable=str(final_config["tabix_executable"])).run(
        VcfFinalizeRequest(SvCaller.CUTESV, sample, reference, request.raw_vcf, output_vcf, result.tool_version, (result.command.args,), overwrite), stderr_path=log
    )
else:  # pragma: no cover
    raise ValueError(f"Unsupported Phase 4 caller: {caller}")
