"""Snakemake bridge delegating Phase 8 SVIM-asm to SvimAsmWrapper."""

from pathlib import Path

from hifivar.assembly_reuse import workflow_assemblies
from hifivar.assembly_sv import AssemblySvCaller, AssemblySvRequest, AssemblySvResources
from hifivar.attempt_receipts import attempt_step
from hifivar.assembly_source_roles import assembly_role_binding
from hifivar.exceptions import InputValidationError
from hifivar.reference import ReferenceGenome
from hifivar.svim_asm import SvimAsmWrapper

sample = str(snakemake.wildcards.sample)  # type: ignore[name-defined]
section = snakemake.config["assembly_sv"]  # type: ignore[name-defined]
caller = section["svim_asm"]
reference = ReferenceGenome.from_fasta(
    Path(str(snakemake.input.reference)),  # type: ignore[name-defined]
    build=snakemake.config["reference"].get("build"),  # type: ignore[name-defined]
)
input_paths = ([Path(str(path)) for path in snakemake.input.assemblies]  # type: ignore[name-defined]
               if hasattr(snakemake.input, "assemblies") else
               [Path(str(snakemake.input.hap1)), Path(str(snakemake.input.hap2))])  # type: ignore[name-defined]
assemblies = workflow_assemblies(sample, input_paths, snakemake.config)  # type: ignore[name-defined]
request = AssemblySvRequest(
    sample, AssemblySvCaller.SVIM_ASM, reference, assemblies,
    Path(str(snakemake.params.work)), Path(str(snakemake.output.vcf)),  # type: ignore[name-defined]
    AssemblySvResources(int(snakemake.threads), int(snakemake.resources.mem_mb), int(snakemake.resources.runtime_min)),  # type: ignore[name-defined]
    bool(section.get("overwrite", False)),
)
wrapper = SvimAsmWrapper(
    executable=str(caller.get("executable", "svim-asm")),
    minimap2_executable=str(caller.get("minimap2_executable", "minimap2")),
    samtools_executable=str(caller.get("samtools_executable", "samtools")),
    bgzip_executable=str(caller.get("bgzip_executable", "bgzip")),
    tabix_executable=str(caller.get("tabix_executable", "tabix")),
)
receipt_outputs = {"vcf": request.output_vcf, "index": request.output_index}
with attempt_step("assembly_sv.svim_asm", (sample,), outputs=receipt_outputs) as receipt:
    input_role_binding = assembly_role_binding(sample, assemblies, receipt.context_path)
    result = wrapper.run(request, stderr_path=Path(str(snakemake.log[0])))  # type: ignore[name-defined]
    if assembly_role_binding(sample, assemblies, receipt.context_path) != input_role_binding:
        raise InputValidationError("Assembly request inputs or upstream role receipt changed during execution.")
    receipt.complete(
        outputs=receipt_outputs, tool_versions=result.tool_versions,
        validation={"validated": True, "scope": "WRAPPER_HAPLOTYPE_ALIGNMENT_SORT_INDEX_AND_FINAL_SV_VCF_HEADER_BGZF_TBI_VALIDATION"},
        metadata={"reference_build": reference.build, "reference_build_source": "EXPLICIT_REFERENCE_CONTEXT",
                  "evidence_source": "ASSEMBLY", "assembly_role_binding": input_role_binding, "tool_version_sources": {tool: "COMMAND_PROBED" for tool in result.tool_versions}},
    )
