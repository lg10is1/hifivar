"""Snakemake bridge delegating Phase 8 PAV execution to PavWrapper."""

from pathlib import Path

from hifivar.assembly_reuse import workflow_assemblies
from hifivar.assembly_sv import AssemblySvCaller, AssemblySvRequest, AssemblySvResources
from hifivar.attempt_receipts import attempt_step
from hifivar.assembly_source_roles import assembly_role_binding
from hifivar.exceptions import InputValidationError
from hifivar.pav import PavWrapper
from hifivar.reference import ReferenceGenome

sample = str(snakemake.wildcards.sample)  # type: ignore[name-defined]
section = snakemake.config["assembly_sv"]  # type: ignore[name-defined]
caller = section["pav"]
reference = ReferenceGenome.from_fasta(
    Path(str(snakemake.input.reference)),  # type: ignore[name-defined]
    build=snakemake.config["reference"].get("build"),  # type: ignore[name-defined]
)
input_paths = ([Path(str(path)) for path in snakemake.input.assemblies]  # type: ignore[name-defined]
               if hasattr(snakemake.input, "assemblies") else
               [Path(str(snakemake.input.hap1)), Path(str(snakemake.input.hap2))])  # type: ignore[name-defined]
assemblies = workflow_assemblies(sample, input_paths, snakemake.config)  # type: ignore[name-defined]
request = AssemblySvRequest(
    sample, AssemblySvCaller.PAV, reference, assemblies,
    Path(str(snakemake.params.work)), Path(str(snakemake.output.vcf)),  # type: ignore[name-defined]
    AssemblySvResources(int(snakemake.threads), int(snakemake.resources.mem_mb), int(snakemake.resources.runtime_min)),  # type: ignore[name-defined]
    bool(section.get("overwrite", False)),
)
wrapper = PavWrapper(
    snakefile=str(caller["snakefile"]), executable=str(caller.get("executable", "snakemake")),
    bgzip_executable=str(caller.get("bgzip_executable", "bgzip")),
    tabix_executable=str(caller.get("tabix_executable", "tabix")),
    pav_version=str(caller.get("version", "VERSION_PENDING_LINUX_VERIFICATION")),
)
receipt_outputs = {"vcf": request.output_vcf, "index": request.output_index}
with attempt_step("assembly_sv.pav", (sample,), outputs=receipt_outputs) as receipt:
    input_role_binding = assembly_role_binding(sample, assemblies, receipt.context_path)
    result = wrapper.run(request, stderr_path=Path(str(snakemake.log[0])))  # type: ignore[name-defined]
    if assembly_role_binding(sample, assemblies, receipt.context_path) != input_role_binding:
        raise InputValidationError("Assembly request inputs or upstream role receipt changed during execution.")
    versions = {"pav": result.pav_version, "snakemake": result.adapter_version,
                **result.finalizer_versions}
    receipt.complete(
        outputs=receipt_outputs, tool_versions=versions,
        validation={"validated": True, "scope": "WRAPPER_PAV_2_4_6_SV_SELECTION_FIELD_PRESERVATION_BGZF_TBI_AND_HEADER_VALIDATION"},
        metadata={"reference_build": reference.build, "reference_build_source": "EXPLICIT_REFERENCE_CONTEXT",
                  "evidence_source": "ASSEMBLY", "assembly_role_binding": input_role_binding, "pav_version_source": result.pav_version_source,
                  "tool_version_sources": {**{tool: "COMMAND_PROBED" for tool in versions},
                                           "pav": "USER_DECLARED_NOT_NATIVE_PROBED"}},
    )
