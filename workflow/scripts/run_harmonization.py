"""Snakemake bridge for Phase 9 Jasmine harmonization and Truvari comparison."""

from pathlib import Path
import os

from hifivar.assembly_source_roles import source_haplotype_binding
from hifivar.attempt_receipts import ATTEMPT_CONTEXT_ENV

from hifivar.assembly_sv import SVEvidenceSource
from hifivar.attempt_receipts import attempt_step
from hifivar.context import AnalysisContext
from hifivar.exceptions import InputValidationError
from hifivar.harmonization import EvidenceRunStatus, SVEvidenceSourceArtifact
from hifivar.jasmine import JasmineWrapper
from hifivar.phase9 import Phase9Settings, phase9_bnd_paths, phase9_identity_paths, run_phase9
from hifivar.truvari import TruvariWrapper
from hifivar.sv_support_selection import support_selection_paths


sample = str(snakemake.wildcards.sample)  # type: ignore[name-defined]
config = snakemake.config  # type: ignore[name-defined]
context = AnalysisContext.from_config(config)
records = tuple(record for record in context.samples if record.sample.sample_id == sample)
if len(records) != 1:
    raise InputValidationError(f"Phase 9 rule sample '{sample}' is not unique in the sample sheet.")
# A Snakemake wildcard owns exactly one sample, even in a multi-sample run.
# Keep the full sheet/config reference in provenance without requesting results
# for unrelated samples or combining their evidence.
context = AnalysisContext(context.reference, records, context.config, context.source_sample_sheet)
section = config["sv"]["harmonization"]
source_paths = [Path(str(item)) for item in snakemake.input.sources]  # type: ignore[name-defined]

if len(source_paths) % 2:
    raise ValueError("Phase 9 workflow inputs must be VCF/index pairs.")

sources = []
role_bindings = {}
role_context = Path(os.environ[ATTEMPT_CONTEXT_ENV]) if os.environ.get(ATTEMPT_CONTEXT_ENV) else None
external_callers = section.get("input_vcfs", {}).get(sample, {})
for vcf, index in zip(source_paths[::2], source_paths[1::2]):
    assembly_derived = ".assembly.sv.vcf.gz" in vcf.name
    suffix = ".assembly.sv.vcf.gz" if assembly_derived else ".sv.vcf.gz"
    caller = vcf.name.removeprefix(f"{sample}.").removesuffix(suffix)
    role_bindings[caller] = source_haplotype_binding(sample, caller, external=caller in external_callers, context_path=role_context)
    sources.append(
        SVEvidenceSourceArtifact(
            sample_id=sample,
            source=SVEvidenceSource.ASSEMBLY if assembly_derived else SVEvidenceSource.READ,
            caller=caller,
            vcf_path=vcf,
            index_path=index,
            status=EvidenceRunStatus.COMPLETED,
            haplotypes=tuple(role_bindings[caller]["roles"]),
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

receipt_outputs = {
    "vcf": Path(str(snakemake.output.vcf)), "index": Path(str(snakemake.output.index)),
    "evidence": Path(str(snakemake.output.evidence)),
    "provenance": Path(str(snakemake.output.provenance)),
}
settings = Phase9Settings.from_config(config)
if settings.minimum_caller_support > 1:
    support_paths = support_selection_paths(receipt_outputs["vcf"].parent.parent, sample, settings.minimum_caller_support)
    if Path(str(snakemake.output.support_bundle)).absolute() != support_paths["directory"]:
        raise InputValidationError("Caller-support rule output differs from the bound per-sample bundle.")
    receipt_outputs.update({f"support_{role}": path for role, path in support_paths.items() if role != "directory"})
bnd_policy = settings.bnd_policy
if bnd_policy is not None:
    bnd_paths = phase9_bnd_paths(receipt_outputs["vcf"].parent.parent, sample)
    if Path(str(snakemake.output.bnd_bundle)).absolute() != bnd_paths["directory"]:  # type: ignore[name-defined]
        raise InputValidationError("Phase 9 BND rule output differs from the bound per-sample bundle.")
    receipt_outputs.update({f"bnd_{role}": path for role, path in bnd_paths.items() if role != "directory"})
runnable_sources = tuple(item for item in sources if item.vcf_path is not None
                         and item.status in {EvidenceRunStatus.COMPLETED, EvidenceRunStatus.NO_CALLS})
identity_policy = Phase9Settings.from_config(config).input_identity_policy
if identity_policy is not None:
    identity_paths = phase9_identity_paths(receipt_outputs["vcf"].parent.parent, sample,
                                          tuple(item.caller for item in runnable_sources))
    if Path(str(snakemake.output.input_identity_bundle)).absolute() != identity_paths["directory"]:  # type: ignore[name-defined]
        raise InputValidationError("Phase 9 identity rule output differs from the bound per-sample bundle.")
    receipt_outputs.update({f"source_identity_{role}": path for role, path in identity_paths.items() if role != "directory"})
for source in runnable_sources:
    comparison_root = receipt_outputs["vcf"].parent / "truvari" / source.caller
    receipt_outputs[f"comparison_{source.caller}_summary"] = comparison_root / "summary.json"
    for role, filename in (("tp_base", "tp-base.vcf.gz"), ("tp_comp", "tp-comp.vcf.gz"),
                           ("fp", "fp.vcf.gz"), ("fn", "fn.vcf.gz")):
        receipt_outputs[f"comparison_{source.caller}_{role}"] = comparison_root / filename
with attempt_step("sv.harmonization", (sample,), outputs=receipt_outputs,
                  inputs={item.caller: {"vcf": item.vcf_path, "index": item.index_path}
                          for item in runnable_sources},
                  input_manifest=Path(config["paths"]["outdir"]).parent / "effective.yaml") as receipt:
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
    sample_result = report.sample_results[0]
    selected_result = getattr(sample_result, "support_selection", None)
    if settings.minimum_caller_support > 1 and selected_result is None:
        raise InputValidationError("Configured caller-support selection did not produce its bound result.")
    versions = dict(sample_result.jasmine.tool_versions)
    versions.update({f"truvari_{source.caller}": comparison.version
                     for source, comparison in zip(runnable_sources, sample_result.truvari)})
    receipt.complete(
        outputs=receipt_outputs, tool_versions=versions,
        validation={"validated": True, "scope": "PHASE9_WRAPPER_HARMONIZED_VCF_HEADER_INDEX_EVIDENCE_MEMBERSHIP_AND_PER_SOURCE_TRUVARI_OUTPUTS_NOT_TRUTH"},
        metadata={"reference_build": context.reference.build, "reference_build_source": "EXPLICIT_REFERENCE_CONTEXT",
                  "source_order": [item.caller for item in runnable_sources],
                  "source_haplotype_roles": role_bindings,
                  "bnd_representation": bnd_policy.to_dict() if bnd_policy else {"enabled": False},
                  "input_identity": identity_policy.to_dict() if identity_policy else {"enabled": False},
                  "canonical_source_files": [str(item.vcf_path) for item in runnable_sources],
                  "analysis_source_files": [str(path) for path in sample_result.jasmine.source_identity.views] if identity_policy else [],
                  "minimum_caller_support": settings.minimum_caller_support,
                  "support_selection": selected_result.to_dict() if selected_result is not None else None,
                  "native_output_role_unchanged": True,
                  "harmonization_is_truth": False, "caller_count_is_confidence": False,
                  "tool_version_sources": {tool: "COMMAND_PROBED" for tool in versions}},
    )
