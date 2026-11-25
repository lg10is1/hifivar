"""Opt-in sequential alignment preparation before the existing calling DAG.

This is a two-stage launcher, not a fabricated single Snakemake alignment DAG.
Completed alignments are sealed by checksums; failed attempts remain separate.
"""
from __future__ import annotations

import csv
import io
import json
import os
from copy import deepcopy
from pathlib import Path
from uuid import uuid4

from hifivar.alignment import AlignmentOutputFormat, AlignmentRequest, AlignmentResources, AlignmentTool
from hifivar.alignment_postprocess import AlignmentArtifact, AlignmentIndexFormat, AlignmentIndexRequest, AlignmentSortOrder, AlignmentSource
from hifivar.command import CommandRunner
from hifivar.command_observation import observe_commands
from hifivar.alignment_preparation import publish_preparation, preparation_path, seal_handoff
from hifivar.context import AnalysisContext
from hifivar.exceptions import ConfigurationError, InputValidationError, OutputValidationError
from hifivar.launch_inputs import IndexedBamPreflight, parse_header_records, unique_bam_index
from hifivar.minimap2 import Minimap2Options, Minimap2Wrapper
from hifivar.pbmm2 import Pbmm2Options, Pbmm2Wrapper
from hifivar.sample import InputType
from hifivar.samtools import SamtoolsWrapper
from hifivar.serialization import utc_now_iso8601, write_json_atomic, write_yaml_atomic
from hifivar.ubam import Pbmm2UbamWrapper


def validate_alignment_budget(config: dict) -> None:
    """Reject silent thread/resource changes and unsupported output contracts."""
    a = config["alignment"]
    if a["overwrite"] or a["output_format"] != "bam":
        raise ConfigurationError("Raw launcher requires alignment.overwrite=false and output_format=bam.")
    if a["threads"] > config["runtime"]["threads"] or a["memory_mb"] > config["workflow"]["memory_mb"]:
        raise ConfigurationError("Alignment resource request exceeds launcher budget.")
    if a["index_threads"] + 1 > config["runtime"]["threads"]:
        raise ConfigurationError("alignment.index_threads plus the main thread exceeds runtime.threads.")
    if a["tool"] == "pbmm2" and a["pbmm2_sort_threads"] >= a["threads"]:
        raise ConfigurationError("pbmm2 sort threads must leave at least one mapping thread.")
    if a["tool"] == "minimap2" and a["bam_index_format"] == "csi":
        raise ConfigurationError("This minimap2 wrapper produces BAI only; CSI is not routed.")


class ScopedAlignmentRunner:
    """Delegate all processes to CommandRunner with owned cwd/tmp and an audit log."""

    def __init__(self, runner: CommandRunner, root: Path) -> None:
        self.runner, self.root, self.sequence = runner, root, 0
        self.environment = {k: v for k, v in os.environ.items()
                            if k in {"PATH", "HOME", "USER", "LOGNAME", "LANG", "LC_ALL", "TZ", "CONDA_PREFIX", "CONDA_EXE"}}
        self.environment.update(TMPDIR=str(root), TMP=str(root), TEMP=str(root), PYTHONDONTWRITEBYTECODE="1")

    def require_executable(self, executable):
        return self.runner.require_executable(executable)

    def run(self, command, **kwargs):
        self.sequence += 1
        entry = self.root / f"command_{self.sequence:04d}.json"
        record = {"command": list(map(str, command)), "started": utc_now_iso8601(), "status": "RUNNING"}
        write_json_atomic(record, entry)
        kwargs.update(cwd=self.root, env=self.environment, inherit_env=False)
        try:
            result = self.runner.run(command, **kwargs)
            record.update(status="COMPLETED" if result.returncode == 0 else "FAILED",
                          returncode=result.returncode, duration_seconds=result.duration_seconds)
            for name in ("stdout", "stderr"):
                captured = getattr(result, name, None)
                path = getattr(result, name + "_path", None)
                if path is None and captured is not None:
                    path = entry.with_suffix("." + name + ".log")
                    path.write_text(captured, encoding="utf-8")
                record[name + "_path"] = str(path) if path else None
            return result
        except BaseException as error:
            record.update(status="FAILED", error_type=type(error).__name__, error=str(error))
            raise
        finally:
            record["finished"] = utc_now_iso8601()
            write_json_atomic(record, entry, overwrite=True)


def classify_launch_inputs(context: AnalysisContext, config: dict, runner: CommandRunner,
                          logs: Path) -> list[dict]:
    """Check every input before aligning the first sample."""
    routes = []
    executable = config["workflow"]["samtools_executable"]
    for record in context.samples:
        sample = record.sample
        files = [p.resolve() for p in sample.input.files]
        route = {"sample_id": sample.sample_id, "inputs": list(map(str, files)), "action": "ALIGN", "kind": "FASTQ"}
        if sample.input.input_type is InputType.BAM:
            header = runner.run([executable, "view", "--no-PG", "-H", str(files[0])], timeout=60)
            records = parse_header_records(header.stdout or "")
            if records["@SQ"]:
                check = IndexedBamPreflight(runner, executable).run(files[0], unique_bam_index(files[0]), sample.sample_id, context.reference, logs)
                route.update(action="REUSE", kind="ALIGNED_BAM", validation=check)
            else:
                if config["alignment"]["tool"] != "pbmm2":
                    raise ConfigurationError("Native uBAM routing requires pbmm2; implicit FASTQ conversion is disabled.")
                wrapper = Pbmm2UbamWrapper(runner, config["alignment"]["pbmm2_executable"], executable, _pbmm2_options(config))
                check = wrapper.inspect(files[0], sample.sample_id, timeout=config["alignment"]["runtime_minutes"] * 60)
                route.update(kind="UBAM", validation=check)
        elif sample.input.input_type is not InputType.FASTQ:
            raise InputValidationError("Raw launcher supports FASTQ/uBAM/aligned BAM, not CRAM.")
        routes.append(route)
    return routes


def _pbmm2_options(config: dict) -> Pbmm2Options:
    a = config["alignment"]
    return Pbmm2Options(preset=a["pbmm2_preset"], log_level=a["pbmm2_log_level"], sort_threads=a["pbmm2_sort_threads"])


def prepare_aligned_samples(context: AnalysisContext, config: dict, root: Path,
                            attempt: Path, runner: CommandRunner, *, dry_run: bool,
                            rerun_incomplete: bool) -> Path | None:
    """Publish only fully validated alignment receipts, then materialize calling YAML."""
    from hifivar.launcher import file_identity

    logs = attempt / "alignment_preflight"
    logs.mkdir()
    scoped = ScopedAlignmentRunner(runner, logs)
    routes = classify_launch_inputs(context, config, scoped, logs)
    write_json_atomic({"stage": "ALIGNMENT_PREPARATION", "routes": routes,
                       "calling_dag": "NOT_PLANNED_UNTIL_ALIGNMENTS_EXIST" if dry_run else "PENDING"}, attempt / "alignment_plan.json")
    if dry_run:
        return None
    a = config["alignment"]
    timeout = a["runtime_minutes"] * 60
    output = io.StringIO(newline="")
    writer = csv.writer(output, delimiter="\t", lineterminator="\n")
    writer.writerow(["sample_id", "input", "input_type", "sex", "father", "mother", "phenotype", "group"])
    for record, route in zip(context.samples, routes, strict=True):
        sample = record.sample
        with observe_commands() as commands:
            result, completion, prior = None, None, None
            action = "VALIDATE_EXTERNAL_REUSE" if route["action"] == "REUSE" else "ALIGN_NEW"
            if route["action"] == "REUSE":
                bam = sample.input.files[0].resolve()
                current_logs = attempt / "alignment_validation_logs" / sample.sample_id
                current_logs.mkdir(parents=True)
                local = ScopedAlignmentRunner(runner, current_logs)
                check = IndexedBamPreflight(local, config["workflow"]["samtools_executable"]).run(
                    bam, unique_bam_index(bam), sample.sample_id, context.reference, current_logs)
            else:
                owner = root / "results" / "alignment" / sample.sample_id
                owner.mkdir(parents=True, exist_ok=True)
                receipt = owner / "completed.json"
                if receipt.exists():
                    action, prior = "VALIDATE_PRIOR_ALIGNMENT_REUSE", receipt
                    saved = json.loads(receipt.read_text(encoding="utf-8"))
                    if saved["sample_id"] != sample.sample_id or saved["inputs"] != route["inputs"]:
                        raise InputValidationError("ALIGNMENT_RECEIPT_IDENTITY_MISMATCH")
                    for item in saved["artifacts"]:
                        path = Path(item["path"])
                        if not path.resolve().is_relative_to(owner) or file_identity(path) != item:
                            raise InputValidationError("ALIGNMENT_ARTIFACT_CHANGED: completed BAM/index no longer matches receipt.")
                    bam = Path(saved["bam"])
                    if bam not in [Path(item["path"]) for item in saved["artifacts"]]:
                        raise InputValidationError("ALIGNMENT_RECEIPT_BAM_NOT_SEALED")
                    current_logs = attempt / "alignment_validation_logs" / sample.sample_id
                    current_logs.mkdir(parents=True)
                    local = ScopedAlignmentRunner(runner, current_logs)
                    check = IndexedBamPreflight(local, config["workflow"]["samtools_executable"]).run(
                        bam, unique_bam_index(bam), sample.sample_id, context.reference, current_logs)
                else:
                    if any(owner.iterdir()) and not rerun_incomplete:
                        raise OutputValidationError("Incomplete alignment retained; explicit --resume --rerun-incomplete required.")
                    job = owner / ("attempt_" + uuid4().hex)
                    job.mkdir()
                    local = ScopedAlignmentRunner(runner, job)
                    bam = job / f"{sample.sample_id}.aligned.bam"
                    tool = AlignmentTool(a["tool"])
                    if route["kind"] == "UBAM":
                        result = Pbmm2UbamWrapper(local, a["pbmm2_executable"], config["workflow"]["samtools_executable"], _pbmm2_options(config)).run(
                            sample.input.files[0].resolve(), sample.sample_id, context.reference, bam,
                            threads=a["threads"] - a["pbmm2_sort_threads"], timeout=timeout)
                    else:
                        threads = a["threads"] - a["pbmm2_sort_threads"] if tool is AlignmentTool.PBMM2 else a["threads"]
                        request = AlignmentRequest(sample=sample, reference=context.reference, output_path=bam,
                                                   tool=tool, output_format=AlignmentOutputFormat.BAM,
                                                   resources=AlignmentResources(threads=threads, memory_mb=a["memory_mb"], runtime_minutes=a["runtime_minutes"]))
                        if tool is AlignmentTool.PBMM2:
                            backend = Pbmm2Wrapper(runner=local, executable=a["pbmm2_executable"], options=_pbmm2_options(config))
                        else:
                            if any(c.length >= 2**29 for c in context.reference.contigs):
                                raise ConfigurationError("minimap2 launcher requires a BAI-compatible reference.")
                            backend = Minimap2Wrapper(runner=local, executable=a["minimap2_executable"],
                                                      samtools_executable=config["workflow"]["samtools_executable"], options=Minimap2Options())
                        result = backend.run(request, timeout=timeout, stderr_path=job / "alignment.stderr.log").to_dict()
                    if tool is AlignmentTool.PBMM2:
                        artifact = AlignmentArtifact(sample.sample_id, bam, AlignmentOutputFormat.BAM, context.reference,
                                                     AlignmentSource.GENERATED, AlignmentSortOrder.COORDINATE, tool=tool)
                        selected = None if a["bam_index_format"] == "auto" else AlignmentIndexFormat(a["bam_index_format"])
                        request = AlignmentIndexRequest.create(artifact, index_format=selected, threads=a["index_threads"])
                        index_result = SamtoolsWrapper(runner=local, executable=config["workflow"]["samtools_executable"]).run_index(request, timeout=timeout, stderr_path=job / "index.stderr.log")
                        result["indexing"] = index_result.to_dict()
                    index = unique_bam_index(bam)
                    check = IndexedBamPreflight(local, config["workflow"]["samtools_executable"]).run(bam, index, sample.sample_id, context.reference, job)
                    completion = {"sample_id": sample.sample_id, "inputs": route["inputs"], "kind": route["kind"],
                                       "bam": str(bam), "artifacts": [file_identity(bam), file_identity(index)],
                                  "execution": result, "validation": check}
            publish_preparation(root, attempt, route, bam, unique_bam_index(bam), check, commands,
                                action=action, execution=result, prior=prior)
            if completion is not None:
                completion.update(producer_attempt_id=attempt.name,
                                  producer_receipt=file_identity(preparation_path(attempt, sample.sample_id)))
                write_json_atomic(completion, receipt)
        writer.writerow([sample.sample_id, str(bam), "bam", record.sex.value if record.sex else "",
                         record.father or "", record.mother or "", record.phenotype or "", record.group or ""])
    sheet = root / "calling.samples.tsv"
    effective = deepcopy(config)
    effective["samples"]["sheet"] = str(sheet)
    calling = root / "calling.yaml"
    if sheet.exists():
        if sheet.read_text(encoding="utf-8") != output.getvalue():
            raise InputValidationError("CALLING_SAMPLE_SHEET_CHANGED")
    else:
        with sheet.open("x", encoding="utf-8", newline="") as stream:
            stream.write(output.getvalue())
    if calling.exists():
        import yaml
        try:
            loaded = yaml.safe_load(calling.read_text(encoding="utf-8"))
        except yaml.YAMLError as error:
            raise InputValidationError("CALLING_CONFIG_CHANGED: malformed YAML") from error
        if loaded != effective:
            raise InputValidationError("CALLING_CONFIG_CHANGED")
    else:
        write_yaml_atomic(effective, calling)
    seal_handoff(root, attempt)
    return calling
