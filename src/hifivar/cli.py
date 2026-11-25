"""Command-line interface foundation for HiFiVar."""

from __future__ import annotations

import argparse
import json
import platform
import sys
from collections.abc import Mapping, Sequence
from contextlib import ExitStack, contextmanager
from importlib.resources import as_file, files
from pathlib import Path
from typing import Iterator

import yaml

from hifivar import __version__
from hifivar.config import (
    HiFiVarConfig,
    load_config,
    validate_config,
    write_effective_config,
)
from hifivar.exceptions import ConfigurationError, HiFiVarError
from hifivar.logging_utils import (
    SUPPORTED_LOG_LEVELS,
    configure_logging,
    get_logger,
    parse_log_level,
)
from hifivar.serialization import redact_sensitive_data


_LOGGER = get_logger(__name__)
_PRESET_NAMES = ("fast", "standard", "comprehensive", "cohort", "trio")


def build_parser() -> argparse.ArgumentParser:
    """Build and return the testable top-level argument parser."""
    parser = argparse.ArgumentParser(
        prog="hifivar",
        description="PacBio HiFi whole-genome variant analysis framework",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )
    parser.add_argument(
        "--config",
        dest="user_config",
        type=Path,
        metavar="CONFIG.yaml",
        help="user YAML configuration overriding preset and defaults",
    )
    parser.add_argument(
        "--preset",
        choices=_PRESET_NAMES,
        default="standard",
        help="packaged workflow preset (default: %(default)s)",
    )
    parser.add_argument(
        "--log-level",
        type=_cli_log_level,
        metavar="LEVEL",
        help="override logging level",
    )
    parser.add_argument(
        "--log-file",
        type=Path,
        metavar="PATH",
        help="override UTF-8 HiFiVar log file",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="plan workflow execution without running analysis jobs",
    )

    subparsers = parser.add_subparsers(dest="command", metavar="COMMAND")

    config_parser = subparsers.add_parser(
        "config",
        help="inspect and validate effective configuration",
    )
    config_subparsers = config_parser.add_subparsers(
        dest="config_command",
        metavar="ACTION",
        required=True,
    )
    config_subparsers.add_parser("show", help="write effective YAML to stdout")
    config_subparsers.add_parser("validate", help="validate effective config")
    dump_parser = config_subparsers.add_parser(
        "dump-effective",
        help="write effective YAML to a file",
    )
    dump_parser.add_argument(
        "--output",
        required=True,
        type=Path,
        metavar="PATH",
        help="destination effective YAML file",
    )

    doctor_parser = subparsers.add_parser(
        "doctor",
        help="report the local HiFiVar foundation environment",
    )
    doctor_parser.add_argument("--preflight", action="store_true", help="check the fresh calling-launcher scope without running analysis")
    doctor_parser.add_argument("--run-dir", type=Path, help="intended new calling directory; preflight does not create it")
    doctor_parser.add_argument("--probe-tools", action="store_true", help="explicitly allow bounded native version/BAM metadata probes")
    doctor_parser.add_argument("--json", action="store_true", help="print a machine-readable preflight report")
    doctor_parser.add_argument("--timeout-seconds", type=float, default=30, help="per-probe timeout (default: 30)")
    run_parser = subparsers.add_parser(
        "run", help="launch batch calling; optionally prepare FASTQ/uBAM with alignment.enabled",
    )
    run_parser.add_argument("--run-dir", required=True, type=Path)
    run_parser.add_argument("--resume", action="store_true", help="resume an identity-matched launcher directory")
    run_parser.add_argument("--rerun-incomplete", action="store_true", help="explicitly retry incomplete jobs with --resume")
    review_parser = subparsers.add_parser("review", help="prepare IGV review bundles or import explicit manual decisions")
    review_actions = review_parser.add_subparsers(dest="review_command", metavar="ACTION", required=True)
    prepare = review_actions.add_parser("prepare", help="write an IGV batch and pending manifest without executing IGV")
    prepare.add_argument("--output-dir", required=True, type=Path)
    decisions = review_actions.add_parser("import-decisions", help="store a supplied reviewer TSV as independent evidence")
    decisions.add_argument("--bundle-dir", required=True, type=Path)
    decisions.add_argument("--decisions", required=True, type=Path)
    decisions.add_argument("--output-dir", required=True, type=Path)
    snapshots = review_actions.add_parser("snapshots", help="execute IGV for explicit existing external source-plan candidates; human review remains pending")
    snapshots.add_argument("--output-dir", required=True, type=Path)
    snapshots.add_argument("--settings", type=Path, help="optional bounded YAML for view windows, panel height, delay and timeout")
    tr_plot = review_actions.add_parser("tr-plot", help="plot one explicitly selected existing TRGT locus; no calling or human review")
    tr_plot.add_argument("--plan", required=True, type=Path, metavar="PLAN.yaml", help="strict external-input plotting YAML plan")
    benchmark_parser = subparsers.add_parser("benchmark", help="benchmark explicitly listed existing variant results")
    benchmark_actions = benchmark_parser.add_subparsers(dest="benchmark_command", metavar="ACTION", required=True)
    sv_benchmark = benchmark_actions.add_parser("sv", help="sequential, truth-bound SV callset benchmarks without calling")
    sv_benchmark.add_argument("--manifest", required=True, type=Path)
    sv_benchmark.add_argument("--output-dir", required=True, type=Path)
    sv_depth = benchmark_actions.add_parser("sv-depth", help="call and benchmark SVs from an explicit existing depth/replicate BAM manifest")
    sv_depth.add_argument("--manifest", required=True, type=Path)
    sv_depth.add_argument("--output-dir", required=True, type=Path)
    sv_depth.add_argument("--unit-id", action="append", dest="unit_ids", help="limit execution to declared unit; repeatable")
    sv_depth.add_argument("--caller", action="append", dest="callers", choices=("sawfish", "sniffles2", "pbsv", "cutesv"))
    sv_depth.add_argument("--resume", action="store_true", help="reuse only verified complete receipts; failed outputs are never retried in place")
    sv_depth_support = benchmark_actions.add_parser("sv-depth-support", help="compare caller-support thresholds from four complete read-SV depth outputs; no calling")
    sv_depth_support.add_argument("--manifest", required=True, type=Path)
    sv_depth_support.add_argument("--output-dir", required=True, type=Path)
    sv_depth_support.add_argument("--unit-id", action="append", dest="unit_ids", help="limit analysis to an explicitly mapped unit; repeatable")
    functional = subparsers.add_parser("functional", help="predict functional effects for explicitly selected variants; not call confidence")
    functional_actions = functional.add_subparsers(dest="functional_command", required=True)
    alpha = functional_actions.add_parser("alphagenome", help="opt-in official AlphaGenome API using the user's own environment credential")
    alpha.add_argument("--plan", required=True, type=Path)
    alpha.add_argument("--output-dir", required=True, type=Path)
    alpha.add_argument("--allow-network", action="store_true", help="explicitly allow sending the selected genomic variants to the official service")
    report_parser = subparsers.add_parser("report", help="generate local reports for an explicitly selected existing launcher attempt")
    report_parser.add_argument("--run-dir", required=True, type=Path)
    report_parser.add_argument("--output-dir", required=True, type=Path)
    report_parser.add_argument("--attempt-id", help="required when the original run contains multiple attempts")
    primers = subparsers.add_parser("primers", help="prepare local experimental-validation candidates for sequence-resolved variants")
    primer_actions = primers.add_subparsers(dest="primer_command", required=True)
    primer_design = primer_actions.add_parser("design", help="design PCR/Sanger or simple DEL/INS flanking candidates; no experimental success claim")
    primer_design.add_argument("--plan", required=True, type=Path)
    primer_design.add_argument("--output-dir", required=True, type=Path)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Run the HiFiVar CLI and return a process exit code."""
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command is None:
        parser.print_help()
        return 0

    try:
        with _packaged_config_paths(args.preset) as (
            default_config,
            preset_config,
        ):
            config = load_config(
                default_config=default_config,
                preset=preset_config,
                user_config=args.user_config,
            )
        config = _apply_cli_overrides(config, args)
        if args.command in {"run", "review", "benchmark", "report", "functional", "primers"} or (args.command == "doctor" and (
            args.preflight or args.run_dir is not None or args.probe_tools or args.json or args.timeout_seconds != 30
        )):
            # Do not open a user-configured log outside the owned run directory.
            configure_logging(level=config["logging"]["level"])
        else:
            _configure_cli_logging(config)
        _LOGGER.info(
            "HiFiVar CLI initialized (command=%s, dry_run=%s)",
            args.command,
            args.dry_run,
        )
        return _dispatch(args, config)
    except HiFiVarError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1


def _cli_log_level(value: str) -> str:
    """Validate and normalize one argparse logging-level value."""
    try:
        parse_log_level(value)
    except ConfigurationError as error:
        supported = ", ".join(SUPPORTED_LOG_LEVELS)
        raise argparse.ArgumentTypeError(
            f"invalid logging level {value!r}; choose from {supported}"
        ) from error
    return value.strip().upper()


@contextmanager
def _packaged_config_paths(preset: str) -> Iterator[tuple[Path, Path]]:
    """Materialize packaged default and preset resources for ``load_config``."""
    config_root = files("hifivar").joinpath("resources", "configs")
    default_resource = config_root.joinpath("default.yaml")
    preset_resource = config_root.joinpath("presets", f"{preset}.yaml")
    if not default_resource.is_file():
        raise ConfigurationError(
            "Packaged default configuration resource is missing: default.yaml"
        )
    if not preset_resource.is_file():
        raise ConfigurationError(
            f"Packaged preset configuration resource is missing: {preset}.yaml"
        )

    with ExitStack() as stack:
        default_path = Path(stack.enter_context(as_file(default_resource)))
        preset_path = Path(stack.enter_context(as_file(preset_resource)))
        yield default_path, preset_path


def _apply_cli_overrides(
    config: HiFiVarConfig,
    args: argparse.Namespace,
) -> HiFiVarConfig:
    """Apply explicit CLI settings without changing any source YAML."""
    effective = config.to_dict()
    logging_config = effective.get("logging")
    if not isinstance(logging_config, dict):
        raise ConfigurationError("Effective logging config must be a mapping.")

    if args.log_level is not None:
        logging_config["level"] = args.log_level
    if args.log_file is not None:
        logging_config["file"] = str(args.log_file.expanduser())

    validate_config(effective, require_complete=True)
    return HiFiVarConfig(effective, config.sources)


def _configure_cli_logging(config: HiFiVarConfig) -> None:
    """Configure logging after config loading and CLI override application."""
    logging_config = config["logging"]
    if not isinstance(logging_config, Mapping):
        raise ConfigurationError("Effective logging config must be a mapping.")

    level = logging_config["level"]
    log_file = logging_config["file"]
    if not isinstance(level, str):
        raise ConfigurationError("logging.level must be a string.")
    if log_file is not None and not isinstance(log_file, str):
        raise ConfigurationError("logging.file must be a string path or null.")

    try:
        configure_logging(
            level=level,
            log_file=Path(log_file) if log_file is not None else None,
        )
    except OSError as error:
        raise ConfigurationError(
            f"Unable to configure CLI logging: {error}"
        ) from error


def _dispatch(args: argparse.Namespace, config: HiFiVarConfig) -> int:
    """Dispatch the registered configuration, calling and review entries."""
    if args.command == "config":
        return _dispatch_config(args, config)
    if args.command == "doctor":
        if args.preflight:
            if args.run_dir is None:
                raise ConfigurationError("doctor --preflight requires --run-dir for a fresh calling run.")
            if args.dry_run and args.probe_tools:
                raise ConfigurationError("--dry-run cannot be combined with --probe-tools.")
            from hifivar.preflight import run_preflight

            report = run_preflight(config, args.run_dir, probe_tools=args.probe_tools,
                                   timeout_seconds=args.timeout_seconds)
            if args.json:
                print(json.dumps(report, ensure_ascii=False))
            else:
                print(f"HiFiVar calling preflight: {report['status']}")
                for check in report["checks"]:
                    print(f"{check['status']}: {check['key']} — {check['message']}")
                print("Scope only; no analysis, universal readiness or accuracy claim.")
            return report["exit_code"]
        if args.run_dir is not None or args.probe_tools or args.json or args.timeout_seconds != 30:
            raise ConfigurationError("Preflight options require doctor --preflight.")
        return _run_doctor()
    if args.command == "run":
        from hifivar.launcher import WorkflowLauncher

        outcome = WorkflowLauncher().run(
            config, run_dir=args.run_dir, resume=args.resume,
            rerun_incomplete=args.rerun_incomplete, dry_run=args.dry_run,
        )
        print(f"{outcome}: {args.run_dir.resolve()}")
        return 0
    if args.command == "review":
        if args.review_command == "prepare":
            from hifivar.review_cli import prepare_review

            result = prepare_review(config, args.output_dir, dry_run=args.dry_run)
        elif args.review_command == "import-decisions":
            from hifivar.review_cli import import_review_decisions

            result = import_review_decisions(args.bundle_dir, args.decisions, args.output_dir,
                                             dry_run=args.dry_run)
        elif args.review_command == "tr-plot":
            from hifivar.trgt_plot_plan import run_trgt_plot_plan

            result = run_trgt_plot_plan(args.plan, dry_run=args.dry_run)
        elif args.review_command == "snapshots":
            from hifivar.review_snapshots import run_review_snapshots

            result = run_review_snapshots(config, args.output_dir, settings_file=args.settings, dry_run=args.dry_run)
        else:
            raise ConfigurationError(f"Unsupported review action: {args.review_command}")
        displayed = ({key:result[key] for key in ('status','target_count','expected_screenshots',
                    'manual_review_status','variant_accuracy','output_directory')}
                    if args.review_command == 'snapshots' else result)
        print(json.dumps(displayed, ensure_ascii=False))
        if args.review_command == "tr-plot":
            return 0 if result["status"] in {"COMPLETE_TECHNICAL_SVG", "PLANNED"} else 1
        return 0
    if args.command == "primers" and args.primer_command == "design":
        from hifivar.primer_design_plan import run_primer_design_plan

        result = run_primer_design_plan(args.plan, args.output_dir, dry_run=args.dry_run)
        fields = ("status", "target_count", "candidate_pair_count", "specificity_status", "experimental_status",
                  "output_directory", "primer_sequences_file")
        print(json.dumps({key: result[key] for key in fields if key in result}, ensure_ascii=False))
        return 0
    if args.command == "benchmark" and args.benchmark_command == "sv":
        from hifivar.sv_benchmark import run_sv_benchmark

        result = run_sv_benchmark(args.manifest, args.output_dir, dry_run=args.dry_run)
        print(json.dumps(result, ensure_ascii=False))
        return 0 if result["status"] in {"COMPLETE", "PLANNED"} else 1
    if args.command == "benchmark" and args.benchmark_command == "sv-depth":
        from hifivar.sv_depth_study import run_sv_depth_study

        result = run_sv_depth_study(args.manifest, args.output_dir, dry_run=args.dry_run,
                                   resume=args.resume, unit_ids=args.unit_ids, callers=args.callers)
        print(json.dumps(result, ensure_ascii=False))
        return 0 if result.get("selection_status", result.get("status")) in {"COMPLETE", "PLANNED"} else 1
    if args.command == "benchmark" and args.benchmark_command == "sv-depth-support":
        from hifivar.sv_depth_support import run_sv_depth_support

        result = run_sv_depth_support(args.manifest, args.output_dir, dry_run=args.dry_run,
                                      unit_ids=args.unit_ids)
        print(json.dumps(result, ensure_ascii=False))
        return 0 if result.get("selection_status", result.get("status")) in {"COMPLETE", "PLANNED"} else 1
    if args.command == "functional" and args.functional_command == "alphagenome":
        from hifivar.functional_cli import run_alphagenome_plan

        result = run_alphagenome_plan(args.plan, args.output_dir, dry_run=args.dry_run,
                                     allow_network=args.allow_network)
        # Full score rows stay in the result file, not the console or command log.
        print(json.dumps({"status": result["status"], "selected_variant_count": result["selected_variant_count"],
                          "output_directory": str(args.output_dir.resolve()),
                          "network_execution": result["network_execution"]}, ensure_ascii=False))
        return 0 if result["status"] in {"COMPLETE", "PLANNED"} else 1
    if args.command == "report":
        if args.dry_run:
            raise ConfigurationError("report --dry-run is not supported; no report files were created. Omit --dry-run to generate a fresh local report.")
        from hifivar.workflow_report import build_launcher_report

        result = build_launcher_report(args.run_dir, args.output_dir, attempt_id=args.attempt_id)
        print(json.dumps(result, ensure_ascii=False))
        print(f"Report generated: report status={result['status']}; selected attempt execution={result['execution_status']}; scientific accuracy={result['scientific_accuracy']}.", file=sys.stderr)
        if result["status"] == "PARTIAL":
            print("PARTIAL may reflect missing per-track receipts; it does not by itself mean calling failed. Read the selected execution status separately.", file=sys.stderr)
        elif result["status"] == "FAILED":
            print("The selected attempt remains FAILED. Report generation did not rerun or repair the workflow.", file=sys.stderr)
        return 0
    raise ConfigurationError(f"Unsupported CLI command: {args.command}")


def _dispatch_config(
    args: argparse.Namespace,
    config: HiFiVarConfig,
) -> int:
    """Dispatch configuration inspection actions."""
    if args.config_command == "show":
        displayed = redact_sensitive_data(config.to_dict())
        sys.stdout.write(
            yaml.safe_dump(
                displayed,
                allow_unicode=True,
                default_flow_style=False,
                sort_keys=False,
            )
        )
        return 0
    if args.config_command == "validate":
        validate_config(config, require_complete=True)
        print("Configuration valid.")
        return 0
    if args.config_command == "dump-effective":
        write_effective_config(config, args.output)
        print(f"Effective configuration written to: {args.output}")
        return 0
    raise ConfigurationError(
        f"Unsupported config action: {args.config_command}"
    )


def _run_doctor() -> int:
    """Report only the currently implemented HiFiVar foundation environment."""
    print("HiFiVar doctor")
    print(f"Version: {__version__}")
    print(f"Python: {platform.python_version()}")
    print(f"Platform: {platform.platform()}")
    print("Configuration: OK")
    return 0


__all__ = ["build_parser", "main"]
