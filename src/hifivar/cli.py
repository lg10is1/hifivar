"""Command-line interface foundation for HiFiVar."""

from __future__ import annotations

import argparse
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
        help="retain dry-run state for future workflow commands",
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

    subparsers.add_parser(
        "doctor",
        help="report the local HiFiVar foundation environment",
    )
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
    """Dispatch only the useful Phase 0.7 commands."""
    if args.command == "config":
        return _dispatch_config(args, config)
    if args.command == "doctor":
        return _run_doctor()
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
