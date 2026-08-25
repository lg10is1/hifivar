"""Shared logging infrastructure for the HiFiVar namespace."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TextIO

from hifivar.exceptions import ConfigurationError


HIFIVAR_LOGGER_NAME = "hifivar"
DEFAULT_LOG_FORMAT = "%(asctime)s | %(levelname)s | %(name)s | %(message)s"
DEFAULT_DATE_FORMAT = "%Y-%m-%d %H:%M:%S"

_MANAGED_HANDLER_ATTRIBUTE = "_hifivar_managed_handler"
_LOG_LEVELS = {
    "DEBUG": logging.DEBUG,
    "INFO": logging.INFO,
    "WARNING": logging.WARNING,
    "ERROR": logging.ERROR,
    "CRITICAL": logging.CRITICAL,
}
SUPPORTED_LOG_LEVELS = tuple(_LOG_LEVELS)


def configure_logging(
    level: str = "INFO",
    log_file: str | Path | None = None,
    *,
    stream: TextIO | None = None,
) -> logging.Logger:
    """Configure console and optional file logging for HiFiVar.

    Only handlers installed by this function are replaced on subsequent calls.
    The Python root logger and third-party logger configurations are left intact.

    Args:
        level: Minimum HiFiVar log level. Matching is case-insensitive.
        log_file: Optional UTF-8 log file. Missing parent directories are created.
        stream: Optional console stream, primarily useful for embedding and tests.

    Returns:
        The configured ``hifivar`` namespace logger.

    Raises:
        ConfigurationError: If ``level`` is not a supported logging level.
        OSError: If the log directory or file cannot be created.
    """
    numeric_level = parse_log_level(level)
    formatter = logging.Formatter(DEFAULT_LOG_FORMAT, datefmt=DEFAULT_DATE_FORMAT)

    new_handlers: list[logging.Handler] = [logging.StreamHandler(stream)]
    try:
        if log_file is not None:
            log_path = Path(log_file)
            log_path.parent.mkdir(parents=True, exist_ok=True)
            new_handlers.append(
                logging.FileHandler(log_path, encoding="utf-8")
            )
    except OSError:
        for handler in new_handlers:
            handler.close()
        raise

    for handler in new_handlers:
        handler.setFormatter(formatter)
        setattr(handler, _MANAGED_HANDLER_ATTRIBUTE, True)

    namespace_logger = logging.getLogger(HIFIVAR_LOGGER_NAME)
    _remove_managed_handlers(namespace_logger)
    namespace_logger.setLevel(numeric_level)
    namespace_logger.propagate = False
    namespace_logger.disabled = False

    for handler in new_handlers:
        namespace_logger.addHandler(handler)

    return namespace_logger


def get_logger(name: str | None = None) -> logging.Logger:
    """Return a logger within the ``hifivar`` namespace.

    Fully qualified HiFiVar names are preserved. Other module-like names are
    prefixed so their records use the shared HiFiVar configuration.
    """
    if not name or name == HIFIVAR_LOGGER_NAME:
        logger_name = HIFIVAR_LOGGER_NAME
    elif name.startswith(f"{HIFIVAR_LOGGER_NAME}."):
        logger_name = name
    else:
        logger_name = f"{HIFIVAR_LOGGER_NAME}.{name}"

    return logging.getLogger(logger_name)


def parse_log_level(level: str) -> int:
    """Translate a supported level name into its standard-library value."""
    if not isinstance(level, str):
        raise ConfigurationError(
            "Logging level must be a string: "
            + ", ".join(_LOG_LEVELS)
        )

    normalized_level = level.strip().upper()
    try:
        return _LOG_LEVELS[normalized_level]
    except KeyError as error:
        supported = ", ".join(_LOG_LEVELS)
        raise ConfigurationError(
            f"Invalid logging level {level!r}. Expected one of: {supported}."
        ) from error


def _remove_managed_handlers(logger: logging.Logger) -> None:
    """Detach and close handlers previously installed by HiFiVar."""
    for handler in list(logger.handlers):
        if getattr(handler, _MANAGED_HANDLER_ATTRIBUTE, False):
            logger.removeHandler(handler)
            handler.close()


__all__ = [
    "DEFAULT_DATE_FORMAT",
    "DEFAULT_LOG_FORMAT",
    "HIFIVAR_LOGGER_NAME",
    "SUPPORTED_LOG_LEVELS",
    "configure_logging",
    "get_logger",
    "parse_log_level",
]
