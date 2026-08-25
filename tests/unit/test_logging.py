"""Tests for the shared HiFiVar logging infrastructure."""

from __future__ import annotations

import io
import logging
import re
from collections.abc import Iterator
from pathlib import Path

import pytest

from hifivar.exceptions import ConfigurationError
from hifivar.logging_utils import (
    HIFIVAR_LOGGER_NAME,
    configure_logging,
    get_logger,
)


@pytest.fixture(autouse=True)
def isolate_hifivar_logging() -> Iterator[None]:
    """Restore the HiFiVar namespace logger after every test."""
    namespace_logger = logging.getLogger(HIFIVAR_LOGGER_NAME)
    original_handlers = list(namespace_logger.handlers)
    original_level = namespace_logger.level
    original_propagate = namespace_logger.propagate
    original_disabled = namespace_logger.disabled

    for handler in original_handlers:
        namespace_logger.removeHandler(handler)
    namespace_logger.setLevel(logging.NOTSET)
    namespace_logger.propagate = True
    namespace_logger.disabled = False

    yield

    for handler in list(namespace_logger.handlers):
        namespace_logger.removeHandler(handler)
        handler.close()
    for handler in original_handlers:
        namespace_logger.addHandler(handler)
    namespace_logger.setLevel(original_level)
    namespace_logger.propagate = original_propagate
    namespace_logger.disabled = original_disabled


def test_get_logger_returns_namespaced_logger() -> None:
    """Module loggers should be standard loggers under ``hifivar``."""
    logger = get_logger("core")

    assert isinstance(logger, logging.Logger)
    assert logger.name == "hifivar.core"


def test_info_level_filters_debug_and_emits_info(
    capsys: pytest.CaptureFixture[str],
) -> None:
    """INFO configuration should suppress DEBUG console records."""
    configure_logging(level="info")
    logger = get_logger("hifivar.info_test")

    logger.debug("hidden debug message")
    logger.info("visible info message")

    console_output = capsys.readouterr().err
    assert "hidden debug message" not in console_output
    assert "visible info message" in console_output


def test_debug_level_emits_debug(capsys: pytest.CaptureFixture[str]) -> None:
    """DEBUG configuration should emit DEBUG console records."""
    configure_logging(level="DeBuG")

    get_logger("debug_test").debug("visible debug message")

    assert "visible debug message" in capsys.readouterr().err


@pytest.mark.parametrize(
    ("level_name", "expected_level"),
    (
        ("DEBUG", logging.DEBUG),
        ("INFO", logging.INFO),
        ("WARNING", logging.WARNING),
        ("ERROR", logging.ERROR),
        ("CRITICAL", logging.CRITICAL),
    ),
)
def test_supported_levels_are_configured(
    level_name: str,
    expected_level: int,
) -> None:
    """Every documented level should map to its standard logging value."""
    namespace_logger = configure_logging(level=level_name, stream=io.StringIO())

    assert namespace_logger.level == expected_level


def test_file_logging_creates_log_file(tmp_path: Path) -> None:
    """Supplying a file path should create the log file."""
    log_file = tmp_path / "hifivar.log"
    configure_logging(log_file=log_file, stream=io.StringIO())

    get_logger("file_creation").info("create the file")

    assert log_file.is_file()


def test_file_logging_contains_message(tmp_path: Path) -> None:
    """File records should preserve the emitted message."""
    log_file = tmp_path / "hifivar.log"
    configure_logging(log_file=log_file, stream=io.StringIO())

    get_logger("file_content").info("analysis started")

    assert "analysis started" in log_file.read_text(encoding="utf-8")


def test_file_logging_supports_utf8(tmp_path: Path) -> None:
    """Unicode sample information should round-trip through a log file."""
    log_file = tmp_path / "hifivar.log"
    message = "样本 HG002 分析开始"
    configure_logging(log_file=log_file, stream=io.StringIO())

    get_logger("unicode").info(message)

    assert message in log_file.read_text(encoding="utf-8")


def test_repeated_configuration_does_not_duplicate_messages() -> None:
    """Reconfiguration should replace, rather than accumulate, handlers."""
    stream = io.StringIO()
    configure_logging(stream=stream)
    configure_logging(stream=stream)

    get_logger("duplicate_check").info("single record")

    assert stream.getvalue().count("single record") == 1
    assert len(logging.getLogger(HIFIVAR_LOGGER_NAME).handlers) == 1


def test_reconfiguration_preserves_unmanaged_namespace_handlers() -> None:
    """HiFiVar should replace only handlers that it installed itself."""
    namespace_logger = logging.getLogger(HIFIVAR_LOGGER_NAME)
    custom_handler = logging.NullHandler()
    namespace_logger.addHandler(custom_handler)

    configure_logging(stream=io.StringIO())
    configure_logging(stream=io.StringIO())

    assert custom_handler in namespace_logger.handlers
    assert len(namespace_logger.handlers) == 2


def test_configuration_does_not_modify_root_logger() -> None:
    """Configuring HiFiVar must leave global root logging state intact."""
    root_logger = logging.getLogger()
    original_handlers = list(root_logger.handlers)
    original_level = root_logger.level

    configure_logging(stream=io.StringIO())

    assert root_logger.handlers == original_handlers
    assert root_logger.level == original_level


def test_invalid_logging_level_raises_configuration_error() -> None:
    """Unsupported levels should fail instead of silently using INFO."""
    with pytest.raises(ConfigurationError, match="INVALID_LEVEL"):
        configure_logging(level="INVALID_LEVEL")


def test_file_logging_creates_missing_parent_directories(tmp_path: Path) -> None:
    """Only the requested log file's missing parents should be created."""
    log_file = tmp_path / "nested" / "logs" / "hifivar.log"
    configure_logging(log_file=log_file, stream=io.StringIO())

    get_logger("nested_path").warning("nested path ready")

    assert log_file.is_file()


def test_default_format_contains_timestamp_level_name_and_message() -> None:
    """The default formatter should expose the required record fields."""
    stream = io.StringIO()
    configure_logging(level="WARNING", stream=stream)

    get_logger("format").warning("formatted record")

    assert re.search(
        r"\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}"
        r" \| WARNING \| hifivar\.format \| formatted record",
        stream.getvalue(),
    )
