"""Tests for the HiFiVar base exception hierarchy."""

import pytest

from hifivar.exceptions import (
    CommandExecutionError,
    ConfigurationError,
    HiFiVarError,
    InputValidationError,
    OutputValidationError,
    ReferenceError,
    ToolNotFoundError,
    ToolVersionError,
    WorkflowError,
)


EXCEPTION_TYPES = (
    ConfigurationError,
    InputValidationError,
    ReferenceError,
    ToolNotFoundError,
    ToolVersionError,
    CommandExecutionError,
    OutputValidationError,
    WorkflowError,
)


@pytest.mark.parametrize("exception_type", EXCEPTION_TYPES)
def test_exceptions_inherit_from_hifivar_error(
    exception_type: type[HiFiVarError],
) -> None:
    """Every public domain exception should share the package base class."""
    assert issubclass(exception_type, HiFiVarError)


@pytest.mark.parametrize("exception_type", EXCEPTION_TYPES)
def test_exceptions_can_be_raised_and_caught_by_base_class(
    exception_type: type[HiFiVarError],
) -> None:
    """Concrete errors should be catchable as :class:`HiFiVarError`."""
    with pytest.raises(HiFiVarError, match="example failure"):
        raise exception_type("example failure")


@pytest.mark.parametrize("exception_type", (HiFiVarError, *EXCEPTION_TYPES))
def test_exception_message_is_preserved(
    exception_type: type[HiFiVarError],
) -> None:
    """Exception construction should preserve its diagnostic message."""
    error = exception_type("diagnostic message")

    assert str(error) == "diagnostic message"
