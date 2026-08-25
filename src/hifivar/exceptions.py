"""Base exception hierarchy for HiFiVar."""


class HiFiVarError(Exception):
    """Base class for all expected HiFiVar errors."""


class ConfigurationError(HiFiVarError):
    """Raised when HiFiVar configuration is invalid or incomplete."""


class InputValidationError(HiFiVarError):
    """Raised when an analysis input fails validation."""


class ReferenceError(HiFiVarError):
    """Raised when reference data are missing, invalid, or incompatible."""


class ToolNotFoundError(HiFiVarError):
    """Raised when a required external executable cannot be found."""


class ToolVersionError(HiFiVarError):
    """Raised when an external tool version is unsupported or unreadable."""


class CommandExecutionError(HiFiVarError):
    """Raised when an external command fails to execute successfully."""


class OutputValidationError(HiFiVarError):
    """Raised when an expected output is missing or invalid."""


class WorkflowError(HiFiVarError):
    """Raised when workflow construction or execution fails."""
