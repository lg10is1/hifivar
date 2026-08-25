import pytest
from hifivar.config import validate_config
from hifivar.exceptions import ConfigurationError


def test_benchmark_disabled_defaults_and_strict_nested_schema() -> None:
    validate_config({"benchmark":{"enabled":False,"small_variants":{"enabled":False,"threads":2,"stratifications":[]}}})
    with pytest.raises(ConfigurationError,match="Unknown"):
        validate_config({"benchmark":{"enabled":False,"sv":{"enabled":False,"mystery":1}}})


def test_enabled_benchmark_requires_explicit_truth_versions() -> None:
    with pytest.raises(ConfigurationError,match="incomplete"):
        validate_config({"benchmark":{"enabled":True,"benchmark_id":"B","sample_id":"S","small_variants":{"enabled":True,"query_vcf":"q","truth_vcf":"t"}}})


def test_happy_metrics_compression_contract_is_strict() -> None:
    validate_config({"benchmark":{"enabled":False,"small_variants":{"enabled":False,"happy_version":"0.3.15","metrics_compression":"gzip"}}})
    validate_config({"benchmark":{"enabled":False,"small_variants":{"enabled":False,"metrics_compression":"plain"}}})
    with pytest.raises(ConfigurationError,match="metrics_compression"):
        validate_config({"benchmark":{"enabled":False,"small_variants":{"enabled":False,"metrics_compression":"auto"}}})
