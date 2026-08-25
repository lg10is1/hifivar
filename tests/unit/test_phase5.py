import pytest

from hifivar.exceptions import ConfigurationError
from hifivar.phase5 import resolve_trgt_karyotype
from hifivar.sample_sheet import Sex


def test_auto_karyotype_uses_only_declared_sample_sex() -> None:
    assert resolve_trgt_karyotype("auto", Sex.FEMALE, "F1") == "XX"
    assert resolve_trgt_karyotype("auto", Sex.MALE, "M1") == "XY"
    assert resolve_trgt_karyotype("XX", None, "S1") == "XX"


@pytest.mark.parametrize("sex", (None, Sex.UNKNOWN))
def test_auto_karyotype_refuses_missing_or_unknown_metadata(sex) -> None:
    with pytest.raises(ConfigurationError, match="cannot be inferred"):
        resolve_trgt_karyotype("auto", sex, "S1")
