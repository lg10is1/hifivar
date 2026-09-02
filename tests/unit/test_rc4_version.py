from hifivar import __version__


def test_release_candidate_version_is_rc4() -> None:
    assert __version__ == "0.1.0rc4"
