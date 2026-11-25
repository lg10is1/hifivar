from hifivar import __version__


def test_current_release_candidate_version() -> None:
    assert __version__ == "0.1.0rc5"
