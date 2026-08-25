from pathlib import Path

from hifivar.cohort import CohortDefinition, CohortTrack, CohortTrackResult, SampleCallState
from hifivar.exceptions import OutputValidationError
from hifivar.phase12 import run_phase12
from hifivar.reference import Contig, ReferenceGenome


def test_track_failure_does_not_block_other_tracks(tmp_path: Path) -> None:
    cohort = CohortDefinition("C1", ("S1",), ReferenceGenome(tmp_path / "ref.fa", tmp_path / "ref.fa.fai", "GRCh38", (Contig("chr1", 10),)))
    executed = []
    def fail():
        executed.append("small")
        raise OutputValidationError("bad joint VCF")
    def sv():
        executed.append("sv")
        return CohortTrackResult(CohortTrack.SV, True, SampleCallState.CALLED, "native", None)
    report = run_phase12(
        cohort,
        enabled_tracks={CohortTrack.SMALL_VARIANTS: True, CohortTrack.SV: True, CohortTrack.TR: False},
        operations={CohortTrack.SMALL_VARIANTS: fail, CohortTrack.SV: sv},
    )
    assert executed == ["small", "sv"]
    assert report.failed_tracks == (CohortTrack.SMALL_VARIANTS,)
    assert report.successful_tracks == (CohortTrack.SV,)
    assert report.manifest.tracks[2].state is SampleCallState.DISABLED
    payload = report.manifest.to_dict()
    serialized = str(payload).lower()
    assert "pathogenic" not in serialized
    assert "disease_causality" not in serialized
    assert "truth" not in serialized
