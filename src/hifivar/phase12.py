"""Phase 12 orchestration that isolates cohort-track failures."""

from __future__ import annotations

from collections.abc import Callable, Mapping
from dataclasses import dataclass

from hifivar import __version__
from hifivar.cohort import (
    CohortDefinition,
    CohortManifest,
    CohortTrack,
    CohortTrackResult,
    SampleCallState,
)
from hifivar.exceptions import HiFiVarError, InputValidationError


TrackOperation = Callable[[], CohortTrackResult]


@dataclass(frozen=True, slots=True)
class Phase12Report:
    manifest: CohortManifest

    @property
    def failed_tracks(self) -> tuple[CohortTrack, ...]:
        return tuple(item.track for item in self.manifest.tracks if item.state is SampleCallState.FAILED)

    @property
    def successful_tracks(self) -> tuple[CohortTrack, ...]:
        return tuple(item.track for item in self.manifest.tracks if item.state in {SampleCallState.CALLED, SampleCallState.NO_CALLS})


def run_phase12(
    cohort: CohortDefinition,
    *,
    enabled_tracks: Mapping[CohortTrack, bool],
    operations: Mapping[CohortTrack, TrackOperation],
) -> Phase12Report:
    """Run tracks independently; one expected failure cannot erase another result."""
    results: list[CohortTrackResult] = []
    for track in CohortTrack:
        enabled = enabled_tracks.get(track, False)
        if not isinstance(enabled, bool):
            raise InputValidationError(f"Enabled state for {track.value} must be boolean.")
        if not enabled:
            results.append(CohortTrackResult(track, False, SampleCallState.DISABLED, None, None))
            continue
        operation = operations.get(track)
        if operation is None:
            results.append(CohortTrackResult(track, True, SampleCallState.NOT_RUN, None, None, message="No track operation was supplied."))
            continue
        try:
            result = operation()
        except HiFiVarError as error:
            results.append(CohortTrackResult(track, True, SampleCallState.FAILED, None, None, message=str(error)))
            continue
        if result.track is not track:
            raise InputValidationError(f"Phase 12 operation for {track.value} returned {result.track.value}.")
        results.append(result)
    return Phase12Report(CohortManifest(cohort, tuple(results), __version__))


__all__ = ["Phase12Report", "TrackOperation", "run_phase12"]
