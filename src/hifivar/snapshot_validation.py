"""Bounded PNG integrity checks, explicitly separate from biological review."""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
from pathlib import Path
import struct
from typing import Sequence
import warnings
import zlib

from hifivar.exceptions import ConfigurationError, OutputValidationError

PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"
VALIDATION_SCOPE = "PNG_CRC_AND_FULL_DECODE_NOT_TRACK_OR_HUMAN_QA"


@dataclass(frozen=True, slots=True)
class SnapshotValidationPolicy:
    """Technical limits; passing them does not establish image interpretability."""

    minimum_width: int = 64
    minimum_height: int = 64
    maximum_pixels: int = 50_000_000
    maximum_bytes: int = 100 * 1024 * 1024

    def __post_init__(self) -> None:
        for name in ("minimum_width", "minimum_height", "maximum_pixels", "maximum_bytes"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
                raise ConfigurationError(f"Snapshot {name} must be a positive integer.")
        if self.minimum_width * self.minimum_height > self.maximum_pixels:
            raise ConfigurationError("Snapshot minimum dimensions exceed maximum_pixels.")


def require_png_decoder():
    """Resolve the optional review dependency before running the GUI."""
    try:
        from PIL import Image, ImageFile
    except ImportError as error:
        raise ConfigurationError(
            "IGV screenshot validation requires Pillow; install HiFiVar with the [review] extra."
        ) from error
    if ImageFile.LOAD_TRUNCATED_IMAGES:
        raise ConfigurationError("IGV snapshot validation requires Pillow truncated-image loading to be disabled.")
    return Image


def _chunk_integrity(path: Path, policy: SnapshotValidationPolicy) -> tuple[int, int]:
    """Validate every chunk CRC, required structure, IEND and absence of trailers."""
    size = path.stat().st_size
    if not 0 < size <= policy.maximum_bytes:
        raise OutputValidationError("IGV snapshot is empty or exceeds maximum_bytes.")
    with path.open("rb") as stream:
        if stream.read(8) != PNG_SIGNATURE:
            raise OutputValidationError("IGV snapshot does not have a PNG signature.")
        width = height = 0
        saw_idat = False
        index = 0
        while True:
            header = stream.read(8)
            if len(header) != 8:
                raise OutputValidationError("IGV snapshot has a truncated PNG chunk or missing IEND.")
            length, kind = struct.unpack(">I4s", header)
            if length > policy.maximum_bytes or length + 4 > size - stream.tell():
                raise OutputValidationError("IGV snapshot has an invalid or truncated PNG chunk length.")
            if index == 0 and (kind != b"IHDR" or length != 13):
                raise OutputValidationError("IGV snapshot must start with a 13-byte IHDR.")
            if index > 0 and kind == b"IHDR":
                raise OutputValidationError("IGV snapshot contains a duplicate IHDR.")
            crc = zlib.crc32(kind)
            payload_head = b""
            remaining = length
            while remaining:
                block = stream.read(min(remaining, 1024 * 1024))
                if not block:
                    raise OutputValidationError("IGV snapshot has truncated PNG data.")
                if kind == b"IHDR":
                    payload_head += block
                crc = zlib.crc32(block, crc)
                remaining -= len(block)
            stored_crc = stream.read(4)
            if len(stored_crc) != 4 or struct.unpack(">I", stored_crc)[0] != crc & 0xFFFFFFFF:
                raise OutputValidationError("IGV snapshot PNG chunk CRC mismatch.")
            if kind == b"IHDR":
                width, height = struct.unpack(">II", payload_head[:8])
                if (width < policy.minimum_width or height < policy.minimum_height
                        or width * height > policy.maximum_pixels):
                    raise OutputValidationError("IGV snapshot dimensions are outside the technical policy.")
            elif kind == b"IDAT":
                saw_idat = True
            elif kind == b"IEND":
                if length or not saw_idat or stream.read(1):
                    raise OutputValidationError("IGV snapshot has invalid IEND, missing IDAT or trailing data.")
                break
            index += 1
    return width, height


def validate_png_snapshot(
    path: str | Path, *, policy: SnapshotValidationPolicy | None = None,
) -> dict[str, object]:
    """Verify structural integrity and decode every pixel; do not score evidence."""
    policy = SnapshotValidationPolicy() if policy is None else policy
    if not isinstance(policy, SnapshotValidationPolicy):
        raise ConfigurationError("Snapshot validation requires SnapshotValidationPolicy.")
    if not isinstance(path, (str, Path)) or not str(path).strip():
        raise OutputValidationError("IGV snapshot path must be nonempty.")
    selected = Path(path)
    Image = require_png_decoder()
    try:
        before = selected.stat()
        width, height = _chunk_integrity(selected, policy)
        with warnings.catch_warnings():
            warnings.simplefilter("error", Image.DecompressionBombWarning)
            with Image.open(selected) as image:
                if image.format != "PNG" or image.size != (width, height):
                    raise OutputValidationError("IGV snapshot decoder/header dimensions differ.")
                if getattr(image, 'is_animated', False):
                    raise OutputValidationError("IGV screenshot validation requires a static PNG, not APNG.")
                image.verify()
            with Image.open(selected) as image:
                image.load()
                mode = image.mode
        digest = hashlib.sha256()
        with selected.open("rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(block)
        after = selected.stat()
        if (before.st_size, before.st_mtime_ns, before.st_ino) != (after.st_size, after.st_mtime_ns, after.st_ino):
            raise OutputValidationError("IGV snapshot changed during validation.")
    except OutputValidationError:
        raise
    except (OSError, ValueError, SyntaxError, struct.error, Image.DecompressionBombWarning,
            Image.DecompressionBombError) as error:
        raise OutputValidationError(f"IGV snapshot cannot be fully decoded: '{selected}': {error}") from error
    return {
        "path": str(selected), "width": width, "height": height,
        "bytes": after.st_size, "sha256": digest.hexdigest(), "mode": mode,
        "validation": VALIDATION_SCOPE, "human_review": "NOT_REVIEWED",
        "track_loading": "NOT_ASSESSED", "biological_accuracy": "NOT_ASSESSED",
    }


def validate_png_snapshots(
    paths: Sequence[Path], *, policy: SnapshotValidationPolicy | None = None,
) -> tuple[dict[str, object], ...]:
    """Validate each expected screenshot once, rejecting duplicate destinations."""
    ordered = tuple(Path(path) for path in paths)
    if len({path.resolve() for path in ordered}) != len(ordered):
        raise OutputValidationError("IGV snapshot list contains duplicate output paths.")
    return tuple(validate_png_snapshot(path, policy=policy) for path in ordered)


__all__ = ["SnapshotValidationPolicy", "VALIDATION_SCOPE", "require_png_decoder",
           "validate_png_snapshot", "validate_png_snapshots"]
