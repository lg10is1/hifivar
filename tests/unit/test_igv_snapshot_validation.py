"""Image integrity is technical QA, never variant correctness or native IGV QA."""
import builtins
import hashlib
from pathlib import Path
import struct
import zlib

import pytest

from hifivar.exceptions import ConfigurationError, InputValidationError, OutputValidationError
from hifivar.igv import IgvWrapper
from hifivar.review import VariantClass
from hifivar.snapshot_validation import (
    SnapshotValidationPolicy, VALIDATION_SCOPE, validate_png_snapshot, validate_png_snapshots,
)
from test_review import FakeIgvRunner, PNG, png_fixture, target


def save(tmp_path, data=PNG):
    path = tmp_path / 'snapshot.png'
    path.write_bytes(data)
    return path


def chunk(kind, payload):
    return struct.pack('>I', len(payload)) + kind + payload + struct.pack('>I', zlib.crc32(kind + payload) & 0xffffffff)


def test_full_decode_metadata_and_source_immutability(tmp_path):
    path = save(tmp_path)
    result = validate_png_snapshot(path)
    assert (result['width'], result['height']) == (64, 64)
    assert result['sha256'] == hashlib.sha256(PNG).hexdigest()
    assert result['validation'] == VALIDATION_SCOPE
    assert result['human_review'] == 'NOT_REVIEWED'
    assert result['track_loading'] == result['biological_accuracy'] == 'NOT_ASSESSED'
    assert path.read_bytes() == PNG


@pytest.mark.parametrize('data', [b'', b'PNG', PNG[:24], PNG[:-1], PNG[:-12], PNG + b'trailer'])
def test_truncated_or_non_png_or_trailing_bytes_are_rejected(tmp_path, data):
    with pytest.raises(OutputValidationError):
        validate_png_snapshot(save(tmp_path, data))


def test_bad_crc_is_rejected_even_when_pixels_decode(tmp_path):
    payload = bytearray(PNG)
    payload[29] ^= 1  # IHDR CRC only, unchanged pixel stream.
    with pytest.raises(OutputValidationError, match='CRC mismatch'):
        validate_png_snapshot(save(tmp_path, payload))


def test_valid_crc_but_invalid_compressed_pixel_stream_is_rejected(tmp_path):
    ihdr = chunk(b'IHDR', struct.pack('>IIBBBBB', 64, 64, 8, 2, 0, 0, 0))
    payload = PNG[:8] + ihdr + chunk(b'IDAT', b'not-zlib') + chunk(b'IEND', b'')
    with pytest.raises(OutputValidationError, match='fully decoded'):
        validate_png_snapshot(save(tmp_path, payload))


def test_decode_rejects_incomplete_pixel_rows_with_structurally_valid_chunks(tmp_path):
    ihdr = chunk(b'IHDR', struct.pack('>IIBBBBB', 64, 64, 8, 2, 0, 0, 0))
    payload = PNG[:8] + ihdr + chunk(b'IDAT', zlib.compress(b'\x00\xff\xff\xff')) + chunk(b'IEND', b'')
    with pytest.raises(OutputValidationError):
        validate_png_snapshot(save(tmp_path, payload))


@pytest.mark.parametrize('width,height', [(1, 1), (64, 1), (1, 64)])
def test_too_small_images_are_not_accepted_as_screenshots(tmp_path, width, height):
    with pytest.raises(OutputValidationError, match='dimensions'):
        validate_png_snapshot(save(tmp_path, png_fixture(width, height)))


def test_byte_and_pixel_limits_are_enforced_before_decode(tmp_path):
    path = save(tmp_path)
    with pytest.raises(OutputValidationError, match='maximum_bytes'):
        validate_png_snapshot(path, policy=SnapshotValidationPolicy(maximum_bytes=10))
    with pytest.raises(OutputValidationError, match='dimensions'):
        validate_png_snapshot(path, policy=SnapshotValidationPolicy(minimum_width=1, minimum_height=1, maximum_pixels=100))


@pytest.mark.parametrize('value', [0, -1, True, 1.5, '64'])
def test_invalid_policy_is_typed(value):
    with pytest.raises(ConfigurationError):
        SnapshotValidationPolicy(minimum_width=value)


@pytest.mark.parametrize('policy', [False,0,'invalid'])
def test_invalid_policy_objects_are_not_silently_replaced(tmp_path, policy):
    with pytest.raises(ConfigurationError):
        validate_png_snapshot(save(tmp_path), policy=policy)


def test_global_truncated_loading_is_rejected_without_mutation(tmp_path, monkeypatch):
    from PIL import ImageFile
    monkeypatch.setattr(ImageFile, 'LOAD_TRUNCATED_IMAGES', True)
    with pytest.raises(ConfigurationError, match='truncated-image'):
        validate_png_snapshot(save(tmp_path))
    assert ImageFile.LOAD_TRUNCATED_IMAGES is True


def test_duplicate_output_destinations_and_missing_files(tmp_path):
    path = save(tmp_path)
    with pytest.raises(OutputValidationError, match='duplicate'):
        validate_png_snapshots((path, path.parent / '.' / path.name))
    with pytest.raises(OutputValidationError, match='fully decoded'):
        validate_png_snapshot(tmp_path / 'missing.png')


def test_all_screenshots_and_count_recorded_in_wrapper_result(tmp_path):
    item = target(tmp_path, VariantClass.BND, end=100)
    wrapper = IgvWrapper(runner=FakeIgvRunner())
    plan = wrapper.plan((item,), batch_path=item.output_directory/'review.igv.batch', snapshot_directory=item.output_directory/'screenshots')
    result = wrapper.run(plan)
    record = result.to_dict()['snapshot_validation']
    assert record['expected_count'] == record['validated_count'] == 2
    assert all(image['validation'] == VALIDATION_SCOPE for image in record['images'])
    assert record['human_review'] == 'NOT_REVIEWED'


def test_dependency_failure_precedes_native_version_or_gui_execution(tmp_path, monkeypatch):
    item = target(tmp_path, VariantClass.SNV)
    runner = FakeIgvRunner()
    wrapper = IgvWrapper(runner=runner)
    plan = wrapper.plan((item,), batch_path=item.output_directory/'review.igv.batch', snapshot_directory=item.output_directory/'screenshots')
    original = builtins.__import__
    def unavailable(name, *args, **kwargs):
        if name == 'PIL' or name.startswith('PIL.'):
            raise ImportError('synthetic optional dependency absent')
        return original(name, *args, **kwargs)
    monkeypatch.setattr(builtins, '__import__', unavailable)
    with pytest.raises(ConfigurationError, match='Pillow'):
        wrapper.run(plan)
    assert not runner.calls and not plan.batch_path.exists()
    planned = wrapper.run(plan, dry_run=True)
    assert planned.status.value == 'planned' and not plan.batch_path.exists()


def test_corrupt_generated_snapshot_never_becomes_completed_evidence(tmp_path):
    item = target(tmp_path, VariantClass.SNV)
    class CorruptRunner(FakeIgvRunner):
        def run(self, argv, **kwargs):
            result = super().run(argv, **kwargs)
            if '--batch' in argv:
                item.screenshot_paths[0].write_bytes(PNG[:24])
            return result
    wrapper = IgvWrapper(runner=CorruptRunner())
    plan = wrapper.plan((item,), batch_path=item.output_directory/'review.igv.batch', snapshot_directory=item.output_directory/'screenshots')
    with pytest.raises(OutputValidationError):
        wrapper.run(plan)
    assert item.screenshot_paths[0].read_bytes() == PNG[:24]  # Failure evidence retained.


def test_display_settings_are_recorded_and_native_calls_are_bounded(tmp_path):
    item = target(tmp_path, VariantClass.SNV)
    runner = FakeIgvRunner()
    wrapper = IgvWrapper(runner=runner)
    plan = wrapper.plan((item,), batch_path=item.output_directory/'review.igv.batch',
                        snapshot_directory=item.output_directory/'screenshots',
                        max_panel_height=1200, sleep_interval_ms=500)
    assert 'maxPanelHeight 1200' in plan.batch_text and 'setSleepInterval 500' in plan.batch_text
    # IGV BatchRunner decides default-genome loading from this first command.
    assert plan.batch_text.splitlines()[0].startswith('genome ')
    assert '"' not in plan.batch_text.splitlines()[0]
    assert plan.to_dict()['display_settings'] == {'max_panel_height':1200, 'sleep_interval_ms':500}
    wrapper.run(plan, timeout_seconds=60)
    assert runner.calls[0][1]['timeout'] == 30
    assert runner.calls[-1][1]['timeout'] == 60


@pytest.mark.parametrize('settings', [dict(max_panel_height=True),dict(max_panel_height=1),
                                    dict(sleep_interval_ms=-1),dict(sleep_interval_ms=1.5)])
def test_invalid_display_settings_are_rejected(tmp_path, settings):
    item = target(tmp_path, VariantClass.SNV)
    with pytest.raises(InputValidationError, match='must be an integer'):
        IgvWrapper().plan((item,), batch_path=item.output_directory/'review.igv.batch',
                          snapshot_directory=item.output_directory/'screenshots', **settings)


@pytest.mark.parametrize('timeout', [True,0,-1,float('inf'),float('nan'),'60'])
def test_invalid_native_timeout_is_rejected_before_execution(tmp_path, timeout):
    item = target(tmp_path, VariantClass.SNV)
    runner = FakeIgvRunner()
    wrapper = IgvWrapper(runner=runner)
    plan = wrapper.plan((item,), batch_path=item.output_directory/'review.igv.batch',snapshot_directory=item.output_directory/'screenshots')
    with pytest.raises(InputValidationError, match='timeout_seconds'):
        wrapper.run(plan, timeout_seconds=timeout)
    assert not runner.calls
