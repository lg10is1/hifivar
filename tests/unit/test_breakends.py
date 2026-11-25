"""Normative grammar examples and fail-closed BND policy tests."""
import hashlib
import json

import pytest

from hifivar.breakends import BndPolicy, POLICY, POLICY_SHA256, parse_breakend
from hifivar.exceptions import ConfigurationError, InputValidationError


LENGTHS = {"chr1": 1000, "chr2": 2000, "aux:segment": 100}


@pytest.mark.parametrize("alt,local,mate", [
    ("A[chr2:12[", "LEFT", "RIGHT"),
    ("A]chr2:12]", "LEFT", "LEFT"),
    ("]chr2:12]A", "RIGHT", "LEFT"),
    ("[chr2:12[A", "RIGHT", "RIGHT"),
])
@pytest.mark.parametrize("position", [10, 12, 15])
def test_four_orientations_do_not_imply_an_event_span(alt, local, mate, position):
    junction = parse_breakend("chr1", position, "A", alt, LENGTHS, svtype="BND")
    assert junction is not None
    assert (junction.local_side, junction.mate_side) == (local, mate)
    assert (junction.mate_contig, junction.mate_position) == ("chr2", 12)
    assert junction.to_dict()["coordinate_semantics"] == "VCF_ANCHORS_NOT_EVENT_SPAN"


@pytest.mark.parametrize("alt,side", [("AC.", "LEFT"), (".AC", "RIGHT")])
def test_single_breakend_does_not_invent_a_mate(alt, side):
    junction = parse_breakend("chr1", 20, "AC", alt, LENGTHS)
    assert junction.local_side == side
    assert junction.mate_contig is junction.mate_position is None


@pytest.mark.parametrize("position", [0, 1001])
@pytest.mark.parametrize("alt", [".[chr2:12[", ".]chr2:12]", "]chr2:12].", "[chr2:12[."])
def test_virtual_local_anchor_is_retained_without_clamping(position, alt):
    junction = parse_breakend("chr1", position, "N", alt, LENGTHS)
    assert junction.position == position
    assert junction.local_virtual


def test_exact_official_telomere_example():
    junction = parse_breakend("1", 0, "N", ".[13:123457[", {"1": 1000, "13": 200000})
    assert junction.local_virtual and junction.mate_position == 123457


@pytest.mark.parametrize("mate", [0, 2001])
def test_virtual_mate_and_colon_in_contig_name(mate):
    junction = parse_breakend("chr1", 1, "A", f"A[chr2:{mate}[", LENGTHS)
    assert junction.mate_position == mate and junction.mate_virtual
    assert parse_breakend("chr1", 1, "A", "A[aux:segment:10[", LENGTHS).mate_contig == "aux:segment"


@pytest.mark.parametrize("alt", [
    "A[chr2:10]", "A[chr2:10[A", "[chr2:10[", "A[chr2:10[A[chr2:11[",
    "A[chr2:-1[", "A[chr2:2002[", "A[chrX:12[", "A[chr2:1.2[",
    "A..", "..A", "A.C.", "A[chr2:10[.", ".[chr2:10[", "A,C",
])
def test_malformed_or_unsupported_alt_fails(alt):
    with pytest.raises(InputValidationError):
        parse_breakend("chr1", 2, "A", alt, LENGTHS, svtype="BND")


@pytest.mark.parametrize("position,ref,alt", [
    (-1, "N", ".[chr2:10["), (1002, "N", ".[chr2:10["),
    (0, "A", ".[chr2:10["), (0, "N", "A[chr2:10["), (0, "N", "N."),
    (1000, "AC", "AC[chr2:10["),
])
def test_bad_local_anchor_fails(position, ref, alt):
    with pytest.raises(InputValidationError):
        parse_breakend("chr1", position, ref, alt, LENGTHS)


def test_non_bnd_and_svtype_disagreement():
    assert parse_breakend("chr1", 0, "N", "<INS>", LENGTHS, svtype="INS") is None
    with pytest.raises(InputValidationError):
        parse_breakend("chr1", 2, "A", "A[chr2:10[", LENGTHS, svtype="DEL")
    with pytest.raises(InputValidationError):
        parse_breakend("chr1", 2, "A", "<BND>", LENGTHS, svtype="BND")


def test_policy_identity_is_semantic_json_not_line_endings():
    expected = hashlib.sha256(json.dumps(POLICY, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("ascii")).hexdigest()
    assert expected == POLICY_SHA256
    policy = BndPolicy(POLICY["policy_id"], 1, POLICY_SHA256, POLICY["mode"])
    assert BndPolicy.from_config(policy.to_dict()) == policy
    assert BndPolicy.from_config({"enabled": False}) is None


@pytest.mark.parametrize("config", [
    {"enabled": 1}, {"enabled": True}, {"enabled": False, "policy_id": "active"},
    {"enabled": False, "timeout_seconds": False}, {"typo": True},
])
def test_bad_config_fails(config):
    with pytest.raises(ConfigurationError):
        BndPolicy.from_config(config)


@pytest.mark.parametrize("key,value", [("policy_sha256", "0" * 64), ("policy_version", True), ("mode", "rewrite")])
def test_policy_drift_is_rejected(key, value):
    config = {"enabled": True, "policy_id": POLICY["policy_id"], "policy_version": 1,
              "policy_sha256": POLICY_SHA256, "mode": POLICY["mode"]}
    config[key] = value
    with pytest.raises(ConfigurationError):
        BndPolicy.from_config(config)
