"""Fake-tool tests; native BGZF/TBI behavior requires separate Linux acceptance."""
from dataclasses import replace
import gzip
import hashlib
import json
from pathlib import Path
import struct
import zlib

import pytest

from hifivar.assembly_sv import SVEvidenceSource
from hifivar.bnd_representation import BndRepresentationRequest, BndRepresentationWrapper
from hifivar.breakends import BndPolicy, POLICY, POLICY_SHA256
from hifivar.command import CommandResult
from hifivar.exceptions import CommandExecutionError, InputValidationError, OutputValidationError, ToolVersionError
from hifivar.harmonization import EvidenceRunStatus, HarmonizedSvArtifact, SVEvidenceSourceArtifact, SVHarmonizationRequest
from hifivar.jasmine import JasmineCommandPlan, JasmineResult, JasmineResultStatus
from hifivar.reference import ReferenceGenome


TEXT = ("##fileformat=VCFv4.4\n##contig=<ID=chr1,length=1000>\n"
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End">\n'
        '##INFO=<ID=SVLEN,Number=A,Type=Integer,Description="Length">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "chr1\t20\tb1\tA\tA[chr1:10[\t.\tPASS\tSVTYPE=BND;END=10;SVLEN=0\tGT\t0/1\n")


def bgzf(data):
    compressor = zlib.compressobj(wbits=-15)
    payload = compressor.compress(data) + compressor.flush()
    size = 18 + len(payload) + 8
    assert size <= 65536
    return (b"\x1f\x8b\x08\x04" + b"\x00" * 4 + b"\x00\xff\x06\x00BC\x02\x00"
            + struct.pack("<H", size - 1) + payload + struct.pack("<II", zlib.crc32(data), len(data)))


class FakeRunner:
    def __init__(self, directory, fail=None, version="1.24"):
        self.directory, self.fail, self.version, self.calls = directory, fail, version, []
        self.executables = {}
        for name in ("bgzip", "tabix"):
            path = directory / name
            path.write_bytes(b"FAKE UNIT TEST TOOL " + name.encode())
            self.executables[name] = path

    def require_executable(self, executable, **kwargs):
        return self.executables[Path(executable).name]

    def run(self, args, **kwargs):
        self.calls.append((args, kwargs))
        if self.fail == Path(args[0]).name and "--version" not in args:
            raise CommandExecutionError("deliberate fake tool failure")
        stdout = b""
        if "--version" in args:
            stdout = (f"{Path(args[0]).name} (htslib) {self.version}\n").encode()
        elif Path(args[0]).name == "bgzip":
            stdout = bgzf(Path(args[-1]).read_bytes())
        elif "-p" in args:
            Path(args[-1] + ".tbi").write_bytes(b"FAKE INDEX NOT NATIVE VALIDATION")
        elif "-l" in args:
            stdout = b"chr1\n"
        else:
            with gzip.open(args[1], "rb") as reader:
                stdout = b"".join(line for line in reader if not line.startswith(b"#"))
        stderr = kwargs.get("stderr_path")
        path = kwargs.get("stdout_path") or Path(str(stderr) + ".stdout.log")
        path.write_bytes(stdout)
        if stderr:
            stderr.write_bytes(b"")
        return CommandResult(tuple(args), 0, None, None, 0.001, kwargs.get("cwd"), True,
                             stdout_path=path, stderr_path=stderr)


@pytest.fixture
def native(tmp_path):
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">chr1\n" + "A" * 1000 + "\n")
    Path(str(fasta) + ".fai").write_text("chr1\t1000\t6\t1000\t1001\n")
    reference = ReferenceGenome.from_fasta(fasta, build="GRCh38")
    source = tmp_path / "caller.vcf.gz"
    source.write_bytes(gzip.compress(TEXT.encode(), mtime=0))
    source_index = Path(str(source) + ".tbi")
    source_index.write_bytes(b"source index")
    sources = (SVEvidenceSourceArtifact("S1", SVEvidenceSource.READ, "sawfish", source, source_index, EvidenceRunStatus.COMPLETED),)
    path, evidence = tmp_path / "S1.harmonized.sv.vcf.gz", tmp_path / "S1.sv.evidence.tsv"
    request = SVHarmonizationRequest("S1", reference, sources, tmp_path, path, evidence)
    path.write_bytes(gzip.compress(TEXT.encode(), mtime=0))
    request.output_index.write_bytes(b"native index")
    evidence.write_text("native evidence remains separate\n")
    raw = tmp_path / "S1.jasmine.raw.vcf"
    raw.write_bytes(TEXT.encode())
    commands = (JasmineCommandPlan("jasmine", ("jasmine", "file_list=bound_inputs", "out_file=bound_raw")),)
    artifact = HarmonizedSvArtifact("S1", fasta, "GRCh38", path, request.output_index, evidence,
        sources, "1.1.5", tuple(item.args for item in commands), intermediate_files=(raw,))
    return JasmineResult(request, JasmineResultStatus.COMPLETED, commands, {"jasmine": "1.1.5"}, 1.0, artifact)


def request_for(native, path):
    return BndRepresentationRequest(native, path, BndPolicy(POLICY["policy_id"], 1, POLICY_SHA256, POLICY["mode"]))


def test_fake_tool_bundle_is_explicit_reversible_and_non_mutating(tmp_path, native, monkeypatch):
    runner = FakeRunner(tmp_path)
    monkeypatch.setenv("ALPHAGENOME_API_KEY", "unit-test-sentinel-never-forward")
    request = request_for(native, tmp_path / "bundle")
    before = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in request.inputs()}
    result = BndRepresentationWrapper(runner=runner).run(request)
    manifest = json.loads(result.manifest_path.read_text())
    assert manifest["status"] == "completed" and manifest["source_order"] == ["sawfish"]
    assert manifest["summary"]["round_trip"]["status"] == "PASS"
    assert manifest["summary"]["compressed_round_trip"]["status"] == "PASS"
    assert manifest["summary"]["local_index_queries"] == [{"ordinal": 1, "region": "chr1:20-20", "status": "PASS"}]
    assert result.to_dict()["native_output_role_unchanged"]
    assert not result.to_dict()["standalone_derived_VCF_is_lossless"]
    assert native.artifact.vcf_path == native.request.output_vcf
    assert {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in before} == before
    assert (request.output_directory / manifest["native_binary_copies"]["vcf"]).read_bytes() == native.artifact.vcf_path.read_bytes()
    for _, kwargs in runner.calls:
        assert kwargs["inherit_env"] is False
        assert "ALPHAGENOME_API_KEY" not in kwargs["env"]
        assert set(kwargs["env"]) <= {"PATH", "SYSTEMROOT", "WINDIR", "HOME", "TMPDIR", "TMP", "TEMP", "XDG_CONFIG_HOME", "XDG_CACHE_HOME", "LC_ALL", "TZ"}
    assert "unit-test-sentinel" not in result.manifest_path.read_text()
    with pytest.raises(OutputValidationError):
        request_for(native, tmp_path / "bundle")


def test_dry_run_creates_nothing_and_claims_no_success(tmp_path, native):
    output = tmp_path / "planned"
    result = BndRepresentationWrapper().run(request_for(native, output), dry_run=True)
    assert result["status"] == "planned" and not result["will_create_outputs"]
    assert not output.exists()


@pytest.mark.parametrize("kind", ["planned", "version", "source_order", "artifact_path", "command"])
def test_unbound_or_incomplete_native_result_fails(tmp_path, native, kind):
    if kind == "planned":
        native = replace(native, status=JasmineResultStatus.PLANNED, artifact=None)
    elif kind == "version":
        native = replace(native, tool_versions={"jasmine": "2.0"})
    elif kind == "source_order":
        native = replace(native, artifact=replace(native.artifact, sources=()))
    elif kind == "artifact_path":
        native = replace(native, artifact=replace(native.artifact, evidence_table=tmp_path / "unbound.tsv"))
    else:
        native = replace(native, artifact=replace(native.artifact, commands=(("unbound",),)))
    with pytest.raises(InputValidationError):
        request_for(native, tmp_path / "bundle")


@pytest.mark.parametrize("failure", ["bgzip", "tabix", "old_version"])
def test_failure_has_no_success_sentinel_and_preserves_native(tmp_path, native, failure):
    runner = FakeRunner(tmp_path, fail=None if failure == "old_version" else failure,
                        version="1.21" if failure == "old_version" else "1.24")
    request = request_for(native, tmp_path / "failed_bundle")
    before = {p: p.read_bytes() for p in request.inputs()}
    with pytest.raises((CommandExecutionError, ToolVersionError)):
        BndRepresentationWrapper(runner=runner).run(request)
    assert not (request.output_directory / "manifest.json").exists()
    failed = json.loads((request.output_directory / "FAILED.json").read_text())
    assert failed["status"] == "failed" and not failed["success_manifest_published"]
    if failure != "old_version":
        assert failed["commands_completed"][-1]["status"] == "failed"
    assert {p: p.read_bytes() for p in before} == before


def test_race_to_existing_directory_is_rejected_before_file_writes(tmp_path, native):
    request = request_for(native, tmp_path / "late")
    request.output_directory.mkdir()
    marker = request.output_directory / "user.txt"
    marker.write_bytes(b"user data")
    with pytest.raises(FileExistsError):
        BndRepresentationWrapper().run(request)
    assert list(request.output_directory.iterdir()) == [marker]
    assert marker.read_bytes() == b"user data"


def test_missing_native_raw_and_relative_output_are_rejected(tmp_path, native):
    with pytest.raises(InputValidationError):
        request_for(native, Path("relative_bundle"))
    native = replace(native, artifact=replace(native.artifact, intermediate_files=()))
    request = request_for(native, tmp_path / "bundle")
    with pytest.raises(InputValidationError):
        request.inputs()
