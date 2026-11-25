"""Disabled producers must not claim ready inputs in their output namespace.

These are real Snakemake DAG tests, not biological execution or file-format
validation. Structural placeholders are never passed to a caller.
"""
from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys

import pytest
import yaml

from hifivar.config import load_config, write_effective_config


ROOT = Path(__file__).resolve().parents[2]
DEFAULT = ROOT / "src/hifivar/resources/configs/default.yaml"
PRODUCER_OUTPUTS = (
    ("small/S1.small.vcf.gz", "deepvariant_small"),
    ("small/S1.g.vcf.gz", "deepvariant_small"),
    ("sv/S1.sawfish.sv.vcf.gz", "read_based_sv"),
    ("tr/S1.tr.vcf.gz", "tandem_repeat"),
    ("phasing/S1.phased.vcf.gz", "hiphase_phasing"),
    ("assembly/S1/S1.hap1.fa", "hifiasm_assemble"),
    ("assembly_sv/S1/S1.pav.assembly.sv.vcf.gz", "pav_assembly_sv"),
    ("assembly_sv/S1/S1.svim_asm.assembly.sv.vcf.gz", "svim_asm_call"),
    ("sv_harmonized/S1/S1.harmonized.sv.vcf.gz", "harmonize_sv"),
)


def _configuration(tmp_path: Path, **overrides) -> Path:
    (tmp_path / "ref.fa").write_text(">chr1\nACGT\n", encoding="utf-8")
    (tmp_path / "ref.fa.fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")
    bam = tmp_path / "S1.bam"
    bam.write_bytes(b"structural-DAG-fixture-not-a-biological-BAM")
    Path(f"{bam}.bai").write_bytes(b"structural-index-fixture")
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(f"sample_id\tinput\tinput_type\nS1\t{bam}\tbam\n", encoding="utf-8")
    settings = {
        "reference": {"fasta": str(tmp_path / "ref.fa"), "build": "GRCh38"},
        "samples": {"sheet": str(sheet)},
        "paths": {"workdir": str(tmp_path / "work"), "outdir": str(tmp_path / "results")},
    }
    settings.update(overrides)
    user = tmp_path / "user.yaml"
    user.write_text(yaml.safe_dump(settings), encoding="utf-8")
    effective = tmp_path / "effective.yaml"
    write_effective_config(load_config(DEFAULT, user_config=user), effective)
    return effective


def _scheduler(tmp_path: Path, config: Path, *arguments: str) -> subprocess.CompletedProcess:
    pytest.importorskip("snakemake")
    environment = os.environ.copy()
    environment.update(
        PYTHONPATH=str(ROOT / "src"), PYTHONDONTWRITEBYTECODE="1",
        TMP=str(tmp_path), TEMP=str(tmp_path), TMPDIR=str(tmp_path),
        XDG_CACHE_HOME=str(tmp_path / "cache"), XDG_CONFIG_HOME=str(tmp_path / "config"),
        LOCALAPPDATA=str(tmp_path / "local"), APPDATA=str(tmp_path / "app"),
    )
    return subprocess.run(
        [sys.executable, "-B", "-m", "snakemake", "--snakefile", str(tmp_path / "Snakefile"),
         "--configfile", str(config), "--cores", "2", "--shared-fs-usage", "input-output",
         "persistence", "software-deployment", "software-deployment-cache", "sources",
         "storage-local-copies", *arguments],
        cwd=tmp_path, env=environment, capture_output=True, text=True,
        encoding="utf-8", errors="replace", timeout=60, check=False,
    )


def _consumer(tmp_path: Path, source: Path) -> None:
    (tmp_path / "Snakefile").write_text(
        f"include: {json.dumps((ROOT / 'workflow/Snakefile').as_posix())}\n"
        "rule consume_ready:\n"
        f"    input: {json.dumps(source.as_posix())}\n"
        "    output: 'consumer.done'\n"
        "    shell: 'touch {output}'\n",
        encoding="utf-8",
    )


@pytest.mark.parametrize("relative,producer", PRODUCER_OUTPUTS)
def test_disabled_producer_does_not_claim_ready_in_namespace(tmp_path, relative, producer):
    effective = _configuration(tmp_path)
    source = tmp_path / "results" / relative
    source.parent.mkdir(parents=True)
    source.write_bytes(b"ready-input-structural-placeholder")
    before = source.read_bytes()
    _consumer(tmp_path, source)
    result = _scheduler(tmp_path, effective, "--dry-run", "--", "consume_ready")
    output = result.stdout + result.stderr
    assert result.returncode == 0, output
    assert "rule consume_ready:" in output
    assert f"rule {producer}:" not in output
    assert source.read_bytes() == before
    assert not (tmp_path / "consumer.done").exists()
    registered = _scheduler(tmp_path, effective, "--list-rules")
    assert registered.returncode == 0, registered.stdout + registered.stderr
    assert producer not in registered.stdout.split()


@pytest.mark.parametrize("relative", ["small/S1.small.vcf.gz", "sv/S1.sawfish.sv.vcf.gz"])
def test_missing_ready_input_fails_without_reenabling_producer(tmp_path, relative):
    effective = _configuration(tmp_path)
    source = tmp_path / "results" / relative
    _consumer(tmp_path, source)
    result = _scheduler(tmp_path, effective, "--dry-run", "--", "consume_ready")
    output = result.stdout + result.stderr
    assert result.returncode != 0
    assert "MissingInputException" in output
    assert "InputFunctionException" not in output
    assert not source.exists()


@pytest.mark.parametrize("ready", [False, True])
def test_disabled_individual_sv_caller_is_not_eligible(tmp_path, ready):
    effective = _configuration(tmp_path, sv={
        "enabled": True, "sawfish": {"enabled": True},
        "sniffles2": {"enabled": False}, "pbsv": {"enabled": False},
        "cutesv": {"enabled": False},
    })
    source = tmp_path / "results/sv/S1.pbsv.sv.vcf.gz"
    if ready:
        source.parent.mkdir(parents=True)
        source.write_bytes(b"ready-disabled-caller-structural-placeholder")
    _consumer(tmp_path, source)
    result = _scheduler(tmp_path, effective, "--dry-run", "--", "consume_ready")
    output = result.stdout + result.stderr
    if ready:
        assert result.returncode == 0, output
        assert "rule read_based_sv:" not in output
    else:
        assert result.returncode != 0
        assert "MissingInputException" in output
    assert not (tmp_path / "consumer.done").exists()


def test_external_harmonization_still_works_with_read_calling_disabled(tmp_path):
    source = tmp_path / "results/sv/S1.sawfish.sv.vcf.gz"
    source.parent.mkdir(parents=True)
    source.write_bytes(b"external-harmonization-structural-placeholder")
    Path(f"{source}.tbi").write_bytes(b"structural-index-fixture")
    effective = _configuration(tmp_path, sv={"enabled": False, "harmonization": {
        "enabled": True, "input_vcfs": {"S1": {"sawfish": str(source)}},
    }})
    (tmp_path / "Snakefile").write_text(
        f"include: {json.dumps((ROOT / 'workflow/Snakefile').as_posix())}\n", encoding="utf-8",
    )
    result = _scheduler(tmp_path, effective, "--dry-run", "--", "all")
    output = result.stdout + result.stderr
    assert result.returncode == 0, output
    assert "rule harmonize_sv:" in output
    assert "rule read_based_sv:" not in output

@pytest.mark.parametrize('policy_kind',['bnd','identity','both'])
def test_generated_consumer_snakemake_disables_standalone_harmonization_without_deleting_bundle(tmp_path,policy_kind):
    from hifivar.breakends import BndPolicy,POLICY,POLICY_SHA256
    from hifivar.sv_input_identity import POLICY as IDENTITY_POLICY,POLICY_SHA256 as IDENTITY_SHA256
    from hifivar.downstream_handoff import consumer_config
    source=tmp_path/'results/sv/S1.sawfish.sv.vcf.gz';source.parent.mkdir(parents=True);source.write_bytes(b'structural existing source')
    section={'enabled':True,'input_vcfs':{'S1':{'sawfish':str(source)}}}
    if policy_kind in ('bnd','both'):
        section['bnd_representation']=BndPolicy(POLICY['policy_id'],1,POLICY_SHA256,POLICY['mode']).to_dict()
    if policy_kind in ('identity','both'):
        section['input_identity']={'enabled':True,'policy_id':IDENTITY_POLICY['policy_id'],'policy_version':1,
            'policy_sha256':IDENTITY_SHA256,'mode':IDENTITY_POLICY['mode'],'timeout_seconds':300}
    effective=_configuration(tmp_path,sv={'enabled':False,'harmonization':section})
    base=tmp_path/'results/sv_harmonized/S1'
    bnd=base/'bndrepr-v1';bnd.mkdir(parents=True)
    identity=base/'source-identity-v1';identity.mkdir()
    ready=bnd/'S1.bndrepr-v1.sv.vcf.gz';ready.write_bytes(b'ready structural derived source')
    protected={source:source.read_bytes(),ready:ready.read_bytes()}
    _consumer(tmp_path,ready)
    previous=_scheduler(tmp_path,effective,'--dry-run','--','consume_ready')
    assert previous.returncode!=0 and 'fresh per-sample bundle' in previous.stdout+previous.stderr
    original=yaml.safe_load(effective.read_text(encoding='utf-8'))
    safe=consumer_config(original)
    assert original['sv']['harmonization']['enabled'] is True
    assert safe['sv']['harmonization']['enabled'] is False
    assert safe['sv']['harmonization']['bnd_representation']==safe['sv']['harmonization']['input_identity']=={'enabled':False}
    prepared=tmp_path/'consumer.effective.yaml';prepared.write_text(yaml.safe_dump(safe),encoding='utf-8')
    result=_scheduler(tmp_path,prepared,'--dry-run','--','consume_ready')
    assert result.returncode==0,result.stdout+result.stderr
    assert 'rule consume_ready:' in result.stdout+result.stderr
    assert 'rule harmonize_sv:' not in result.stdout+result.stderr
    assert not (tmp_path/'consumer.done').exists()
    assert all(p.read_bytes()==v for p,v in protected.items()) and bnd.is_dir() and identity.is_dir()
