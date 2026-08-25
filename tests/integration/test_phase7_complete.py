from __future__ import annotations

import json
from pathlib import Path

import yaml

from hifivar.command import CommandResult
from hifivar.config import load_config
from hifivar.context import AnalysisContext
from hifivar.hifiasm import HifiasmWrapper
from hifivar.phase7 import Phase7RunStatus, run_phase7


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = ROOT / "src" / "hifivar" / "resources" / "configs" / "default.yaml"


class Phase7Runner:
    def require_executable(self, executable: str) -> Path:
        return Path(executable)

    def run(self, command, **kwargs):
        args = tuple(str(value) for value in command)
        if "--version" in args:
            return CommandResult(
                args,
                0,
                "hifiasm version 0.25.0-r726",
                "",
                0.01,
                None,
                True,
            )
        prefix = Path(args[args.index("-o") + 1])
        prefix.parent.mkdir(parents=True, exist_ok=True)
        for suffix, sequence in (
            (".bp.p_ctg.gfa", "ACGT"),
            (".bp.hap1.p_ctg.gfa", "AAAA"),
            (".bp.hap2.p_ctg.gfa", "TTTT"),
        ):
            Path(f"{prefix}{suffix}").write_text(
                f"H\tVN:Z:1.0\nS\tctg1\t{sequence}\n",
                encoding="utf-8",
            )
        return CommandResult(args, 0, "", "", 1.0, None, True)


def test_fastq_context_to_haplotype_assemblies_and_provenance(tmp_path: Path) -> None:
    fasta = tmp_path / "reference.fa"
    fasta.write_text(">chr1\nACGTACGT\n", encoding="utf-8")
    Path(f"{fasta}.fai").write_text("chr1\t8\t6\t8\t9\n", encoding="utf-8")
    first = tmp_path / "reads1.fastq"
    second = tmp_path / "reads2.fastq"
    first.write_text("@one\nACGT\n+\n!!!!\n", encoding="utf-8")
    second.write_text("@two\nTGCA\n+\n!!!!\n", encoding="utf-8")
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        f"sample_id\tinput\tinput_type\nS1\t{first};{second}\tfastq\n",
        encoding="utf-8",
    )
    user = tmp_path / "phase7.yaml"
    user.write_text(
        yaml.safe_dump(
            {
                "reference": {"fasta": str(fasta), "build": "GRCh38"},
                "samples": {"sheet": str(sheet)},
                "assembly": {"enabled": True},
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )
    config = load_config(DEFAULT_CONFIG, user_config=user)
    context = AnalysisContext.from_config(config)
    report = run_phase7(
        context,
        output_directory=tmp_path / "results" / "assembly",
        work_directory=tmp_path / "work" / "hifiasm",
        config=config,
        wrapper=HifiasmWrapper(runner=Phase7Runner()),  # type: ignore[arg-type]
    )
    assert report.status is Phase7RunStatus.COMPLETED
    artifact = report.sample_results[0].assembly.artifact
    assert artifact is not None
    assert len(artifact.raw_gfas) == 3
    assert len(artifact.assemblies) == 3
    assert report.sample_results[0].assembly.command.args[-2:] == (
        str(first.absolute()),
        str(second.absolute()),
    )
    json_path = report.write_json(tmp_path / "phase7-report.json")
    yaml_path = report.write_yaml(tmp_path / "phase7-report.yaml")
    assert json.loads(json_path.read_text(encoding="utf-8")) == yaml.safe_load(
        yaml_path.read_text(encoding="utf-8")
    )
