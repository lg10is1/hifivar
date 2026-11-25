from __future__ import annotations

from dataclasses import replace
import gzip
import hashlib
import os
from pathlib import Path

import pytest
import yaml

from hifivar.config import load_config
from hifivar.exceptions import ConfigurationError, InputValidationError, OutputValidationError, ReferenceError
from hifivar.hiphase import HiPhaseWrapper
from hifivar.joint_phasing import PhasingTrackInput, phasing_options, validate_phasing_producers, validate_record_preservation
from hifivar.phasing import PhasingRequest
from test_hiphase import FakeRunner, make_request, write_bgzf


def header(category: str, sample: str = "S1", length: int = 8) -> str:
    specific = {
        "small": "##source=DeepVariant\n",
        "sv": '##source=pbsv\n##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">\n',
        "tr": ('##source=TRGT\n##INFO=<ID=TRID,Number=1,Type=String,Description="ID">\n'
               '##INFO=<ID=MOTIFS,Number=.,Type=String,Description="Motifs">\n'
               '##INFO=<ID=STRUC,Number=1,Type=String,Description="Structure">\n'),
    }[category]
    return ('##fileformat=VCFv4.2\n##reference=GRCh38\n'
            f'##contig=<ID=chr1,length={length}>\n{specific}'
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + sample + '\n')


def joint_request(tmp_path: Path) -> PhasingRequest:
    base = make_request(tmp_path)
    write_bgzf(base.small_variants.vcf_path,
               (header("small") + 'chr1\t2\tsmall1\tC\tT\t30\tPASS\t.\tGT\t0/1\n').encode())
    tracks = []
    for category, tool in (("sv", "pbsv"), ("tr", "trgt")):
        path = tmp_path / ("S1.pbsv.sv.vcf.gz" if category == "sv" else "S1.tr.vcf.gz")
        info = "SVTYPE=DEL" if category == "sv" else "TRID=L1;MOTIFS=C;STRUC=(C)n"
        write_bgzf(path, (header(category) + f'chr1\t3\t{category}1\tG\tA\t30\tPASS\t{info}\tGT\t0/1\n').encode())
        write_bgzf(Path(f"{path}.tbi"), b"TBI\x01")
        tracks.append(PhasingTrackInput(category, tool, "S1", "GRCh38", path, Path(f"{path}.tbi")))
    write_bgzf(base.small_variants.vcf_index_path, b"TBI\x01")
    return replace(base, additional_inputs=tuple(tracks))


class JointRunner(FakeRunner):
    def __init__(self, corrupt: str = ""):
        super().__init__()
        self.corrupt = corrupt

    def run(self, command, **kwargs):
        args = tuple(str(p) for p in command)
        result = super().run(command, **kwargs)
        if kwargs.get("dry_run") or "--version" in args:
            return result
        if args[0] == "hiphase":
            sources = [args[i + 1] for i, v in enumerate(args) if v == "--vcf"]
            outputs = [args[i + 1] for i, v in enumerate(args) if v == "--output-vcf"]
            for n, (source, output) in enumerate(zip(sources, outputs)):
                rows = []
                with gzip.open(source, "rt") as inp:
                    for line in inp:
                        if line.startswith("#CHROM"):
                            if not (self.corrupt == "no_ps" and n == 2):
                                rows.append('##FORMAT=<ID=PS,Number=1,Type=Integer,Description="PS">\n')
                        if not line.startswith("#"):
                            if self.corrupt == "drop" and n == 1:
                                continue
                            f = line.rstrip().split("\t")
                            f[8] += ":PS"
                            f[9] = f[9].replace("/", "|") + ":2"
                            if self.corrupt == "gt" and n == 1:
                                f[9] = "1|1:2"
                            line = "\t".join(f) + "\n"
                        rows.append(line)
                write_bgzf(Path(output), "".join(rows).encode())
        elif args[0] == "tabix":
            write_bgzf(Path(args[-1] + ".tbi"), b"TBI\x01")
        return result


def test_joint_flags_output_order_global_auto_and_shell_false(tmp_path):
    req = joint_request(tmp_path)
    commands = HiPhaseWrapper().plan_commands(req)
    args = commands[0].args
    assert args.count("--vcf") == args.count("--output-vcf") == 3
    assert [p[0] for p in req.vcf_pairs] == ["small", "sv", "tr"]
    assert [args[i + 1] for i, value in enumerate(args) if value == "--output-vcf"] == [str(p[3].absolute()) for p in req.vcf_pairs]
    assert "--disable-global-realignment" not in args
    assert len(commands) == 4 and all(c.to_dict()["shell"] is False for c in commands)
    assert "--disable-global-realignment" in HiPhaseWrapper().plan_commands(replace(req, global_realignment="disabled"))[0].args


def test_joint_dry_run_writes_no_output_and_needs_no_tools(tmp_path):
    req = joint_request(tmp_path)
    result = HiPhaseWrapper(runner=JointRunner()).run(req, dry_run=True)
    assert result.artifact is None and not req.output_vcf.parent.exists()


def test_joint_fake_execution_preserves_three_sources_and_provenance(tmp_path):
    req = joint_request(tmp_path)
    before = {str(p): hashlib.sha256(p.read_bytes()).hexdigest()
              for _, source, index, _, _, _ in req.vcf_pairs for p in (source, index)}
    result = HiPhaseWrapper(runner=JointRunner()).run(req)
    assert result.artifact.records_identity_and_genotype_checked == 1
    assert len(result.artifact.additional_artifacts) == 2
    assert all(a.records_identity_and_genotype_checked == 1 for a in result.artifact.additional_artifacts)
    assert all(a.to_dict()["phase_accuracy"] == "NOT_ASSESSED" for a in result.artifact.additional_artifacts)
    assert all(hashlib.sha256(Path(p).read_bytes()).hexdigest() == digest for p, digest in before.items())
    assert result.to_dict()["request"]["effective_global_realignment"] is True


@pytest.mark.parametrize("corrupt,match", [("drop", "number"), ("gt", "genotype"), ("no_ps", "FORMAT/PS")])
def test_joint_rejects_partial_or_mutated_outputs(tmp_path, corrupt, match):
    with pytest.raises(OutputValidationError, match=match):
        HiPhaseWrapper(runner=JointRunner(corrupt)).run(joint_request(tmp_path))


@pytest.mark.parametrize("caller", ["sniffles2", "cutesv", "pav", "svim_asm", "unknown"])
def test_unsupported_sv_producers_fail_before_execution(tmp_path, caller):
    with pytest.raises(InputValidationError):
        PhasingTrackInput("sv", caller, "S1", "GRCh38", tmp_path / "x.vcf.gz", tmp_path / "x.vcf.gz.tbi")


def test_cross_sample_build_and_duplicate_categories_fail(tmp_path):
    req = joint_request(tmp_path)
    with pytest.raises(InputValidationError, match="sample"):
        replace(req, additional_inputs=(replace(req.additional_inputs[0], sample_id="S2"),))
    with pytest.raises(ReferenceError, match="build"):
        replace(req, additional_inputs=(replace(req.additional_inputs[0], reference_build="GRCh37"),))
    with pytest.raises(InputValidationError, match="one explicitly"):
        replace(req, additional_inputs=(req.additional_inputs[0], req.additional_inputs[0]))


def test_missing_index_and_wrong_actual_header_fail(tmp_path):
    req = joint_request(tmp_path)
    source = req.additional_inputs[0]
    write_bgzf(source.vcf_path, (header("sv", sample="S2") + 'chr1\t3\tx\tG\tA\t30\tPASS\tSVTYPE=DEL\tGT\t0/1\n').encode())
    with pytest.raises((InputValidationError, OutputValidationError), match="(?i)sample"):
        HiPhaseWrapper().plan_commands(req)
    write_bgzf(source.vcf_path, header("sv").encode())
    source.index_path.unlink()
    with pytest.raises((InputValidationError, OutputValidationError, OSError)):
        HiPhaseWrapper().plan_commands(req)


def test_actual_contig_length_mismatch_fails(tmp_path):
    req = joint_request(tmp_path)
    write_bgzf(req.additional_inputs[0].vcf_path, header("sv", length=9).encode())
    with pytest.raises((InputValidationError, ReferenceError), match="length|LENGTH|contig|CONTIG"):
        HiPhaseWrapper().plan_commands(req)


def test_joint_checks_all_existing_outputs_even_after_request_creation(tmp_path):
    req = joint_request(tmp_path)
    extra_output = req.vcf_pairs[2][3]
    extra_output.parent.mkdir(parents=True)
    extra_output.write_bytes(b"existing")
    with pytest.raises(OutputValidationError, match="already exists"):
        HiPhaseWrapper().run(req)
    assert extra_output.read_bytes() == b"existing"


def test_overwrite_true_never_replaces_input_hardlink(tmp_path):
    req = joint_request(tmp_path)
    out = req.vcf_pairs[1][3]
    out.parent.mkdir(parents=True)
    os.link(req.additional_inputs[0].vcf_path, out)
    with pytest.raises(InputValidationError, match="source artifact"):
        replace(req, overwrite=True)


@pytest.mark.parametrize("bad", [{"sv_caller": "sniffles2"}, {"include_tr": "true"},
                                 {"global_realignment": True}, {"global_realignment": "guess"}])
def test_invalid_options_fail_closed(bad):
    with pytest.raises(ConfigurationError):
        phasing_options(bad)


def test_producer_dependencies_do_not_reenable_disabled_callers():
    config = {"small": {"enabled": True}, "phasing": {"enabled": True, "sv_caller": "pbsv"},
              "sv": {"enabled": True, "pbsv": {"enabled": False}}, "tr": {"enabled": False}}
    with pytest.raises(ConfigurationError, match="enabled SV producer"):
        validate_phasing_producers(config)
    assert config["sv"]["pbsv"]["enabled"] is False
    config["phasing"] = {"enabled": True, "include_tr": True}
    with pytest.raises(ConfigurationError, match="tr.enabled"):
        validate_phasing_producers(config)


def test_config_schema_accepts_new_settings_without_unknown_warning(tmp_path):
    config_path = tmp_path / "joint.yaml"
    config_path.write_text(yaml.safe_dump({"phasing": {"sv_caller": "pbsv", "include_tr": True,
                                                      "global_realignment": "auto"}}))
    default = Path(__file__).resolve().parents[2] / "src/hifivar/resources/configs/default.yaml"
    cfg = load_config(default, user_config=config_path)
    assert cfg["phasing"]["sv_caller"] == "pbsv" and cfg["phasing"]["include_tr"]


def test_record_identity_change_rejected_independently_of_gt(tmp_path):
    req = joint_request(tmp_path)
    source = req.additional_inputs[0].vcf_path
    output = tmp_path / "different.vcf.gz"
    with gzip.open(source, "rt") as handle:
        text = handle.read().replace("\tsv1\t", "\tchanged\t")
    write_bgzf(output, text.encode())
    with pytest.raises(OutputValidationError, match="identity"):
        validate_record_preservation(source, output)


def test_phase6_api_routes_real_artifact_types_and_rejects_missing_selected_track(tmp_path):
    from hifivar.context import AnalysisContext
    from hifivar.phase6 import run_phase6
    from hifivar.sv import StructuralVariantArtifact, SvCaller
    from hifivar.tr import TandemRepeatArtifact, TandemRepeatCatalog

    req = joint_request(tmp_path)
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(f"sample_id\tinput\tinput_type\nS1\t{req.alignment.path}\tbam\n", encoding="utf-8")
    user = tmp_path / "api.yaml"
    user.write_text(yaml.safe_dump({"reference": {"fasta": str(req.alignment.reference.fasta), "build": "GRCh38"},
        "samples": {"sheet": str(sheet)}, "phasing": {"enabled": True, "sv_caller": "pbsv", "include_tr": True}}))
    default = Path(__file__).resolve().parents[2] / "src/hifivar/resources/configs/default.yaml"
    cfg = load_config(default, user_config=user)
    context = AnalysisContext.from_config(cfg)
    sv, tr = req.additional_inputs
    sv_artifact = StructuralVariantArtifact(SvCaller.PBSV, "S1", context.reference.fasta, "GRCh38",
        sv.vcf_path, sv.index_path, "2.11.0", ())
    tr_artifact = TandemRepeatArtifact("S1", TandemRepeatCatalog(tmp_path / "test.bed", "GRCh38"),
        "GRCh38", tr.vcf_path, tr.index_path, req.alignment.path, req.alignment.index_path, "5.1.0", ())
    common = dict(alignment_artifacts={"S1": req.alignment}, small_variant_artifacts={"S1": req.small_variants},
                  config=cfg, wrapper=HiPhaseWrapper(runner=JointRunner()))
    with pytest.raises(InputValidationError, match="selected SV artifact"):
        run_phase6(context, output_directory=tmp_path / "missing", **common)
    report = run_phase6(context, sv_variant_artifacts={"S1": sv_artifact}, tandem_repeat_artifacts={"S1": tr_artifact},
                        output_directory=tmp_path / "api_results", **common)
    assert len(report.sample_results[0].phasing.artifact.additional_artifacts) == 2
    assert report.settings.global_realignment == "auto"
    assert report.to_dict()["settings"]["sv_caller"] == "pbsv"


@pytest.mark.parametrize("role", ["reference", "sv", "tr", "index"])
def test_joint_log_source_alias_is_rejected_before_version_command(tmp_path, role):
    req = joint_request(tmp_path)
    source = {"reference": req.alignment.reference.fasta,
              "sv": req.additional_inputs[0].vcf_path,
              "tr": req.additional_inputs[1].vcf_path,
              "index": req.additional_inputs[0].index_path}[role]
    before = source.read_bytes()
    runner = JointRunner()
    with pytest.raises(InputValidationError, match="log path aliases"):
        HiPhaseWrapper(runner=runner).run(req, stderr_path=source)
    assert runner.calls == [] and source.read_bytes() == before


def test_implicit_numbered_stdout_log_hardlink_is_rejected(tmp_path):
    req = joint_request(tmp_path)
    base = tmp_path / "logs/hiphase.log"
    base.parent.mkdir()
    actual_stdout = base.parent / "hiphase.1.stdout.log"
    os.link(req.alignment.reference.fasta, actual_stdout)
    runner = JointRunner()
    with pytest.raises(InputValidationError, match="log path aliases"):
        HiPhaseWrapper(runner=runner).run(req, stderr_path=base)
    assert runner.calls == []
