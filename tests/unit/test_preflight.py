"""Synthetic preflight contracts, not real tool or scientific acceptance."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace
import pytest
import yaml
from hifivar import preflight
from hifivar.cli import main
from hifivar.config import HiFiVarConfig, load_config
from hifivar.exceptions import CommandExecutionError, ConfigurationError, ToolNotFoundError

DEFAULT = Path(__file__).parents[2] / "src/hifivar/resources/configs/default.yaml"
HEADER = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@RG\tID:g\tSM:SYNTHETIC\tPL:PACBIO\n"


class FakeRunner:
    def __init__(self, *, missing=(), header=HEADER, fail_version=None):
        self.missing, self.header, self.fail_version, self.commands = set(missing), header, fail_version, []

    def find_executable(self, executable):
        return None if executable in self.missing else Path(executable)

    def require_executable(self, executable):
        if executable in self.missing:
            raise ToolNotFoundError("synthetic missing executable")
        return Path(executable)

    def run(self, command, **kwargs):
        self.commands.append((command, kwargs))
        if command[0] == self.fail_version:
            raise CommandExecutionError("synthetic timeout")
        if "snakemake" in command:
            text = "9.27.0\n"
        elif command[0] == "samtools":
            text = self.header if "-H" in command else ("0\n" if "-X" in command else "samtools 1.23.1\n")
        else:
            outputs = {"run_deepvariant": "DeepVariant version 1.10.0", "pbmm2": "pbmm2 1.17.0",
                       "minimap2": "2.31-r1302", "sawfish": "sawfish 2.2.1", "sniffles": "Sniffles2 2.8.0",
                       "pbsv": "pbsv 2.11.0", "cuteSV": "cuteSV 2.1.4", "bgzip": "bgzip (htslib) 1.23.1",
                       "tabix": "tabix (htslib) 1.23.1", "trgt": "trgt 1.0.0", "bcftools": "bcftools 1.23.1",
                       "hiphase": "HiPhase 1.4.0"}
            outputs.update(jasmine="Jasmine 1.1.5", truvari="Truvari 5.4.0",
                           hifiasm="0.25.0-r726", **{"svim-asm": "SVIM-asm 1.0.3",
                           "vep": "ensembl-vep version 116.3", "glnexus_cli": "glnexus_cli 1.4.1"})
            text = outputs.get(command[0], "DeepVariant version 1.10.0")
        return SimpleNamespace(args=tuple(command), returncode=0, stdout=text, stderr="", duration_seconds=0.01)


@pytest.fixture
def case(tmp_path, monkeypatch):
    (tmp_path / "ref.fa").write_text(">chr1\n" + "A" * 1000 + "\n", encoding="utf-8")
    (tmp_path / "ref.fa.fai").write_text("chr1\t1000\t6\t1000\t1001\n", encoding="utf-8")
    (tmp_path / "SYNTHETIC.bam").write_bytes(b"synthetic, not real BAM")
    (tmp_path / "SYNTHETIC.bam.bai").write_bytes(b"synthetic, not real index")
    (tmp_path / "samples.tsv").write_text("sample_id\tinput\tinput_type\tsex\nSYNTHETIC\tSYNTHETIC.bam\tbam\tfemale\n", encoding="utf-8")
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump({"reference": {"fasta": "ref.fa", "build": "CUSTOM_TEST"},
        "samples": {"sheet": "samples.tsv"}, "runtime": {"threads": 32},
        "small": {"enabled": True, "threads": 4, "memory_mb": 100}}), encoding="utf-8")
    monkeypatch.setattr(preflight.platform, "system", lambda: "Linux")
    monkeypatch.setattr(preflight.importlib.util, "find_spec", lambda name: object())
    return load_config(default_config=DEFAULT, user_config=path), tmp_path, path


def change(config, **sections):
    values = config.to_dict()
    for section, changes in sections.items():
        values[section].update(changes)
    return HiFiVarConfig(values, config.sources)


def hashes(root):
    return {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest() for p in root.rglob("*") if p.is_file()}


def check(report, key):
    return next(row for row in report["checks"] if row["key"] == key)


def test_metadata_never_executes_or_creates(case):
    cfg, root, _ = case
    before = hashes(root); runner = FakeRunner()
    result = preflight.run_preflight(cfg, root / "planned", runner=runner)
    assert result["status"] == "INCOMPLETE" and result["exit_code"] == 2
    assert not runner.commands and hashes(root) == before and not (root / "planned").exists()
    assert check(result, "input.SYNTHETIC")["status"] == "NOT_CHECKED"


def test_explicit_probes_use_existing_wrappers_without_files(case, monkeypatch):
    cfg, root, _ = case
    before = hashes(root); runner = FakeRunner()
    monkeypatch.setenv("PREFLIGHT_PRIVATE_TOKEN", "must-not-be-inherited")
    result = preflight.run_preflight(cfg, root / "planned", probe_tools=True, runner=runner, timeout_seconds=7)
    assert result["status"] == "CHECKED_SCOPE_PASS"
    assert check(result, "tool.deepvariant.version")["details"]["version"] == "1.10.0"
    assert check(result, "input.SYNTHETIC")["details"]["indexed_access"] == "PASS"
    assert hashes(root) == before and not (root / "planned").exists()
    for command, options in runner.commands:
        assert options["timeout"] == 7 and not options.get("stdout_path") and not options.get("stderr_path")
        assert options["inherit_env"] is False and "PREFLIGHT_PRIVATE_TOKEN" not in options["env"]
        assert "align" not in command and "genotype" not in command
    assert check(result, "resource_headroom")["status"] == "NOT_CHECKED"
    assert not result["ready_for_all_workflows"] and not result["analysis_executed"]
    assert not result["scientific_accuracy_verified"] and not result["raw_content_immutability_verified"]


@pytest.mark.parametrize("missing", ["samtools", "run_deepvariant"])
def test_missing_enabled_tool_fail_disabled_not_probed(case, missing):
    cfg, root, _ = case
    runner = FakeRunner(missing=[missing, "sawfish", "trgt", "pbmm2"])
    result = preflight.run_preflight(cfg, root / "planned", probe_tools=True, runner=runner)
    assert result["status"] == "BLOCKED"
    assert not any(c[0][0] in {"sawfish", "trgt", "pbmm2"} for c in runner.commands)


def test_probe_failure_isolated(case):
    cfg, root, _ = case
    result = preflight.run_preflight(cfg, root / "planned", probe_tools=True, runner=FakeRunner(fail_version="run_deepvariant"))
    assert check(result, "tool.deepvariant.version")["status"] == "FAIL"
    assert check(result, "input.SYNTHETIC")["status"] == "PASS"


@pytest.mark.parametrize("header", [HEADER.replace("SM:SYNTHETIC", "SM:OTHER"), HEADER.replace("LN:1000", "LN:999"),
                                   HEADER.replace("SN:chr1", "SN:1"), HEADER.replace("coordinate", "unsorted")])
def test_bam_mismatch_no_repair(case, header):
    cfg, root, _ = case
    before = hashes(root)
    result = preflight.run_preflight(cfg, root / "planned", probe_tools=True, runner=FakeRunner(header=header))
    assert check(result, "input.SYNTHETIC")["status"] == "FAIL" and hashes(root) == before


@pytest.mark.parametrize("section", ["assembly", "assembly_sv", "review", "annotation", "cohort", "benchmark", "low_depth"])
def test_unsupported_scope_before_probes(case, section):
    cfg, root, _ = case
    runner = FakeRunner()
    result = preflight.run_preflight(change(cfg, **{section: {"enabled": True}}), root / "planned", probe_tools=True, runner=runner)
    assert result["status"] == "BLOCKED" and not runner.commands
    assert section in check(result, "launcher_contract")["message"]


@pytest.mark.parametrize("values", [{"threads": 33}, {"memory_mb": 64001}, {"overwrite": True}])
def test_budget_fail_no_scaling(case, values):
    cfg, root, _ = case
    assert preflight.run_preflight(change(cfg, small=values), root / "planned", runner=FakeRunner())["status"] == "BLOCKED"


@pytest.mark.parametrize("value", [0, -1, float("nan"), float("inf"), True, "30"])
def test_bad_timeout(case, value):
    cfg, root, _ = case
    with pytest.raises(ConfigurationError):
        preflight.run_preflight(cfg, root / "planned", timeout_seconds=value)


def test_existing_destination_preserved(case):
    cfg, root, _ = case
    output = root / "existing"; output.mkdir(); (output / "user.txt").write_text("keep", encoding="utf-8")
    before = hashes(root)
    assert preflight.run_preflight(cfg, output, runner=FakeRunner())["status"] == "BLOCKED"
    assert hashes(root) == before


def test_windows_not_relabelled_linux(case, monkeypatch):
    cfg, root, _ = case
    monkeypatch.setattr(preflight.platform, "system", lambda: "Windows")
    assert check(preflight.run_preflight(cfg, root / "planned", runner=FakeRunner()), "platform")["status"] == "FAIL"


def test_container_not_started_or_pulled(case):
    cfg, root, _ = case
    image = root / "existing.sif"; image.write_bytes(b"synthetic placeholder")
    cfg = change(cfg, small={"execution_mode": "apptainer", "deepvariant_image": str(image)})
    runner = FakeRunner()
    result = preflight.run_preflight(cfg, root / "planned", probe_tools=True, runner=runner)
    assert check(result, "tool.deepvariant.version")["status"] == "NOT_CHECKED"
    assert result["status"] == "INCOMPLETE" and not any(c[0][0] == "apptainer" for c in runner.commands)


def test_fastq_route_and_pbmm2_version(case):
    cfg, root, _ = case
    (root / "reads.fastq").write_text("@movie/1/ccs\nACGT\n+\nIIII\n", encoding="utf-8")
    (root / "samples.tsv").write_text("sample_id\tinput\tinput_type\nSYNTHETIC\treads.fastq\tfastq\n", encoding="utf-8")
    result = preflight.run_preflight(change(cfg, alignment={"enabled": True}), root / "planned", probe_tools=True, runner=FakeRunner())
    assert result["status"] == "CHECKED_SCOPE_PASS"
    assert check(result, "tool.pbmm2.version")["status"] == "PASS"


def test_ubam_not_fully_scanned(case):
    cfg, root, _ = case
    (root / "SYNTHETIC.bam.bai").unlink()
    runner = FakeRunner(header="@HD\tVN:1.6\tSO:unknown\n@RG\tID:g\tSM:SYNTHETIC\tPL:PACBIO\n")
    result = preflight.run_preflight(change(cfg, alignment={"enabled": True}), root / "planned", probe_tools=True, runner=runner)
    assert result["status"] == "INCOMPLETE" and check(result, "input.SYNTHETIC")["details"]["kind"] == "UBAM"
    assert not any("-c" in command for command, _ in runner.commands)


def test_all_native_tools_and_tr_metadata(case):
    cfg, root, _ = case
    catalog = root / "repeats.bed"; catalog.write_text("chr1\t10\t20\tID=L1;MOTIFS=A;STRUC=(A)n\n", encoding="utf-8")
    cfg = change(cfg, sv={"enabled": True}, phasing={"enabled": True},
                 tr={"enabled": True, "catalog": str(catalog), "catalog_reference_build": "CUSTOM_TEST"})
    result = preflight.run_preflight(cfg, root / "planned", probe_tools=True, runner=FakeRunner())
    assert result["status"] == "CHECKED_SCOPE_PASS", [c for c in result["checks"] if c["status"] == "FAIL"]
    for name in ["sawfish", "sniffles2", "pbsv", "cutesv", "trgt", "hiphase", "sv_finalization"]:
        assert check(result, f"tool.{name}.version")["status"] == "PASS"


def test_missing_snakemake(case, monkeypatch):
    cfg, root, _ = case
    monkeypatch.setattr(preflight.importlib.util, "find_spec", lambda name: None)
    assert check(preflight.run_preflight(cfg, root / "planned", runner=FakeRunner()), "snakemake_package")["status"] == "FAIL"


def test_cli_json_zero_write(case, capsys, monkeypatch):
    _, root, path = case
    monkeypatch.setattr(preflight, "CommandRunner", FakeRunner)
    before = hashes(root)
    assert main(["--config", str(path), "doctor", "--preflight", "--run-dir", str(root / "planned"), "--json"]) == 2
    assert json.loads(capsys.readouterr().out)["status"] == "INCOMPLETE" and hashes(root) == before


@pytest.mark.parametrize("args", [["--probe-tools"], ["--json"], ["--preflight"],
                                  ["--preflight", "--run-dir", "planned", "--timeout-seconds", "nan"]])
def test_cli_options_concise(args, capsys):
    assert main(["doctor", *args]) == 1
    assert "ERROR:" in capsys.readouterr().err


def test_dry_run_no_probe(case, capsys):
    _, root, path = case
    assert main(["--config", str(path), "--dry-run", "doctor", "--preflight", "--probe-tools", "--run-dir", str(root / "planned")]) == 1
    assert "cannot be combined" in capsys.readouterr().err


@pytest.mark.parametrize("build", [None, "OTHER_BUILD"])
def test_catalog_build_must_be_declared_and_match(case, build):
    cfg, root, _ = case
    catalog = root / "repeats.bed"
    catalog.write_text("chr1\t10\t20\tID=L1;MOTIFS=A;STRUC=(A)n\n", encoding="utf-8")
    cfg = change(cfg, tr={"enabled": True, "catalog": str(catalog), "catalog_reference_build": build})
    assert check(preflight.run_preflight(cfg, root / "planned", runner=FakeRunner()), "tr_catalog")["status"] == "FAIL"


def test_reference_build_not_inferred(case):
    cfg, root, _ = case
    result = preflight.run_preflight(change(cfg, reference={"build": None}), root / "planned", runner=FakeRunner())
    assert check(result, "reference_build")["status"] == "FAIL"


def test_no_implicit_harmonization_route(case):
    cfg, root, _ = case
    values = cfg.to_dict(); values["sv"]["harmonization"]["enabled"] = True
    result = preflight.run_preflight(HiFiVarConfig(values, cfg.sources), root / "planned", runner=FakeRunner())
    assert result["status"] == "BLOCKED"


@pytest.mark.parametrize("missing", ["jasmine", "truvari", "bgzip", "tabix"])
def test_harmonization_missing_dependencies_block(case, missing):
    cfg, root, _ = case
    values = cfg.to_dict()
    values["sv"]["enabled"] = True
    values["sv"]["harmonization"]["enabled"] = True
    result = preflight.run_preflight(HiFiVarConfig(values, cfg.sources), root / "planned",
        probe_tools=True, runner=FakeRunner(missing=[missing]))
    assert result["status"] == "BLOCKED"
    assert not (root / "planned").exists()


def test_harmonization_probe_versions_only_no_analysis(case):
    cfg, root, _ = case
    values = cfg.to_dict()
    values["sv"]["enabled"] = True
    values["sv"]["harmonization"]["enabled"] = True
    before = hashes(root)
    runner = FakeRunner()
    result = preflight.run_preflight(HiFiVarConfig(values, cfg.sources), root / "planned",
        probe_tools=True, runner=runner)
    assert result["status"] == "CHECKED_SCOPE_PASS", result["checks"]
    assert check(result, "tool.jasmine.version")["details"]["version"]["jasmine"] == "1.1.5"
    assert check(result, "tool.truvari.version")["details"]["version"] == "5.4.0"
    assert check(result, "tool.jasmine.vcf44_indexer")["status"] == "PASS"
    assert hashes(root) == before
    for cmd, _ in runner.commands:
        if cmd[0] in {"jasmine", "truvari", "bgzip", "tabix"}:
            assert cmd[-1] in {"--version", "version"}


@pytest.mark.parametrize("version,state", [("1.21", "BLOCKED"), ("1.23.0", "BLOCKED"),
                                          ("1.23.1", "CHECKED_SCOPE_PASS"), ("1.24", "CHECKED_SCOPE_PASS")])
def test_harmonization_indexer_floor_before_analysis(case, version, state):
    cfg, root, _ = case
    values = cfg.to_dict()
    values["sv"]["enabled"] = True
    values["sv"]["harmonization"]["enabled"] = True

    class IndexerVersionRunner(FakeRunner):
        def run(self, command, **kwargs):
            result = super().run(command, **kwargs)
            if command == ("tabix", "--version") or command == ["tabix", "--version"]:
                result.stdout = f"tabix (htslib) {version}"
            return result

    before = hashes(root)
    result = preflight.run_preflight(HiFiVarConfig(values, cfg.sources), root / "planned",
        probe_tools=True, runner=IndexerVersionRunner())
    assert result["status"] == state
    assert check(result, "tool.jasmine.vcf44_indexer")["details"]["observed"] == version
    assert hashes(root) == before and not result["analysis_executed"]


def test_cram_is_not_promised(case):
    cfg, root, _ = case
    (root / "input.cram").write_bytes(b"synthetic")
    (root / "samples.tsv").write_text("sample_id\tinput\tinput_type\nSYNTHETIC\tinput.cram\tcram\n", encoding="utf-8")
    assert preflight.run_preflight(cfg, root / "planned", runner=FakeRunner())["status"] == "BLOCKED"


def test_actual_runner_probe_captures_without_output_files(tmp_path):
    import sys
    from hifivar.command import CommandRunner
    before = hashes(tmp_path)
    runner = preflight.ProbeRunner(CommandRunner(), 3)
    result = runner.run([sys.executable, "-I", "-B", "-c", "print('SYNTHETIC_PROBE')"])
    assert result.stdout.strip() == "SYNTHETIC_PROBE" and hashes(tmp_path) == before


def test_probe_runner_refuses_file_destinations(tmp_path):
    from hifivar.exceptions import InputValidationError
    runner = preflight.ProbeRunner(FakeRunner(), 3)
    with pytest.raises(InputValidationError):
        runner.run(["samtools", "--version"], stderr_path=tmp_path / "unexpected.log")
    assert not list(tmp_path.iterdir())


def test_catalog_out_of_reference_bounds_rejected(case):
    cfg, root, _ = case
    catalog = root / "repeats.bed"
    catalog.write_text(" # retained comment\nchr1\t999\t1001\tID=L1;MOTIFS=A;STRUC=(A)n\n", encoding="utf-8")
    cfg = change(cfg, tr={"enabled": True, "catalog": str(catalog), "catalog_reference_build": "CUSTOM_TEST"})
    assert check(preflight.run_preflight(cfg, root / "planned", runner=FakeRunner()), "tr_catalog")["status"] == "FAIL"


def test_preflight_does_not_open_configured_log(case, capsys, monkeypatch):
    _, root, path = case
    value = yaml.safe_load(path.read_text(encoding="utf-8"))
    log = root / "outside.log"
    value["logging"] = {"file": str(log)}
    path.write_text(yaml.safe_dump(value), encoding="utf-8")
    before = hashes(root)
    assert main(["--config", str(path), "doctor", "--preflight", "--run-dir", str(root / "planned"), "--json"]) == 1
    assert json.loads(capsys.readouterr().out)["status"] == "BLOCKED"
    assert not log.exists() and hashes(root) == before


def test_external_only_harmonization_keeps_own_dependency_groups(case):
    cfg, _, _ = case
    values=cfg.to_dict()
    values['small']['enabled']=False
    values['sv']['enabled']=False
    values['sv']['harmonization']['enabled']=True
    names={name for name, _, _ in preflight._tool_groups(values, object())}
    assert {'jasmine','truvari'} <= names
    assert not {'sawfish','sniffles2','pbsv','cutesv','deepvariant'} & names


@pytest.fixture
def external_harmonization(case):
    import struct,zlib
    cfg,root,path=case
    (root/'reads.fastq').write_text('@movie/1/ccs\nACGT\n+\nIIII\n',encoding='utf-8')
    (root/'samples.tsv').write_text('sample_id\tinput\tinput_type\nSYNTHETIC\treads.fastq\tfastq\n',encoding='utf-8')
    def bgzf(data):
        compressor=zlib.compressobj(wbits=-15);payload=compressor.compress(data)+compressor.flush()
        return b'\x1f\x8b\x08\x04'+b'\0'*4+b'\0\xff\x06\0BC\x02\0'+struct.pack('<H',18+len(payload)+8-1)+payload+struct.pack('<II',zlib.crc32(data),len(data))
    source=root/'SYNTHETIC.pav.assembly.sv.vcf.gz'
    text='##fileformat=VCFv4.3\n##reference=CUSTOM_TEST\n##contig=<ID=chr1,length=1000>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTHETIC\nchr1\t10\tv1\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=70;SVLEN=-60\tGT\t0/1\n'
    source.write_bytes(bgzf(text.encode()))
    # Structural index fixture only; actual tabix readability belongs to native gate.
    Path(str(source)+'.tbi').write_bytes(bgzf(b'TBI\x01'))
    values=cfg.to_dict();values['small']['enabled']=False;values['sv']['enabled']=False
    values['sv']['harmonization'].update(enabled=True,threads=1,memory_mb=100,
        input_vcfs={'SYNTHETIC':{'pav':str(source)}})
    return HiFiVarConfig(values,cfg.sources),root


@pytest.mark.parametrize('missing',['jasmine','truvari','bgzip','tabix'])
def test_external_only_missing_dependency_is_blocked_without_mutation(external_harmonization,missing):
    cfg,root=external_harmonization;before=hashes(root)
    result=preflight.run_preflight(cfg,root/'planned',probe_tools=True,runner=FakeRunner(missing=[missing]))
    assert result['status']=='BLOCKED'
    assert hashes(root)==before and not result['analysis_executed'] and not (root/'planned').exists()


def test_external_only_probes_its_tools_not_disabled_callers_or_unused_bam_tool(external_harmonization):
    cfg,root=external_harmonization;before=hashes(root)
    runner=FakeRunner(missing=['samtools','run_deepvariant','sawfish','sniffles','pbsv','cuteSV'])
    result=preflight.run_preflight(cfg,root/'planned',probe_tools=True,runner=runner)
    assert result['status']=='CHECKED_SCOPE_PASS',result['checks']
    assert check(result,'tool.jasmine.version')['status']=='PASS'
    assert check(result,'tool.truvari.version')['status']=='PASS'
    assert not any(cmd[0] in {'samtools','run_deepvariant','sawfish','sniffles','pbsv','cuteSV'} for cmd,_ in runner.commands)
    assert hashes(root)==before and not result['ready_for_all_workflows']


@pytest.mark.parametrize('version,state',[('1.21','BLOCKED'),('1.23.1','CHECKED_SCOPE_PASS')])
def test_external_only_indexer_compatibility_floor(external_harmonization,version,state):
    cfg,root=external_harmonization
    class VersionRunner(FakeRunner):
        def run(self,command,**kwargs):
            result=super().run(command,**kwargs)
            if list(command)==['tabix','--version']:result.stdout='tabix (htslib) '+version
            return result
    result=preflight.run_preflight(cfg,root/'planned',probe_tools=True,runner=VersionRunner())
    assert result['status']==state


def test_external_only_metadata_probe_remains_incomplete_without_version_execution(external_harmonization):
    cfg,root=external_harmonization;runner=FakeRunner()
    result=preflight.run_preflight(cfg,root/'planned',runner=runner)
    assert result['status']=='INCOMPLETE' and not runner.commands
    assert check(result,'tool.jasmine.presence')['status']=='PASS'
    assert check(result,'tool.truvari.version')['status']=='NOT_CHECKED'


@pytest.mark.parametrize('version,state',[('1.21','BLOCKED'),('1.23.1','CHECKED_SCOPE_PASS')])
def test_identity_bgzip_floor_matches_runtime_wrapper(external_harmonization,version,state):
    from hifivar.sv_input_identity import POLICY,POLICY_SHA256
    cfg,root=external_harmonization;values=cfg.to_dict()
    values['sv']['harmonization']['input_identity'].update(enabled=True,
        policy_id=POLICY['policy_id'],policy_version=POLICY['policy_version'],
        policy_sha256=POLICY_SHA256,mode=POLICY['mode'])
    class VersionRunner(FakeRunner):
        def run(self,command,**kwargs):
            result=super().run(command,**kwargs)
            if list(command)==['bgzip','--version']:result.stdout='bgzip (htslib) '+version
            return result
    result=preflight.run_preflight(HiFiVarConfig(values,cfg.sources),root/'planned',probe_tools=True,runner=VersionRunner())
    assert result['status']==state
    assert check(result,'tool.jasmine.source_identity_compressor')['details']['observed']==version


def test_enabled_hifiasm_uses_its_own_version_probe(external_harmonization):
    cfg,root=external_harmonization;values=cfg.to_dict()
    values['assembly'].update(enabled=True,threads=1,memory_mb=100)
    result=preflight.run_preflight(HiFiVarConfig(values,cfg.sources),root/'planned',probe_tools=True,runner=FakeRunner())
    assert result['status']=='CHECKED_SCOPE_PASS',result['checks']
    item=check(result,'tool.assembly.hifiasm.version')
    assert item['status']=='PASS' and item['details']['version']=='0.25.0-r726'
    assert item['details']['version_source']=='command'
    assert not result['analysis_executed']


# Probe-registry isolation tests deliberately keep launcher/input preparation
# mocked: they do not claim real assembly/cache/gVCF analysis acceptance.
@pytest.fixture
def extended_native_registry(case, monkeypatch):
    import logging
    cfg, root, _ = case
    # Earlier CLI cases own capture-stream handlers. Do not reuse their closed
    # pytest stream while constructing this independent registry-only fixture.
    monkeypatch.setattr(logging.getLogger("hifivar"), "handlers", [])
    effective, context, manifest, inputs = preflight.prepare_launch(cfg, root / 'planned')
    effective['assembly']['enabled'] = True
    effective['assembly_sv']['enabled'] = True
    effective['assembly_sv']['pav']['enabled'] = False
    effective['annotation'].update(enabled=True, vep_enabled=True, annovar_enabled=False)
    effective['cohort']['enabled'] = True
    effective['cohort']['small_variants']['enabled'] = True
    monkeypatch.setattr(preflight, 'prepare_launch', lambda *args, **kwargs: (effective, context, manifest, inputs))
    return cfg, root, effective


def test_extended_native_probes_keep_wrapper_ownership_and_no_outputs(extended_native_registry, monkeypatch):
    cfg, root, _ = extended_native_registry
    before = hashes(root); runner = FakeRunner()
    monkeypatch.setenv('EXTENDED_PROBE_PRIVATE_TOKEN', 'not-to-be-inherited')
    result = preflight.run_preflight(cfg, root / 'planned', probe_tools=True, runner=runner, timeout_seconds=5)
    assert result['status'] == 'CHECKED_SCOPE_PASS', result['checks']
    for name in ['assembly.hifiasm', 'assembly_sv.svim_asm', 'annotation.vep', 'cohort.glnexus']:
        item = check(result, 'tool.' + name + '.version')
        assert item['status'] == 'PASS' and item['details']['version_source'] == 'command'
        assert item['details']['validation_scope'] == 'WRAPPER_NATIVE_VERSION_ONLY_NOT_ANALYSIS_OR_DATABASE_CONTENT'
    native = [tuple(cmd) for cmd, _ in runner.commands if cmd[0] in {'hifiasm','svim-asm','vep','glnexus_cli'}]
    assert native == [('hifiasm','--version'),('svim-asm','--version'),('vep','--help'),('glnexus_cli','--help')]
    for command, kwargs in runner.commands:
        assert kwargs['timeout'] == 5 and kwargs['inherit_env'] is False
        assert 'EXTENDED_PROBE_PRIVATE_TOKEN' not in kwargs['env']
        assert not kwargs.get('stdout_path') and not kwargs.get('stderr_path')
    assert hashes(root) == before and not (root / 'planned').exists()
    assert not result['analysis_executed'] and not result['ready_for_all_workflows']


@pytest.mark.parametrize('tool,group', [('hifiasm','assembly.hifiasm'),('svim-asm','assembly_sv.svim_asm'),
    ('minimap2','assembly_sv.svim_asm'),('samtools','assembly_sv.svim_asm'),('bgzip','assembly_sv.svim_asm'),
    ('tabix','assembly_sv.svim_asm'),('vep','annotation.vep'),('glnexus_cli','cohort.glnexus'),('bcftools','cohort.glnexus')])
def test_extended_missing_dependency_blocks_own_group(extended_native_registry, tool, group):
    cfg, root, _ = extended_native_registry
    before = hashes(root); runner = FakeRunner(missing=[tool])
    result = preflight.run_preflight(cfg, root / 'planned', probe_tools=True, runner=runner)
    assert result['status'] == 'BLOCKED'
    assert check(result, 'tool.' + group + '.presence')['status'] == 'FAIL'
    assert not any(tuple(cmd) == (tool, '--version') or tuple(cmd) == (tool, '--help') for cmd, _ in runner.commands)
    assert not any(x['key'] == 'tool.' + group + '.version' for x in result['checks'])
    assert hashes(root) == before


@pytest.mark.parametrize('tool,group', [('hifiasm','assembly.hifiasm'),('svim-asm','assembly_sv.svim_asm'),
    ('vep','annotation.vep'),('glnexus_cli','cohort.glnexus')])
def test_extended_version_failure_is_scoped(extended_native_registry, tool, group):
    cfg, root, _ = extended_native_registry; before = hashes(root)
    result = preflight.run_preflight(cfg, root / 'planned', probe_tools=True, runner=FakeRunner(fail_version=tool))
    assert result['status'] == 'BLOCKED'
    assert check(result, 'tool.' + group + '.version')['status'] == 'FAIL'
    assert check(result, 'input.SYNTHETIC')['status'] == 'PASS'
    assert hashes(root) == before


def test_extended_metadata_does_not_execute_or_borrow_versions(extended_native_registry):
    cfg, root, _ = extended_native_registry; before = hashes(root); runner = FakeRunner()
    result = preflight.run_preflight(cfg, root / 'planned', runner=runner)
    assert result['status'] == 'INCOMPLETE' and not runner.commands
    for name in ['assembly.hifiasm','assembly_sv.svim_asm','annotation.vep','cohort.glnexus']:
        assert check(result, 'tool.' + name + '.version')['status'] == 'NOT_CHECKED'
    assert hashes(root) == before


@pytest.mark.parametrize('section', ['assembly','assembly_sv','annotation','cohort'])
def test_extended_disabled_branches_do_not_require_tools(extended_native_registry, section):
    _, _, effective = extended_native_registry
    effective[section]['enabled'] = False
    names = {name for name, _, _ in preflight._tool_groups(effective, object())}
    prefixes={'assembly':'assembly.hifiasm','assembly_sv':'assembly_sv.svim_asm','annotation':'annotation.vep','cohort':'cohort.glnexus'}
    assert prefixes[section] not in names


@pytest.mark.parametrize('section,child', [('assembly_sv','svim_asm'),('annotation','vep_enabled'),('cohort','small_variants')])
def test_extended_disabled_child_does_not_probe(extended_native_registry, section, child):
    _, _, effective = extended_native_registry
    if section == 'annotation':effective[section][child] = False
    else:effective[section][child]['enabled'] = False
    names = {name for name, _, _ in preflight._tool_groups(effective, object())}
    assert {'assembly_sv':'assembly_sv.svim_asm','annotation':'annotation.vep','cohort':'cohort.glnexus'}[section] not in names


def test_pav_and_annovar_declarations_do_not_become_native_versions(extended_native_registry):
    cfg, root, effective = extended_native_registry
    effective['assembly_sv']['pav']['enabled'] = True
    effective['assembly_sv']['pav']['version'] = '2.4.6_DECLARED_ONLY'
    effective['annotation']['annovar_enabled'] = True
    effective['annotation']['annovar_version'] = 'DECLARED_ONLY'
    runner = FakeRunner(); before = hashes(root)
    result = preflight.run_preflight(cfg, root / 'planned', probe_tools=True, runner=runner)
    assert result['status'] == 'INCOMPLETE'
    for name in ['assembly_sv.pav','annotation.annovar']:
        item=check(result,'tool.'+name+'.probe_scope')
        assert item['status']=='NOT_CHECKED' and item['required']
        assert not any(x['key']=='tool.'+name+'.version' for x in result['checks'])
    assert not any(cmd[0] in {'snakemake','table_annovar.pl','apptainer','docker'} for cmd,_ in runner.commands)
    assert hashes(root)==before


def test_extended_custom_paths_are_used_without_default_fallback(extended_native_registry):
    _, _, effective = extended_native_registry
    effective['assembly']['executable'] = 'owned/hifiasm'
    effective['assembly_sv']['svim_asm']['executable'] = 'owned/svim-asm'
    effective['annotation']['vep_executable'] = 'owned/vep'
    effective['cohort']['small_variants']['glnexus_executable'] = 'owned/glnexus_cli'
    groups={name:deps for name,deps,_ in preflight._tool_groups(effective,object())}
    for group,tool in [('assembly.hifiasm','hifiasm'),('assembly_sv.svim_asm','svim_asm'),('annotation.vep','vep'),('cohort.glnexus','glnexus')]:
        assert groups[group][tool].startswith('owned/')
