from __future__ import annotations
import os, shutil, subprocess, sys
from pathlib import Path
import pytest, yaml
from hifivar.config import load_config, write_effective_config

ROOT=Path(__file__).resolve().parents[2]; SNAKEFILE=ROOT/"workflow"/"Snakefile"; DEFAULT=ROOT/"src"/"hifivar"/"resources"/"configs"/"default.yaml"

def _snakemake():
    name="snakemake.exe" if sys.platform=="win32" else "snakemake"; adjacent=Path(sys.executable).with_name(name); found=adjacent if adjacent.is_file() else shutil.which("snakemake")
    if found is None: pytest.skip("Snakemake is not installed.")
    return Path(found)

def test_phase13_small_track_is_optional_independent_dry_run(tmp_path: Path) -> None:
    ref=tmp_path/"ref.fa"; ref.write_text(">chr1\nACGT\n"); Path(f"{ref}.fai").write_text("chr1\t4\t6\t4\t5\n")
    bam=tmp_path/"S1.bam"; bam.write_bytes(b"BAM"); samples=tmp_path/"samples.tsv"; samples.write_text(f"sample_id\tinput\tinput_type\nS1\t{bam}\tbam\n")
    query=tmp_path/"q.vcf.gz"; truth=tmp_path/"t.vcf.gz"; bed=tmp_path/"confident.bed"
    for path in (query,truth): path.write_bytes(b"vcf"); Path(f"{path}.tbi").write_bytes(b"idx")
    bed.write_text("chr1\t0\t4\n")
    user=tmp_path/"config.yaml"; user.write_text(yaml.safe_dump({"reference":{"fasta":str(ref),"build":"GRCh38"},"samples":{"sheet":str(samples)},"paths":{"workdir":str(tmp_path/"work"),"outdir":str(tmp_path/"results")},"benchmark":{"enabled":True,"benchmark_id":"B1","sample_id":"S1","small_variants":{"enabled":True,"query_vcf":str(query),"truth_vcf":str(truth),"truth_version":"v1","truth_source":"synthetic","confident_bed":str(bed),"confident_bed_version":"v1"}}},sort_keys=False))
    effective=tmp_path/"effective.yaml"; write_effective_config(load_config(DEFAULT,user_config=user),effective)
    env=os.environ.copy(); env["XDG_CACHE_HOME"]=str(tmp_path/"cache"); env["LOCALAPPDATA"]=str(tmp_path/"local"); env["APPDATA"]=str(tmp_path/"app")
    result=subprocess.run([str(_snakemake()),"--snakefile",str(SNAKEFILE),"--configfile",str(effective),"--cores","1","--dry-run","--shared-fs-usage","input-output","persistence","software-deployment","software-deployment-cache","sources","storage-local-copies"],cwd=tmp_path,env=env,capture_output=True,text=True,encoding="utf-8",errors="replace",timeout=60)
    output=result.stdout+result.stderr; assert result.returncode==0,output
    assert "benchmark_small_variants" in output and "benchmark_manifest" in output
    assert "happy.metrics.json.gz" in output
    assert "benchmark_sv" not in output and "benchmark_tr" not in output

def test_phase13_bridges_never_invoke_subprocess_directly() -> None:
    for path in (ROOT/"workflow"/"scripts").glob("run_benchmark_*.py"):
        text=path.read_text(encoding="utf-8"); assert "subprocess" not in text and "os.system" not in text
