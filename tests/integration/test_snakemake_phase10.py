from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

from hifivar.config import load_config, write_effective_config


ROOT = Path(__file__).resolve().parents[2]
SNAKEFILE = ROOT / "workflow" / "Snakefile"
DEFAULT = ROOT / "src" / "hifivar" / "resources" / "configs" / "default.yaml"


def snakemake() -> Path:
    name = "snakemake.exe" if sys.platform == "win32" else "snakemake"
    adjacent = Path(sys.executable).with_name(name)
    found = adjacent if adjacent.is_file() else shutil.which("snakemake")
    if found is None:
        pytest.skip("Snakemake is not installed.")
    return Path(found)


def test_phase10_review_is_optional_independent_dry_run(tmp_path):
    reference = tmp_path / "reference.fa"
    reference.write_text(">chr1\nACGT\n", encoding="utf-8")
    Path(f"{reference}.fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")
    bam = tmp_path / "S1.bam"
    bam.write_bytes(b"BAM")
    Path(f"{bam}.bai").write_bytes(b"BAI")
    samples = tmp_path / "samples.tsv"
    samples.write_text(f"sample_id\tinput\tinput_type\nS1\t{bam}\tbam\n", encoding="utf-8")
    source = tmp_path / "explicit.vcf"
    source.write_text(
        "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "chr1\t2\tv1\tA\tT\t.\tPASS\t.\tGT\t0/1\n",
        encoding="utf-8",
    )
    selection = tmp_path / "selection.tsv"
    selection.write_text(
        "review_id\tsample\tvariant_id\tvariant_type\tcontig\tstart\tend\t"
        "source_vcf\tsource_caller\tevidence_class\n"
        f"R1\tS1\tv1\tSNV\tchr1\t2\t2\t{source}\texplicit\texplicit\n",
        encoding="utf-8",
    )
    user = tmp_path / "phase10.yaml"
    user.write_text(
        yaml.safe_dump(
            {
                "reference": {"fasta": str(reference), "build": "GRCh38"},
                "samples": {"sheet": str(samples)},
                "paths": {"workdir": str(tmp_path / "work"), "outdir": str(tmp_path / "results")},
                "review": {"enabled": True, "selection_file": str(selection)},
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )
    effective = tmp_path / "effective.yaml"
    write_effective_config(load_config(DEFAULT, user_config=user), effective)
    environment = os.environ.copy()
    environment["XDG_CACHE_HOME"] = str(tmp_path / ".cache")
    if sys.platform == "win32":
        environment["LOCALAPPDATA"] = str(tmp_path / ".local")
        environment["APPDATA"] = str(tmp_path / ".app")
    result = subprocess.run(
        [
            str(snakemake()), "--snakefile", str(SNAKEFILE), "--configfile", str(effective),
            "--cores", "1", "--dry-run", "--printshellcmds", "--shared-fs-usage",
            "input-output", "persistence", "software-deployment",
            "software-deployment-cache", "sources", "storage-local-copies",
        ],
        cwd=tmp_path,
        env=environment,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        check=False,
        timeout=60,
    )
    output = result.stdout + result.stderr
    assert result.returncode == 0, output
    assert "manual_review" in output
    assert "review_manifest.json" in output and "review.igv.batch" in output
    assert "deepvariant_small" not in output and "harmonize_sv" not in output


def test_phase10_bridge_uses_wrappers_not_subprocess():
    bridge = (ROOT / "workflow" / "scripts" / "run_review.py").read_text(encoding="utf-8")
    assert "run_review_step" in bridge
    orchestration = (ROOT / "src/hifivar/review_handoff.py").read_text(encoding="utf-8")
    assert "IgvWrapper" in orchestration and "run_phase10" in orchestration and "attempt_step" in orchestration
    assert "subprocess" not in orchestration and "os.system" not in orchestration
    assert "subprocess" not in bridge and "os.system" not in bridge



@pytest.mark.skipif(sys.platform=='win32',reason='POSIX synthetic executable runtime; Windows has separate fake-runner contracts')
def test_phase10_standalone_snakemake_fake_executable_runtime(tmp_path):
    """Real scheduler/script boundary, synthetic BAM/VCF/IGV; no biological claim."""
    reference=tmp_path/'ref.fa';reference.write_text('>chr1\nACGT\n',encoding='utf-8')
    Path(str(reference)+'.fai').write_text('chr1\t4\t6\t4\t5\n',encoding='utf-8')
    bam=tmp_path/'S1.bam';bam.write_bytes(b'synthetic BAM');Path(str(bam)+'.bai').write_bytes(b'synthetic BAI')
    source=tmp_path/'source.vcf';source.write_text('##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\nchr1\t2\tv1\tC\tT\t.\tPASS\t.\tGT\t0/1\n')
    sheet=tmp_path/'samples.tsv';sheet.write_text(f'sample_id\tinput\tinput_type\nS1\t{bam}\tbam\n')
    selection=tmp_path/'selection.tsv';selection.write_text('review_id\tsample\tvariant_id\tvariant_type\tcontig\tstart\tend\tsource_vcf\tsource_caller\tevidence_class\n'+f'R1\tS1\tv1\tSNV\tchr1\t2\t2\t{source}\tfixture\texplicit\n')
    igv=tmp_path/'fixture_igv.py'
    program=r'''import sys,shlex,struct,zlib
from pathlib import Path
if '--version' in sys.argv:
 print('IGV 9.9.0 synthetic fixture');raise SystemExit(0)
def chunk(kind,data):return struct.pack('>I',len(data))+kind+data+struct.pack('>I',zlib.crc32(kind+data))
png=b'\x89PNG\r\n\x1a\n'+chunk(b'IHDR',struct.pack('>IIBBBBB',64,64,8,2,0,0,0))+chunk(b'IDAT',zlib.compress((b'\x00'+b'\x88\xaa\xcc'*64)*64))+chunk(b'IEND',b'')
directory=None
for raw in Path(sys.argv[-1]).read_text().splitlines():
 fields=shlex.split(raw)
 if fields and fields[0]=='snapshotDirectory':directory=Path(fields[1])
 if fields and fields[0]=='snapshot':
  directory.mkdir(parents=True,exist_ok=True);(directory/fields[1]).write_bytes(png)
'''
    igv.write_text('#!'+sys.executable+'\n'+program,encoding='utf-8');igv.chmod(0o755)
    user=tmp_path/'user.yaml';user.write_text(yaml.safe_dump({'reference':{'fasta':str(reference),'build':'GRCh38'},
        'samples':{'sheet':str(sheet)},'paths':{'workdir':str(tmp_path/'work'),'outdir':str(tmp_path/'results')},
        'review':{'enabled':True,'selection_file':str(selection),'igv_executable':str(igv)}}))
    effective=tmp_path/'effective.yaml';write_effective_config(load_config(DEFAULT,user_config=user),effective)
    env={'PATH':str(Path(sys.executable).parent)+os.pathsep+os.defpath,'HOME':str(tmp_path),
         'TMPDIR':str(tmp_path),'XDG_CACHE_HOME':str(tmp_path/'.cache'),'XDG_CONFIG_HOME':str(tmp_path/'.config'),
         'PYTHONPATH':str(ROOT/'src'),'PYTHONDONTWRITEBYTECODE':'1','LANG':'C','LC_ALL':'C'}
    result=subprocess.run([str(snakemake()),'--snakefile',str(SNAKEFILE),'--configfile',str(effective),'--cores','1',
        '--shared-fs-usage','input-output','persistence','software-deployment','software-deployment-cache','sources','storage-local-copies','--','all'],
        cwd=tmp_path,env=env,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,timeout=120,check=False)
    assert result.returncode==0,result.stdout+result.stderr
    output=tmp_path/'results/review'
    assert (output/'review_provenance.json').is_file()
    assert (output/'screenshots/R1.locus01.png').read_bytes().startswith(b'\x89PNG')
    import json
    provenance=json.loads((output/'review_provenance.json').read_text())
    assert provenance['igv']['tool_version']=='9.9.0'
    assert provenance['manifest']['results'][0]['status']=='NOT_REVIEWED'
    assert source.read_text().endswith('chr1\t2\tv1\tC\tT\t.\tPASS\t.\tGT\t0/1\n')
