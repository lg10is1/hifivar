from __future__ import annotations
import os
import subprocess
import sys
import tarfile
import venv
import zipfile
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]


def test_wheel_contains_runtime_resources_and_clean_install_smoke(tmp_path: Path) -> None:
    subprocess_env={key:value for key,value in os.environ.items() if key not in {"PYTHONPATH","PYTHONHOME"}}
    wheelhouse=tmp_path/"wheelhouse"
    build=subprocess.run([sys.executable,"-m","build","--no-isolation","--outdir",str(wheelhouse)],cwd=ROOT,capture_output=True,text=True,encoding="utf-8",errors="replace",env=subprocess_env)
    assert build.returncode==0,build.stdout+build.stderr
    wheel=next(wheelhouse.glob("hifivar-*.whl"))
    with zipfile.ZipFile(wheel) as archive:
        names=archive.namelist(); payload=b"".join(archive.read(name) for name in names)
    assert any(name.endswith("hifivar/workflow/Snakefile") for name in names)
    assert any(name.endswith("hifivar/workflow/rules/benchmark.smk") for name in names)
    assert any(name.endswith("hifivar/resources/configs/default.yaml") for name in names)
    assert b"C:\\Users\\private-user" not in payload
    assert not any("internal_audit" in name or name.endswith((".sif",".pem",".key")) for name in names)
    private_key_marker = b"BEGIN " + b"PRIVATE KEY"
    assert private_key_marker not in payload and b"/home/private-user" not in payload
    sdist=next(wheelhouse.glob("hifivar-*.tar.gz"))
    with tarfile.open(sdist,"r:gz") as archive: source_names=archive.getnames()
    assert any(name.endswith("docs/phases/phase14.md") for name in source_names)
    assert any(name.endswith("RELEASE_CHECKLIST.md") for name in source_names)
    assert any(name.endswith("environment.yml") for name in source_names)
    assert any(name.endswith("conda-recipe/meta.yaml") for name in source_names)
    assert any(name.endswith("containers/hifivar.def") for name in source_names)
    assert any(name.endswith("examples/minimal/config.yaml") for name in source_names)
    assert any(name.endswith("CITATION.cff") for name in source_names)
    assert not any("internal_audit" in name for name in source_names)
    environment=tmp_path/"clean-env"; venv.EnvBuilder(with_pip=True).create(environment)
    python=environment/("Scripts/python.exe" if os.name=="nt" else "bin/python")
    install=subprocess.run([str(python),"-m","pip","install",str(wheel)],cwd=tmp_path,capture_output=True,text=True,encoding="utf-8",errors="replace",env=subprocess_env)
    assert install.returncode==0,install.stdout+install.stderr
    smoke="from importlib.metadata import version; from importlib.resources import files; import hifivar; from hifivar.package_resources import installed_workflow_root; assert version('hifivar') == hifivar.__version__; assert files('hifivar').joinpath('resources/configs/default.yaml').is_file(); assert installed_workflow_root().joinpath('rules/benchmark.smk').is_file(); print(hifivar.__version__)"
    result=subprocess.run([str(python),"-c",smoke],cwd=tmp_path,capture_output=True,text=True,encoding="utf-8",errors="replace",env=subprocess_env)
    assert result.returncode==0,result.stdout+result.stderr
    executable=environment/("Scripts/hifivar.exe" if os.name=="nt" else "bin/hifivar")
    for args in (("--help",),("--version",),("config","validate")):
        cli=subprocess.run([str(executable),*args],cwd=tmp_path,capture_output=True,text=True,encoding="utf-8",errors="replace",env=subprocess_env)
        assert cli.returncode==0,cli.stdout+cli.stderr
