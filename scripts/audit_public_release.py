"""Audit a HiFiVar source tree for public-release hygiene.

This script performs read-only checks and exits non-zero when it finds internal
workspace references, private path fragments, private-key material, generated
artifacts, or unexpectedly large files.
"""

from __future__ import annotations

import argparse
import re
import tarfile
import zipfile
from email.parser import BytesParser
from pathlib import Path


TEXT_PATTERNS = {
    "internal assistant or audit reference": re.compile(
        r"(?:Kimi|ChatGPT|OpenCode|Codex|00\.hifivar_kimi|01\.hifivar_opencode)",
        re.IGNORECASE,
    ),
    "private identity or path": re.compile(
        r"(?:yangqiangzhen1|809078831|qq\.com|F:\\科研项目|"
        r"C:\\Users\\My|/data/yangqiangzhen1|/home/yangqiangzhen1)",
        re.IGNORECASE,
    ),
    "private key material": re.compile(
        r"BEGIN (?:RSA |OPENSSH |EC |DSA )?PRIVATE KEY",
        re.IGNORECASE,
    ),
    "email address": re.compile(
        r"\b[A-Z0-9._%+-]+@[A-Z0-9.-]+\.[A-Z]{2,}\b",
        re.IGNORECASE,
    ),
    "credential-looking token": re.compile(
        r"\b(?:AKIA[0-9A-Z]{16}|gh[opsu]_[A-Za-z0-9]{30,})\b"
    ),
    "credential embedded in URL": re.compile(
        r"[a-z][a-z0-9+.-]*://[^\s/:]+:[^\s/@]+@",
        re.IGNORECASE,
    ),
}

NAME_PATTERN = re.compile(
    r"(?:kimi|chatgpt|opencode|codex|id_rsa|id_ed25519|private\.key)",
    re.IGNORECASE,
)
GENERATED_PARTS = {".git", ".pytest_cache", ".venv", "build", "dist", "__pycache__"}
TEXT_SUFFIXES = {
    "",
    ".cfg",
    ".ini",
    ".json",
    ".md",
    ".py",
    ".smk",
    ".toml",
    ".tsv",
    ".txt",
    ".yaml",
    ".yml",
}


def audit(root: Path, max_bytes: int) -> list[str]:
    """Return public-release hygiene findings for *root*."""

    findings: list[str] = []
    audit_script = Path(__file__).resolve()
    for path in sorted(root.rglob("*")):
        relative = path.relative_to(root)
        if NAME_PATTERN.search(relative.as_posix()):
            findings.append(f"suspicious filename: {relative}")
        if any(part in GENERATED_PARTS for part in relative.parts):
            continue
        if not path.is_file():
            continue
        size = path.stat().st_size
        if size > max_bytes:
            findings.append(f"large file ({size} bytes): {relative}")
        if path.resolve() == audit_script or path.suffix.lower() not in TEXT_SUFFIXES:
            continue
        text = path.read_text(encoding="utf-8", errors="ignore")
        for label, pattern in TEXT_PATTERNS.items():
            for line_number, line in enumerate(text.splitlines(), start=1):
                if pattern.search(line):
                    findings.append(f"{label}: {relative}:{line_number}")
    return findings


def audit_payload(label: str, name: str, payload: bytes) -> list[str]:
    """Return findings for one archive member without extracting it."""

    findings: list[str] = []
    if NAME_PATTERN.search(name):
        findings.append(f"suspicious archive member: {label}:{name}")
    text = payload.decode("utf-8", errors="ignore")
    for finding_label, pattern in TEXT_PATTERNS.items():
        if pattern.search(text):
            findings.append(f"{finding_label} in archive: {label}:{name}")
    return findings


def audit_artifacts(artifacts_dir: Path) -> list[str]:
    """Audit wheel and source-distribution members in *artifacts_dir*."""

    findings: list[str] = []
    wheels = sorted(artifacts_dir.glob("hifivar-*.whl"))
    sdists = sorted(artifacts_dir.glob("hifivar-*.tar.gz"))
    if len(wheels) != 1 or len(sdists) != 1:
        return [
            "expected exactly one HiFiVar wheel and one source distribution in "
            f"{artifacts_dir}"
        ]
    wheel_version = ""
    with zipfile.ZipFile(wheels[0]) as archive:
        metadata_names = [name for name in archive.namelist() if name.endswith(".dist-info/METADATA")]
        if len(metadata_names) != 1:
            findings.append(f"expected exactly one wheel METADATA member: {wheels[0].name}")
        else:
            metadata = BytesParser().parsebytes(archive.read(metadata_names[0]))
            if metadata.get("Name") != "hifivar":
                findings.append(f"unexpected wheel project name: {metadata.get('Name')!r}")
            wheel_version = metadata.get("Version", "")
        for name in archive.namelist():
            findings.extend(audit_payload(wheels[0].name, name, archive.read(name)))
    with tarfile.open(sdists[0], "r:gz") as archive:
        pkg_info_names = [
            member for member in archive.getmembers() if member.name.count("/") == 1 and member.name.endswith("/PKG-INFO")
        ]
        if len(pkg_info_names) != 1:
            findings.append(f"expected exactly one top-level sdist PKG-INFO: {sdists[0].name}")
        else:
            stream = archive.extractfile(pkg_info_names[0])
            assert stream is not None
            metadata = BytesParser().parsebytes(stream.read())
            if metadata.get("Name") != "hifivar":
                findings.append(f"unexpected sdist project name: {metadata.get('Name')!r}")
            if metadata.get("Version") != wheel_version:
                findings.append("wheel and sdist versions do not match")
        for member in archive.getmembers():
            if member.isfile():
                stream = archive.extractfile(member)
                assert stream is not None
                findings.extend(audit_payload(sdists[0].name, member.name, stream.read()))
    if wheel_version:
        expected_wheel = f"hifivar-{wheel_version}-py3-none-any.whl"
        expected_sdist = f"hifivar-{wheel_version}.tar.gz"
        if wheels[0].name != expected_wheel:
            findings.append(f"unexpected wheel filename: {wheels[0].name}")
        if sdists[0].name != expected_sdist:
            findings.append(f"unexpected sdist filename: {sdists[0].name}")
    return findings


def compare_production_trees(root: Path, reference_root: Path) -> list[str]:
    """Confirm that production Python and workflow files match a reference tree."""

    findings: list[str] = []

    def comparable(path: Path) -> bool:
        return (
            path.is_file()
            and path.suffix.lower() not in {".pyc", ".pyo"}
            and not any(
                part in GENERATED_PARTS or part.endswith(".egg-info")
                for part in path.parts
            )
        )

    def normalized_content(path: Path) -> bytes:
        if path.suffix.lower() in TEXT_SUFFIXES or path.name == "Snakefile":
            return "\n".join(
                path.read_text(encoding="utf-8", errors="strict").splitlines()
            ).encode("utf-8")
        return path.read_bytes()

    for relative_root in (Path("src/hifivar"), Path("workflow")):
        public_files = {
            path.relative_to(root): path
            for path in (root / relative_root).rglob("*")
            if comparable(path)
        }
        reference_files = {
            path.relative_to(reference_root): path
            for path in (reference_root / relative_root).rglob("*")
            if comparable(path)
        }
        for relative in sorted(public_files.keys() | reference_files.keys()):
            if relative not in public_files:
                findings.append(f"production file missing from public tree: {relative}")
            elif relative not in reference_files:
                findings.append(f"unexpected production file in public tree: {relative}")
            elif normalized_content(public_files[relative]) != normalized_content(reference_files[relative]):
                findings.append(f"production file differs from reference tree: {relative}")
    return findings


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", nargs="?", type=Path, default=Path.cwd())
    parser.add_argument("--max-mib", type=int, default=10)
    parser.add_argument("--artifacts-dir", type=Path)
    parser.add_argument("--reference-root", type=Path)
    args = parser.parse_args()

    root = args.root.resolve()
    findings = audit(root, args.max_mib * 1024 * 1024)
    if args.artifacts_dir is not None:
        findings.extend(audit_artifacts(args.artifacts_dir.resolve()))
    if args.reference_root is not None:
        findings.extend(compare_production_trees(root, args.reference_root.resolve()))
    if findings:
        print("PUBLIC RELEASE AUDIT: FAIL")
        for finding in findings:
            print(f"- {finding}")
        return 1
    print("PUBLIC RELEASE AUDIT: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
