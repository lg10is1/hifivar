# HiFiVar release checklist

[简体中文](RELEASE_CHECKLIST.zh-CN.md)

No item in this checklist authorizes an automatic push, tag, PyPI/GitHub
release, container publication, or version bump. A human release owner must
approve those actions explicitly.

## 0.1.0rc5 source/package candidate identity

- Canonical version source: `src/hifivar/__init__.py`
- Expected wheel basename: `hifivar-0.1.0rc5-py3-none-any.whl`
- Expected sdist basename: `hifivar-0.1.0rc5.tar.gz`
- Build in a distinct dated directory below `dist/rc5/`; never overwrite an
  earlier candidate while retaining its identity.
- Preserve all RC1-RC4 tags, releases and artifacts.
- Do not equate a working-tree source manifest with a release commit.
- Candidate preparation does not create a commit/tag or publish an artifact.
- Multiple development builds may share `0.1.0rc5`; record the chosen candidate's
  dated artifact directory, source manifest and exact wheel/sdist SHA256.
- Confirm the chosen package contains the documented launcher before advertising
  `hifivar run`; the initial RC5 security-only build does not contain it.
- The [15 October payload](docs/rc5_candidate.md#validated-development-payload--已验收开发包)
  has dedicated Linux/package and bounded launcher acceptance. Its hashes describe
  that historical payload; later documentation revisions require new measured
  artifact hashes and matching source manifests.

## Source and scope

- [x] Apache License 2.0 selected for HiFiVar source and package metadata.
- [x] Confirm public repository URLs point to `https://github.com/lg10is1/hifivar`.
- [ ] Confirm intended release scope and candidate version.
- [ ] Confirm `git status --short` is clean or every remaining path is reviewed.
- [ ] Review `git diff --check` and the complete staged diff.
- [ ] Confirm no new caller or unreviewed scientific policy entered the release.
- [ ] Confirm CHANGELOG and phase documentation match actual behavior.

## Tests

- [ ] `python -m pytest -p no:cacheprovider` passes.
- [ ] `python -m compileall src` passes.
- [ ] `python -m build` produces both wheel and sdist.
- [ ] All Snakemake dry-run/mock regressions pass.
- [ ] CLI help, version, config validation, errors and Unicode paths pass.
- [ ] Wheel clean-install smoke passes outside the repository.
- [ ] Packaged default/preset YAML and installed workflow resource checks pass.
- [ ] `environment.yml` creates successfully in a clean Linux environment.
- [ ] Local `conda-recipe` build/install/import/CLI/resource tests pass.
- [ ] Docker core image build and non-root smoke tests pass, or remain explicitly NOT_RUN.
- [ ] Apptainer core build/test/run smoke passes, or remains explicitly NOT_RUN.

## Linux/HPC evidence

- [ ] Review every external-tool status in `docs/deployment.md`.
- [ ] Required production tools have exact versions/backends and current real-run evidence.
- [ ] Environment-blocked or NOT_RUN tools remain explicitly labelled and are
      either disabled for the candidate or accepted as documented limitations.
- [ ] Validate tiny real-tool workflows before chromosome/WGS execution.
- [ ] Record effective config, Git SHA, commands, logs and reference/tool/database identities.

## Artifacts and security

- [ ] Complete the maintainer's credential and secret-management review;
      environment isolation and file permissions do not invalidate copied keys.
- [ ] Verify both successful and failing hap.py calls exclude unrelated secrets
      and protect runinfo; retain sensitive originals outside public artifacts.

- [ ] Inspect wheel and sdist file lists.
- [ ] Confirm no audit scratch, raw WGS data, BAM/CRAM/FASTQ, SIF, caches or logs
      were included accidentally.
- [ ] Search release artifacts for credentials, tokens, SSH/private keys,
      private host paths and usernames.
- [ ] Confirm bundle configs/environment summaries are redacted.
- [ ] Confirm large primary data are pointers unless explicitly approved.
- [ ] Confirm permission for every human-data label, result, screenshot or
      manuscript asset proposed for publication; private runtime evidence is not
      implicitly approved by software-release authorization.
- [ ] Confirm selected VCF/BCF/TSV artifacts and indexes are complete.

## Reproducibility and documentation

- [ ] Final JSON/YAML and Markdown/HTML reports agree on track statuses.
- [ ] Reproducibility bundle includes effective config, software versions,
      commands, environment summary, Git SHA, sample sheet and reference metadata.
- [ ] Installation, deployment, known limitations and release docs are current.
- [ ] Benchmark and manual-review statements contain no clinical interpretation.
- [ ] Deferred human review remains explicitly unreviewed; no screenshot or
      automated pre-review is counted as a signed manual-validation result.
- [ ] Historical Sawfish metrics retain their unresolved sample-identity qualifier
      in text, tables and any derived figure.
- [ ] Distinguish the regional 17-job multi-track and 13-job read-SV acceptance
      records from matched-truth accuracy benchmarks and WGS-scale validation.

## Human-authorized release actions

- [x] Preserve the historical approval for `0.1.0rc1`.
- [x] Obtain explicit approval to prepare `0.1.0rc2` from the independently
      validated DAG/runtime and GLnexus sample-order remediation.
- [x] Obtain explicit approval to prepare `0.1.0rc3` from the independently
      validated Jasmine SUPP_VEC evidence-source remediation.
- [x] Prepare `0.1.0rc4` only from the independently Linux-validated Jasmine
      header/HTSlib and PAV 2.4.6 SV-only/six-source harmonization deltas.
- [ ] Commit the reviewed release candidate.
- [ ] Create the approved annotated tag.
- [ ] Push only after explicit approval.
- [ ] Publish GitHub/PyPI/container artifacts only after explicit approval and
      post-upload checksum/install verification.
- [ ] Clone the public repository into a clean Linux path and repeat documented
      source, wheel, Conda and available container installation commands.
