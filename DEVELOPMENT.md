# HiFiVar development verification

HiFiVar is developed cross-platform, but external bioinformatics tools and the
production Snakemake workflow target Linux/HPC. Windows verification uses unit,
fake-executable, CLI, and Snakemake smoke tests; it does not claim real-tool
execution.

## Local regression

```bash
python -m pytest
python -m compileall src
python -m build
```

The test suite contains the Phase 0 infrastructure smoke regression and modular
Snakemake regressions through Phase 13.

## Phase 6 real-tool verification

Prerequisites are indexed aligned BAM, indexed `sample.small.vcf.gz`, the same
reference FASTA/FAI, HiPhase, and tabix. Confirm the official CLI on the installed
build before execution:

```bash
hiphase --version
hiphase --help
tabix --version
python -m pytest -p no:cacheprovider tests/unit/test_hiphase.py \
  tests/integration/test_phase6_complete.py \
  tests/integration/test_snakemake_phase6.py
```

`LINUX_REAL_VERIFICATION: NOT_RUN`

## Phase 7 real-tool verification

Prerequisites are a tiny PacBio HiFi FASTQ and hifiasm. Phase 7 does not require
a reference for assembly execution.

```bash
hifiasm --version
hifiasm --help
hifiasm -o work/assembly/TINY/TINY.asm -t 4 tiny.hifi.fastq.gz
test -s work/assembly/TINY/TINY.asm.bp.p_ctg.gfa
test -s work/assembly/TINY/TINY.asm.bp.hap1.p_ctg.gfa
test -s work/assembly/TINY/TINY.asm.bp.hap2.p_ctg.gfa
python -m pytest -p no:cacheprovider tests/unit/test_hifiasm.py \
  tests/integration/test_phase7_complete.py \
  tests/integration/test_snakemake_phase7.py
```

`LINUX_REAL_VERIFICATION: NOT_RUN`

Before a production run, record the executable path, exact version, environment
or container identity, scheduler resources, effective config, and Git commit.
## Phase 8-9 real-tool verification

Use only tiny assemblies and caller VCFs for the first Linux/HPC smoke. Confirm
the installed command contract before execution:

    svim-asm --version
    svim-asm haploid --help
    jasmine --version
    jasmine file_list=inputs.txt out_file=merged.vcf genome_file=reference.fa max_dist=1000
    truvari --version
    truvari bench -b merged.vcf.gz -c source.vcf.gz -f reference.fa -o truvari-out
    python -m pytest -p no:cacheprovider tests/unit/test_assembly_sv.py \
      tests/unit/test_harmonization.py tests/integration/test_phase8_complete.py \
      tests/integration/test_phase9_complete.py tests/integration/test_snakemake_phase8.py \
      tests/integration/test_snakemake_phase9.py

PAV is invoked through its Snakemake analysis workflow adapter; confirm the
site Snakefile and target against the installed PAV release. Validated runtime
contracts and deployment requirements are documented in `docs/deployment.md`.

### Validated PAV 2.4.6 deployment contract

The validated production backend for PAV is **Apptainer**. The Linux/HPC delta
revalidation completed through `PavWrapper` with PAV 2.4.6 and Snakemake 8.24.1
inside the image. Native PAV and other container frameworks are not part of
this validated contract.

The deployment uses a small executable launcher with this argument-preserving
pattern:

```bash
#!/usr/bin/env bash
set -euo pipefail
exec apptainer exec --bind "<root>:<root>" "<image>" snakemake "$@"
```

Configure `assembly_sv.pav.executable` as the host launcher and
`assembly_sv.pav.version` as `"2.4.6"`. The launcher must be on an executable
host filesystem; on the validated HPC the repository software tree was mounted
`noexec`, so the launcher lived under an exec-mounted deployment root. `<root>`
must contain or expose the analysis work directory, reference, haplotype
assemblies, output directories, launcher, and host PAV workflow copy at the
same absolute paths visible inside the container.

PAV's authoritative workflow entry point in the validated image is
`/opt/pav/Snakefile`. HiFiVar validates `assembly_sv.pav.snakefile` on the host
before launching, so the container's `/opt/pav` workflow tree must be copied to
a host path under `<root>` and that host `Snakefile` path must be configured.
Copy the complete workflow tree, not only a detached Snakefile, so its relative
rules and support files remain together. For example:

```bash
PAV_ROOT=/path/on/exec-mounted-storage/pav-2.4.6
PAV_IMAGE=/path/to/containers/pav_2.4.6.sif
mkdir -p "$PAV_ROOT/workflow"
apptainer exec --bind "$PAV_ROOT:$PAV_ROOT" "$PAV_IMAGE" \
  cp -a /opt/pav/. "$PAV_ROOT/workflow/"
test -f "$PAV_ROOT/workflow/Snakefile"
```

The real validation reused the local SIF artifact `pav_latest.sif`, originally
staged from the PAV Sylabs Library namespace
`library://becklab/pav/pav:latest`; its contained PAV release was verified as
2.4.6. A production deployment must not treat that filename or floating source
tag as identity. Store the source URI in the deployment manifest, give the
validated local artifact a versioned name such as `pav_2.4.6.sif`, record its
SHA-256, and verify the contained release. Pin HiFiVar's configured PAV version
to `2.4.6`; do not re-pull `latest` and assume it is the validated artifact.

The local SIF, executable launcher, and host copy of `/opt/pav` are deployment
artifacts. They are not HiFiVar production source, must not be placed in
`src/`, and are not shipped in the Python wheel. Site-specific absolute paths
belong in deployment/effective configuration. Whenever the container artifact
is upgraded or replaced, refresh the host workflow copy from that same image,
recalculate the image checksum, update the declared PAV version, and repeat the
Linux real-tool validation before production use. Never combine a Snakefile
copied from one image with a different container version.

PAV 2.4.6 APPTAINER REAL-TOOL VALIDATION: PASS

## Phase 10 IGV/manual-review verification

Phase 10 is optional and downstream of calling. Supply an explicit UTF-8 TSV in
`review.selection_file`; HiFiVar does not select variants using a hard-coded
caller-support threshold. Required columns are:

```text
review_id sample variant_id variant_type contig start end source_vcf source_caller evidence_class
```

Optional columns are `mate_contig`, `mate_position`, `flank_bp`, and
`trgt_visualization`. The sample must resolve to an indexed BAM/CRAM, and its
source VCF and reference FASTA/FAI must exist. Enable the branch with:

```yaml
review:
  enabled: true
  selection_file: configs/review_targets.tsv
  igv_executable: igv.sh
  flank_bp: 500
  overwrite: false
```

Windows verification uses a fake IGV executable and Snakemake dry-run. On a
Linux workstation or HPC visualization node, inspect the installed launcher and
then run the generated batch in an environment with the required Java runtime
and graphical display:

```bash
igv.sh --version
igv.sh --help
python -m pytest -p no:cacheprovider tests/unit/test_review.py \
  tests/integration/test_phase10_complete.py \
  tests/integration/test_snakemake_phase10.py
igv.sh --batch results/review/review.igv.batch
test -s results/review/review_manifest.tsv
find results/review/screenshots -type f -name '*.png' -size +0c
```

The batch contract follows IGV's official `genome`, `load`, `goto`,
`snapshotDirectory`, `snapshot`, and `exit` commands. Real execution is
site-dependent and is not claimed by the Windows tests.

`LINUX_REAL_IGV_VERIFICATION: NOT_RUN`

## Phase 11 annotation verification

Phase 11 consumes an explicit TSV with exactly these columns:

```text
sample variant_category source_vcf source_tool source_variant_ids
```

`variant_category` is `small`, `sv`, or `tr`. Semicolon-separated
`source_variant_ids` are optional for ANNOVAR/VEP and required for bounded
SV/TR region-overlap annotation. Source paths are resolved beside the manifest.

ANNOVAR uses the official `table_annovar.pl <input> <humandb> -buildver ...
-out ... -protocol ... -operation ... -vcfinput -polish` contract. ANNOVAR has
no stable machine-readable version flag in its documented quick-start contract,
so `annotation.annovar_version` is a mandatory deployment release identifier.
The database root and `annotation.annovar_database_version` are also mandatory;
HiFiVar never downloads databases.

VEP uses local `--offline --cache --dir_cache` execution with explicit cache
version, species, assembly, reference FASTA, and tabular output. The wrapper
reads its release from the documented version block emitted by `vep --help`.

Linux/HPC smoke verification:

```bash
table_annovar.pl
vep --help

python -m pytest -p no:cacheprovider \
  tests/unit/test_annotation.py \
  tests/integration/test_phase11_complete.py \
  tests/integration/test_snakemake_phase11.py

snakemake --snakefile workflow/Snakefile \
  --configfile effective.phase11.yaml --cores 1 --dry-run
snakemake --snakefile workflow/Snakefile \
  --configfile effective.phase11.yaml --cores 4 results/annotation/S1/small
```

Record the ANNOVAR release, each protocol/database version, VEP release/cache,
reference build/checksum, effective config, commands, Git commit, and logs.

AlphaGenome remains a clean injected-backend boundary. The official API requires
external credentials and is intended for selected variants rather than
unbounded whole-genome requests. Linux/cloud real verification must provide a
credential out of band, explicitly set the model/API version and requested
modalities, and never store the credential in config or provenance.

`ANNOVAR_LINUX_REAL_VERIFICATION: NOT_RUN`

`VEP_LINUX_REAL_VERIFICATION: NOT_RUN`

`ALPHAGENOME_CLOUD_VERIFICATION: NOT_RUN`

## Phase 12 cohort verification

Phase 12 consumes a long-form TSV with columns:

```text
sample track state source_path index_path source_tool source_version reference_build catalog_id
```

Each enabled track must contain exactly one row for every cohort sample in
sample-sheet order. Small variants require indexed DeepVariant gVCFs; SV
prefers Phase 9 harmonized VCFs; TR requires one explicit catalog identity.
Missing records remain `NOT_OBSERVED`, never `0/0`.

GLnexus BCF stdout is persisted directly by `CommandRunner`, followed by
separate bcftools commands for BCF/VCF indexing and compressed VCF conversion;
no shell pipeline is used.

```bash
conda env create -f workflow/envs/glnexus.yaml
conda activate hifivar-glnexus-1.4.1
glnexus_cli --help
bcftools --version
python -m pytest -p no:cacheprovider tests/unit/test_cohort.py \
  tests/unit/test_glnexus.py tests/unit/test_cohort_tracks.py \
  tests/unit/test_phase12.py tests/integration/test_phase12_complete.py \
  tests/integration/test_snakemake_phase12.py
snakemake --snakefile workflow/Snakefile --configfile effective.phase12.yaml \
  --cores 8 --use-conda --keep-going --dry-run
```

`GLNEXUS_LINUX_REAL_VERIFICATION: NOT_RUN`

## Phase 13 benchmark verification

Provision tiny, reference-compatible query/truth VCFs and confident regions;
HiFiVar does not download benchmark data. Record exact resource versions and
checksums. Confirm the installed CLIs before real execution:

```bash
hap.py --version
hap.py --help
truvari version
truvari bench --help
python -m pytest -p no:cacheprovider tests/unit/test_benchmark.py \
  tests/unit/test_happy.py tests/unit/test_truvari_benchmark.py \
  tests/unit/test_benchmark_config.py tests/integration/test_snakemake_phase13.py
snakemake --snakefile workflow/Snakefile --configfile effective.phase13.yaml \
  --cores 8 --dry-run
snakemake --snakefile workflow/Snakefile --configfile effective.phase13.yaml \
  --cores 8 results/benchmark/TINY/benchmark_manifest.json
```

PASS requires non-empty validated hap.py/Truvari outputs, parsed metrics with
the declared truth/region versions, unchanged query checksums, and manifest
provenance. Missing truth resources must remain `NOT_RUN`, not zero metrics.

`HAPPY_LINUX_REAL_VERIFICATION: NOT_RUN`

`TRUVARI_PHASE13_LINUX_REAL_VERIFICATION: NOT_RUN`

## Phase 14 release-hardening verification

Run the release candidate from a clean source state and inspect all artifacts:

```bash
python -m pytest -p no:cacheprovider
python -m compileall src
python -m build
python -m zipfile -l dist/hifivar-0.0.1.dev0-py3-none-any.whl
```

Then follow `docs/installation.md` for an out-of-tree clean virtual environment
and verify packaged YAML/workflow access. Review `docs/deployment.md` rather
than assuming every wrapper has a real-tool PASS. Complete
`RELEASE_CHECKLIST.md` before requesting any version bump, tag, push or publish.

CI intentionally excludes heavy tools, databases, containers, WGS data and
scientific benchmark runs. Those are explicit Linux/HPC evidence gates.
