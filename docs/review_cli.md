# Offline review preparation and decision import

[简体中文](zh_CN/review_cli.md)

The newer unreleased source also offers a separate
[`review tr-plot --plan` consumer](trgt_plot.md) for one existing TRGT 5.1.0
locus. It does not modify an IGV bundle or perform human review. Its exact
candidate/native acceptance is separate from the older offline CLI below.

This RC5 candidate extension provides an offline entry to the existing Phase 10
review module. It is included in the
[identified 15 October RC5 payload](rc5_candidate.md#validated-development-payload--已验收开发包),
not RC1–RC4 or the early security-only RC5 build. Select a candidate by its
source manifest and artifact SHA256, not its version alone.
[中文指南](zh_CN/review_cli.md).

## 1. Prepare a review bundle

Use an existing reference, sample TSV and an explicitly selected variant TSV.
Reference build must be declared; no contig or sample renaming occurs. Existing
aligned BAM/CRAM inputs need their indexes. Raw FASTQ/uBAM is not review input.
Selection paths resolve beside the user YAML; source VCF/optional TR visualization
paths resolve beside the selection TSV; alignments resolve beside the sample TSV.

```yaml
reference:
  fasta: /your/reference/GRCh38.fa
  build: GRCh38
samples:
  sheet: /your/analysis/samples.tsv
review:
  enabled: false
  selection_file: /your/analysis/selected.tsv
  igv_executable: /your/tools/igv.sh
  flank_bp: 500
  overwrite: false
```

The explicit `review prepare` command does not require `review.enabled=true`:
that switch controls the separate Snakemake review rule. No calling rule or
IGV process is invoked by preparation. Use a new directory in your project:

```bash
hifivar --config review.yaml --dry-run review prepare --output-dir /your/project/review/bundle01
hifivar --config review.yaml review prepare --output-dir /your/project/review/bundle01
```

Dry-run reads/checks inputs and prints a plan but writes no bundle or log file.
Preparation creates the batch, pending JSON/YAML/TSV manifest, effective config,
an editable `manual_decisions.tsv`, an empty screenshots directory and
`review_bundle.json`. All targets remain **NOT_REVIEWED**. No screenshot is
claimed generated. The receipt is written last: an absent receipt means an
incomplete bundle. Existing output directories are refused, with no overwrite,
automatic cleanup or retry. Keep failed preparation evidence and choose a new,
deliberately named run directory after diagnosis.

The required selection columns are:

```text
review_id sample variant_id variant_type contig start end source_vcf source_caller evidence_class
```

Use actual tabs, not the spaces above. Optional columns are `mate_contig`,
`mate_position`, `flank_bp`, `trgt_visualization`. IDs must be unique. Coordinates
are 1-based closed; BND needs mate coordinates and creates two planned screenshots.
SNV/INS windows use the anchor. Existing Phase 10 window semantics are retained;
selected coordinates are checked against FAI lengths, but a right-hand flank can
extend beyond the contig end and must be inspected in the viewer.

## 2. Capture and inspect evidence explicitly

Start IGV in a separately provisioned, suitable display environment and use
**Tools > Run Batch Script** to select `review.igv.batch`, or use the generated
shell-free command preview. IGV documents the `--batch` option and its batch
commands in the [desktop batch guide](https://igv.org/doc/desktop/UserGuide/tools/batch/)
and [command-line guide](https://igv.org/doc/desktop/UserGuide/advanced/command_line/).
Preparation itself does not start IGV, set up displays or verify generated images.

Choose a real reviewer. Do not convert screenshot availability or caller support
into a human decision. Keep reviewer comments and screenshots as private evidence
unless their sharing is authorized.

## 3. Import a reviewer-completed table

Copy the template to a decision TSV and fill the exact columns:

```text
review_id sample variant_id status reviewer reviewed_at comment
```

Allowed statuses are `NOT_REVIEWED`, `SUPPORT`, `NOT_SUPPORT`, `UNCERTAIN`.
These retain existing Phase 10 observation semantics: they are not truth,
clinical classification, ACMG labels or PCR confirmation. `PASS`/`REJECT` are
not silently converted. A decided row requires a nonblank reviewer and comment,
plus a timezone-aware ISO timestamp such as `2026-10-14T10:00:00+08:00`.
This is a format example, not a completed review. NOT_REVIEWED fields stay blank.

```bash
hifivar --dry-run review import-decisions --bundle-dir /your/project/review/bundle01 --decisions /your/project/review/reviewer.tsv --output-dir /your/project/review/adjudicated01
hifivar review import-decisions --bundle-dir /your/project/review/bundle01 --decisions /your/project/review/reviewer.tsv --output-dir /your/project/review/adjudicated01
```

The importer matches `review_id`, sample and variant ID exactly, rejects duplicate
or unknown IDs and invalid status/date/comment fields, and verifies the sealed
bundle metadata. Missing rows remain NOT_REVIEWED; an empty table completes no
reviews. New adjudicated JSON/YAML/TSV and an import receipt are written outside
the original bundle, retaining declared order and source paths. The original
prepared manifest and raw variant/alignment/reference files are not rewritten.

## Validation and privacy boundary

- Preparation reuses lightweight Phase 10 validators. It does not parse full BAM
  records, prove BAM header/sample identity, locate the selected allele/ID in a
  full VCF, or establish call accuracy. Complete those checks before human review.
- Raw inputs are recorded with size/mtime metadata, **not full-file SHA256**.
  Selection and generated metadata are hashed. Metadata checks do not prove
  byte-for-byte raw-data immutability; checksums are not authentication/signatures.
- Decision import does not need the large raw files and does not revalidate them
  or screenshots. It records a supplied reviewer identity without authenticating
  that person or verifying that a review occurred. Never label import success
  as independent human-review verification.
- Effective config redacts semantic secret keys, but paths, comments and sample
  names can still be sensitive. Do not embed credentials in paths/free text or
  upload raw bundles by default. TSV is plain data, not a spreadsheet executable.
- The offline CLI has dedicated bounded Linux, regression and package checks
  for the identified development payload. These do not execute IGV or complete
  human adjudication; earlier screenshot capture is a separate scope. No
  scientific performance result is added.

## Junction-aware BND source validation

Development BND targets now bind their exact selected source record and ALT endpoints. Single/virtual endpoints remain explicit, with no invented screenshot. The optional bnd_bundle_manifest column is required for an explicitly derived BND source; standalone derived VCFs fail closed. See [junction consumers](bnd_junction_consumers.md). Preparation is not manual review.
