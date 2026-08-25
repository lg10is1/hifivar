# Phase 10 — IGV / Manual Variant Review

**Status: COMPLETE**

Phase 11 annotation is **NOT STARTED**.

## Scope

Phase 10 adds an optional downstream review path for explicitly selected
artifacts from Phases 3–9. It does not choose high-confidence variants, alter a
caller result, rewrite a VCF, or make review a prerequisite of variant calling.

```text
explicit selection + BAM/CRAM + reference + source VCF
                         |
                         v
                    ReviewTarget
                         |
                         v
             deterministic IGV batch
                         |
                         v
              screenshots + evidence
                         |
                         v
             manual-review manifest
```

The implementation follows IGV's documented
[batch commands](https://igv.org/doc/desktop/UserGuide/tools/batch/) and
[command-line batch option](https://igv.org/doc/desktop/UserGuide/advanced/command_line/).

## Public API

- `hifivar.review.VariantClass`: `SNV`, `INDEL`, `DEL`, `DUP`, `INV`, `INS`,
  `BND`, and `TR`.
- `hifivar.review.ReviewStatus`: `NOT_REVIEWED`, `SUPPORT`, `NOT_SUPPORT`, and
  `UNCERTAIN`.
- `ReviewLocus`, `ReviewTarget`, and `ReviewEvidence` preserve coordinates,
  sources, and screenshot provenance.
- `ReviewResult` and `ReviewManifest` provide JSON/YAML/TSV-friendly manual
  review records and atomic, no-overwrite writers.
- `read_review_selection()` accepts an explicit TSV selection.
- `hifivar.igv.IgvWrapper` plans or executes deterministic batch mode only
  through `CommandRunner`.
- `hifivar.phase10.run_phase10()` creates the evidence bundle and manifests.

## Locus and artifact contract

Coordinates are 1-based closed intervals. SNVs use their position. Insertions
use their anchor rather than an arbitrary `END`. INDEL, DEL, DUP, INV, and TR
use their declared interval. BND produces a primary and mate locus, each with
its own screenshot. The generic `review.flank_bp` is deliberately simple and
configurable; it is not presented as a benchmark-derived optimum.

Each screenshot has the stable name
`<review_id>.locus<ordinal>.png` under `results/review/screenshots/`. Unique
review IDs and explicit overwrite protection prevent collisions. Each evidence
record links the source VCF, caller/evidence class, alignment, reference, IGV
batch, loci, screenshots, and optional pre-existing TRGT visualization.

The IGV batch contains deterministic `new`, `genome`, `snapshotDirectory`,
`load`, `goto`, `snapshot`, and `exit` operations. Python does not simulate GUI
clicks. A real run checks `igv.sh --version`, persists stdout/stderr logs,
propagates external failures, and verifies every expected non-empty screenshot.
Dry-run produces a reproducible command preview without filesystem writes or an
installed IGV.

## Selection, review, and scientific safeguards

Selection is an explicit/configurable input. There is no hard-coded
`support_count >= N` rule. Empty selections are valid and produce empty
manifests. `SUPPORT` is visual support assigned by a reviewer—not truth,
pathogenicity, clinical significance, or a replacement for benchmarking.
`UNCERTAIN` remains in the manifest and is never silently filtered. Phase 10
does not implement ACMG or any clinical-classification field.

All Phase 3–9 source VCFs, BAM/CRAM files, harmonized evidence, and TRGT outputs
are read-only inputs. Screenshot or IGV failure cannot modify those artifacts.
TRGT-specific plot/deepdive output is represented only as optional existing
evidence metadata; Phase 10 neither reruns TRGT nor implements its plotting
toolchain.

## Configuration and workflow

`review.enabled` defaults to `false`. Enabling it requires
`review.selection_file`. Other keys are `review.igv_executable`,
`review.flank_bp`, resource fields, and `review.overwrite`.

`workflow/rules/review.smk` is a modular optional branch. It resolves selected
sample alignments from `AnalysisContext`, declares selected source artifacts as
inputs, and delegates to `workflow/scripts/run_review.py`. Its deterministic
outputs are:

- `results/review/review.igv.batch`
- `results/review/screenshots/`
- `results/review/review_manifest.json`
- `results/review/review_manifest.yaml`
- `results/review/review_manifest.tsv`

Calling and harmonization targets remain independent when review is disabled or
not requested.

## Verification

Unit and fake end-to-end tests cover serialization, every supported locus type,
BND dual loci, stable screenshots, Unicode paths, overwrite refusal, empty
selection, missing inputs, dry-run, external failure, manifest retention, raw
artifact immutability, and optional Snakemake construction. Windows does not
claim a real IGV GUI/batch run. Linux verification requires an installed IGV
launcher, Java, and a suitable graphical/display environment.

## Definition of Done

- Review contract and serialization: complete.
- Variant-centred loci and BND dual-locus handling: complete.
- Deterministic IGV batch and screenshot evidence: complete.
- Manual-review manifest with non-truth semantics: complete.
- Raw Phase 3–9 artifacts remain immutable: complete.
- Optional Snakemake review branch: complete.
- Fake/mock end-to-end and regressions: complete.
- Annotation, ACMG, clinical interpretation, benchmark, and new callers: not
  implemented by design.
