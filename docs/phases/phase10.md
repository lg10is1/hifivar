# Phase 10 — IGV / Manual Variant Review

[中文阶段说明](../zh_CN/phases/phase10.md)

**Status: COMPLETE**

## Development CLI extension (14 October 2026)

The original Phase 10 history below remains unchanged. The unpublished source
now exposes `hifivar review prepare` and `hifivar review import-decisions`:
offline batch/pending-manifest preparation and separately recorded explicit
reviewer decisions. Preparation does not execute IGV or claim screenshots;
decision import does not authenticate a reviewer or establish truth. This
extension is not in historical release assets and requires separate Linux
revalidation. See [English guide](../review_cli.md) and
[中文指南](../zh_CN/review_cli.md).

Historical phase boundary: annotation was outside the original Phase 10.
Later annotation development is documented separately; the Phase 10 completion
status does not certify annotation or a manual review.

### Standalone TR plotting consumer (unreleased source)

The later explicit [`review tr-plot --plan`](../trgt_plot.md) entry consumes one
existing TRGT 5.1.0 locus and its spanning reads. It supports allele/motifs SVG
only, with exact seven-file identities and a fresh output directory. It is not
part of the IGV rule below, does not rerun TR genotyping and does not implement
deepdive, methylation or parental-phase validation. Its software/native/pixel
acceptance and package identity are separate from historical Phase 10.

## Scope

The unpublished [review source library](../review_sources.md) binds explicit
source records to existing ReviewTarget objects without commands or writes.
The separate development execution bridge now consumes this foundation through
an explicit source plan, review-only DAG, native step receipt and report track.
The library itself remains read-only; historical Phase 10 completion and old
release assets do not establish the new bridge or an actual human decision.

### Development write-path safety increment

The unpublished IGV wrapper rejects direct, resolved/symlink and hardlink
aliases between its write destinations and selected read-only inputs. The
check includes the batch staging path, all screenshots, explicit logs and the
implicit stdout log derived by CommandRunner; existing input indexes, selected
BND native/ledger bundle members, optional TR visualization and the discovered
executable are protected. Conflicting outputs/logs are rejected before mkdir,
unlink, file opening or a version command, with a second check after version
detection. Old screenshots are not removed when version detection fails.
Phase 10 applies the same guard to its JSON/YAML/TSV manifests and fixed TSV
staging path before IGV, and rechecks before writing manifests after IGV.
This finite pre-execution guard is not an OS sandbox or a guarantee against
concurrent external filesystem changes or arbitrary behavior inside IGV. It
does not change calling, review status, ALT coordinates or image rendering,
and is not retroactively present in historical release assets.

未发布源码新增IGV写路径保护：在执行/创建/删除输出之前，检查batch及临时文件、
截图、显式日志和隐式stdout日志是否与输入、索引、参考、选定BND完整来源包、
TR可视化或程序文件重合；同时检查符号链接解析和硬链接身份。版本检查后再次
核对，版本失败不先删除旧截图。这不是操作系统沙箱，不改变calling、人工判定
或断点坐标。Phase 10的JSON/YAML/TSV清单及固定TSV临时路径也在IGV执行前
检查，并在IGV完成后、写清单之前再次检查。不应将历史发布包追认为已包含该修复。

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
The original IGV workflow represents pre-existing TRGT visualization as optional
evidence metadata; it does not rerun TRGT. The separate unreleased consumer above
adds only native allele/motifs plotting. Deepdive remains unsupported, and neither
route completes human review or establishes TR accuracy.

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
