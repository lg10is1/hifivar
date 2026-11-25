# Phase 13 — Benchmark / Truth-set Evaluation

[中文阶段说明](../zh_CN/phases/phase13.md)

Status: implementation complete; Linux/HPC real-tool validation complete for
the validated hap.py, Truvari, and TR comparison boundaries.

## Scope and scientific contract

Phase 13 is an optional downstream evaluation layer. It never modifies query
VCFs, generates truth, changes call confidence, or assigns pathogenicity. Every
result records the sample, reference build, immutable query, explicit truth-set
source/version, tool/version, command, outputs, metrics, and status. Truth and
region versions must be pinned; `latest` is rejected by the public model.

Statuses are `PASS`, `PARTIAL`, `NOT_RUN`, and `UNSUPPORTED`. A missing truth
resource or unavailable tool is not encoded as zero precision/recall. Metrics
are meaningful only for the declared truth set, confident regions, matching
policy, sample and reference.

## Small variants: hap.py

`HappyWrapper` uses only `CommandRunner` and the official truth/query,
`-f` confident BED, `-r` reference, `-o` prefix, threads, and engine boundary.
The parser locates SNP/INDEL and PASS (or explicitly configured filter) rows by
column names. It reads `METRIC.Recall`, `METRIC.Precision`, and
`METRIC.F1_Score`; when older summaries omit F1, it calculates F1 from recall
and precision. It never relies on row numbers. Optional stratification BEDs
carry a name, version, and region class and are written to hap.py's explicit
stratification TSV. HiFiVar does not download GIAB data or stratification BEDs.

Some packaged hap.py 0.3.15 builds return an empty `--version`. HiFiVar accepts
an explicit `benchmark.small_variants.happy_version` only after the executable
successfully runs its version command; provenance distinguishes `command` from
`config`. An absent parseable version and absent explicit version remains a
hard error. The configured value is never inferred or fabricated.

hap.py metrics discovery accepts exactly one of `<prefix>.metrics.json.gz` and
`<prefix>.metrics.json`, validates and parses either JSON form, and rejects
missing or ambiguous dual artifacts. The workflow defaults to the Linux-
observed gzip contract; `metrics_compression: plain` retains explicit support
for a legitimate uncompressed installation.

Production target: hap.py 0.3.15. Linux confirmed the empty version output and
gzip artifact behavior, and the remediation passed real delta revalidation.
`HAPPY_LINUX_REAL_VERIFICATION: PASS`.

## SV and assembly-SV: Truvari

The existing Phase 9 `TruvariWrapper` is reused and remains backward
compatible. `TruvariThresholds` optionally pins `refdist`, sequence/size/overlap
fractions, size bounds, BND distance and PASS-only policy. An explicit region
BED is passed via `--includebed`. Official `summary.json` keys are parsed as
`TP-base`, `TP-comp` (reported by HiFiVar as `tp_call`), `FP`, `FN`, precision,
recall and F1.

Additional SVTYPE and user-supplied size-bin summaries stream Truvari's assigned
tp-base/tp-comp/fp/fn VCFs. This is descriptive aggregation, not a new matching
algorithm. BND/TRA/CPX/CTX/unresolved events remain type-countable but are
explicitly excluded from unsafe length bins, yielding a partial result where
applicable. Region classes are evaluated by separate, explicitly configured
Truvari runs. Assembly-derived VCFs use the same engine but retain the distinct
`assembly_sv` variant class and provenance.

`TRUVARI_PHASE13_LINUX_REAL_VERIFICATION: PASS` for the explicitly documented
truth, region, and threshold contracts. Public performance claims remain
restricted to the validated benchmark scopes.

## Tandem repeats

The TR boundary compares single-sample TRGT-like VCFs only when truth and query
declare the same explicit catalog identity. Loci join by TRID; no locus-matching
algorithm is invented. The streaming comparison reports truth loci, compared
loci, exact allele-field agreement, genotype agreement and query no-call rate.
It does not define disease thresholds or clinical interpretations. A truth set
with no loci is `UNSUPPORTED`, not zero performance.

## Outputs and workflow

`BenchmarkManifest` writes JSON/YAML, long-form TSV metrics, and an optional
Markdown summary with atomic overwrite protection. Small/SV/assembly-SV/TR
files remain separate. `workflow/rules/benchmark.smk` supplies independent,
config-disabled tracks and a manifest rule. Benchmark targets are appended only
when explicitly enabled and are never dependencies of calling, harmonization,
review, annotation, or cohort outputs.

## Explicit limitations

- Windows tests use fake/mocked hap.py and Truvari plus Snakemake dry-run.
- Real truth resources, confident BEDs, hap.py and Truvari are externally
  provisioned on Linux/HPC; core code downloads none of them.
- The implementation does not benchmark phasing, methylation, CNV truth,
  pathogenicity, or clinical utility.
- Complex/BND length stratification is unsupported; records are not silently
  forced into a size bin.
- TR allele comparison supports explicit TRGT `AL` (with `MC` fallback) and GT;
  richer motif/methylation comparison remains future work.
- Phase 14 reporting and release hardening are complete.

## Definition of Done

- benchmark contracts and provenance: complete
- hap.py wrapper/parser/dry-run/validation: complete
- Truvari threshold, summary and stratification reuse: complete
- TR and assembly-SV boundaries: complete
- independent Snakemake tracks: complete
- unit/fake integration/full regression/build: see validation handoff
- Linux/HPC real-tool verification: complete for the documented validation
  scopes; see `docs/benchmarking.md` and `docs/reproducibility.md`
