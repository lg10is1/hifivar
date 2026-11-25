# Plot one existing TRGT locus

[简体中文](zh_CN/trgt_plot.md)

This is an **unreleased development-source consumer**, not a feature in
RC1–RC4 or an earlier package merely reporting the same RC5 version. Verify the
exact source/package manifest and its acceptance scope. It does not run TR
genotyping, change the calling DAG, execute IGV or perform a human review.

## Inputs and supported scope

Use the VCF, its TBI, spanning-read BAM and its BAI from the **same original
TRGT analysis**, together with that analysis's reference FASTA/FAI and repeat
catalog. The command supports native **TRGT 5.1.0**, with an exact version report
of `5.1.0` or the observed deployment build `5.1.0-ec66463`, one exact `TRID:<catalog ID>`,
`allele`/`motifs` plots and SVG only. Methylation, waterfall, deepdive and arbitrary
native arguments are not exposed. There is no threads parameter for this route.

The plan must supply SHA256 for exactly seven roles: `fasta`, `fai`, `catalog`,
`vcf`, `tbi`, `spanning_bam`, `bai`. These hashes bind file bytes; they do **not**
authenticate a producer, establish that the BAM came from that VCF, or prove
sample ancestry. The initial source mode is deliberately `EXTERNAL_FILES` only.
Keep the original producer records separately; do not replace them with a
user-entered attempt identifier.

The consumer checks one declared VCF sample and one unique TRID, matched catalog
coordinates, motifs, structure, complete GT and explicit sequence alleles. It
reads the selected anchor-plus-repeat REF from the actual FASTA and checks it
against the VCF. FAI must be the conventional `FASTA.fai`; VCF/TBI and BAM/BAI
also use their explicit conventional names. Reference and catalog builds must
match. No silent contig rename, allele rewrite or genotype adjustment is made.
These are finite input/native-parser checks, not a full BAM lineage audit or a
matched-truth accuracy benchmark.

## Prepare a strict plan

Copy [the generic example](../examples/review/trgt_plot.yaml) into your project,
then replace **all** paths, the sample/locus and seven checksum placeholders.
Do not run the unedited template. Obtain hashes from your actual files, for
example:

```bash
sha256sum reference.fa reference.fa.fai repeats.bed sample.tr.vcf.gz sample.tr.vcf.gz.tbi sample.tr.spanning.bam sample.tr.spanning.bam.bai
```

The schema has these required sections:

- `schema_version: 1`, `sample_id`, `selector: TRID:<raw catalog ID>`;
- `reference: {fasta, fai, build}` and `catalog: {path, build}`;
- `source: {mode, vcf, tbi, spanning_bam, bai}`;
- `expected_sha256`: exactly the seven lowercase 64-character hashes;
- `runtime: {executable, supported_version}` with version string `"5.1.0"`;
- `plot: {type, show, format, max_allele_reads}` set to allele/motifs/svg and a
  positive read-display limit;
- `execution: {timeout_seconds}`: finite, positive, at most 1200 seconds;
- `output: {directory}`: one new leaf directory whose parent already exists.

Relative paths resolve from the plan's directory, not the terminal's current
directory. Unknown, missing, duplicate, merged, aliased or non-string YAML keys
are rejected. Native flags cannot be embedded in `runtime.executable`. The
loader binds the original plan's SHA separately from the normalized request.
The timeout is per version/plot command, not a total consumer wall-time limit;
input validation and repeated hashing also take time.

## Preview, then render

```bash
hifivar --dry-run review tr-plot --plan trgt_plot.yaml
hifivar review tr-plot --plan trgt_plot.yaml
```

Dry-run fully validates the selected files and hashes but does not probe or run
TRGT and does not create the output directory. Hashing a whole reference can
take time even for a one-locus plot. To actually render, install TRGT 5.1.0
separately and make the plan's executable available.

The actual command uses `CommandRunner` and an exclusively claimed fresh output
directory. It records separate version and plot stdout/stderr, the effective
inputs, native argv, tool identities and file hashes. The version logs retain
their raw bytes; the receipt records normalized `tool_version: 5.1.0` and the
complete observed token in `tool_version_reported`. Unknown suffixes, additional
version-line annotations and conflicting reports across stdout/stderr are rejected, rather than stripped or
silently accepted. The plan's `runtime.supported_version` remains `"5.1.0"`.
Support for this exact deployment token is not native plot acceptance of a
new package; that package still requires its own execution evidence. The plot is `plot.svg`.
The normal successful CLI status is `COMPLETE_TECHNICAL_SVG`; `PLANNED` is only a
preview. Neither means scientific accuracy or manual review is complete.

## Completion and failure

The fixed evidence set includes `effective_input.json`, `plot.svg`, separate
`version.*.log` and `command.*.log`, the retained `.completion.pending.json`,
and the last-published `completion.json`. Source/tool/plan and output identities
are checked again. An error may additionally leave `failure.json`.

Accept a run only when the process returned **exit 0**, its completion receipt
matches the actual files, and **no `failure.json` is present**.
File existence alone is insufficient. An interruption immediately after the
receipt link can leave a receipt even though the command did not return
normally. Do not silently treat that directory as accepted or resume into it.

There is no overwrite/resume switch, automatic retry, deletion or recursive
cleanup. Existing and late-arriving destinations fail closed. Keep a failed
directory for diagnosis and use a new, deliberately selected directory only
after resolving the cause. The directory is not an OS sandbox against hostile
concurrent mutation, and receipt publication is not a multi-file physically
atomic transaction or a power-loss recovery guarantee.

SVG output is structurally checked against a conservative static whitelist;
scripts, dynamic content and external resources are rejected. Structural
validation is not pixel QA: open the SVG in a trusted viewer and inspect text,
read tracks and layout before using it. Keep the native SVG unchanged when
making manuscript presentation derivatives.

## Review remains a separate evidence layer

The receipt always preserves `manual_review: NOT_REVIEWED` and
`accuracy: NOT_ESTABLISHED`. Native allele-panel order is **not** paternal versus
maternal order and does not establish phasing accuracy. Repeat lengths and
motif patterns in the image are not clinical interpretations.

To link the generated SVG into a later offline review bundle, explicitly place
its path in the existing `trgt_visualization` selection column described in the
[review guide](review_cli.md). This does not execute IGV, edit old bundles or
complete the user's manual decisions. The plot/effective plan may contain
sample names and absolute paths; keep them private unless separately cleared
for sharing.
