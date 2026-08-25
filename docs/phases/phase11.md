# Phase 11 — Annotation / Functional Prioritization

**Status: COMPLETE**

Phase 12 is **NOT STARTED**.

## Scope and invariant

Phase 11 annotates independent small-variant, structural-variant, and tandem-
repeat sources without rewriting or combining their raw caller VCFs.

```text
variant call confidence != biological/functional impact
```

Annotation evidence and functional predictions do not establish truth,
pathogenicity, clinical significance, or improved genotype/call confidence.

## Annotation contract

`AnnotationInput` preserves sample, variant category, source tool/VCF, explicit
source variant IDs, `ReferenceGenome`, reference build, and source immutability.
`AnnotationDatabase`, `AnnotationArtifact`, and `AnnotationResult` record the
annotation source, tool/database versions, local database path, command,
runtime, output, and scientific policy. Phase 11 reports have atomic JSON/YAML
writers and retain small/SV/TR sources as separate entries.

## ANNOVAR

`AnnovarWrapper` implements the official
[table_annovar.pl workflow](https://annovar.openbioinformatics.org/en/latest/user-guide/startup/):

- executable availability through `CommandRunner`;
- explicit reference build, protocols, and paired `g`/`gx`/`r`/`f` operations;
- externally managed database root and database version;
- deterministic VCF-input command and native multianno TSV/VCF validation;
- dry-run, logs, external failure propagation, and overwrite protection.

The documented ANNOVAR interface does not provide a stable machine-readable
version command. HiFiVar therefore requires an explicit deployment release ID
and records it as a compatibility assertion. Database download is deliberately
outside the wrapper and workflow.

## Ensembl VEP

`VepWrapper` implements the documented
[offline/cache mode](https://www.ensembl.org/info/docs/tools/vep/script/vep_options.html)
with an explicit local cache directory/version, species, assembly, reference
FASTA, VCF input, tabular output, and threads. It parses the version block from
`vep --help`, executes through `CommandRunner`, validates a non-empty output,
and provides dry-run, logs, overwrite protection, and clear failures. ANNOVAR
and VEP share only the tool-neutral artifact contract; neither depends on the
other.

## SV/TR region annotation

`RegionDatabase` and `annotate_region_overlaps()` support separate versioned BED
sources for gene, exon, regulatory region, repeat, and segmental duplication
overlap. The VCF is streamed to materialize only explicitly listed source
variant IDs. The overlap table retains the source ID/tool/VCF and original
contig/start/end, and states `breakpoint_modified = false`. No breakpoint,
SVTYPE, allele, or source VCF is normalized or rewritten.

This is a bounded overlap abstraction, not AnnotSV and not a replacement for a
fully indexed large-scale interval engine. It also supplies the TR annotation
boundary without rerunning or changing TRGT.

## AlphaGenome boundary

The `FunctionalBackend` protocol and `FunctionalPrioritizationRequest` require a
non-empty explicit candidate selection, model name/version, requested modalities,
source annotation, and source variant identity. Results preserve modality-level
scores rather than reducing them to a pathogenic/benign label.

No production AlphaGenome credential or cloud backend is bundled. This matches
the official [AlphaGenome API guidance](https://github.com/google-deepmind/alphagenome),
which describes API credentials and positions the service for limited selected
variant/region prediction rather than unrestricted high-volume execution.
Windows tests inject a mock backend. Linux/cloud verification remains `NOT_RUN`.

## Configuration and workflow

`annotation.enabled` defaults to `false` and requires an explicit
`annotation.input_manifest`. ANNOVAR, VEP, and region-overlap switches are
independent. Database/cache paths and versions are configuration, never embedded
in Python or the Snakefile.

`workflow/rules/annotation.smk` adds optional deterministic directories:

```text
results/annotation/<sample>/<small|sv|tr>/
```

Each directory owns source-specific outputs plus `phase11.provenance.json` and
YAML. The rule delegates to `workflow/scripts/run_annotation.py`. Upstream VCFs
are inputs only; annotation failure cannot remove or modify them. The
AlphaGenome switch intentionally stops with a clear boundary message until a
credentialed production backend has been independently verified.

## Verification coverage

- ANNOVAR and VEP deterministic commands and provenance;
- missing database/cache, overwrite, dry-run, Unicode paths, output validation,
  version parsing, and external failure;
- small/SV/TR independent integration;
- gene/exon/regulatory/repeat/segdup overlap and immutable breakpoints;
- explicit AlphaGenome selection and exact result traceability;
- functional impact distinct from confidence/truth;
- optional Snakemake dry-run and dedicated wrapper bridge;
- full Phase 0–10 regression, compile, and package build.

## Definition of Done

- Annotation contract: complete.
- ANNOVAR integration: complete with external database provisioning boundary.
- VEP integration: complete with offline/cache boundary.
- SV/TR annotation boundary: complete.
- AlphaGenome explicit-selection interface: complete; real cloud backend not
  bundled by design.
- Provenance and raw-caller immutability: complete.
- Optional Snakemake annotation branch: complete.
- Phase 12: not started.
