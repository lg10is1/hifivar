# Phase 8  Assembly-based structural-variant calling

**Phase 8 status: COMPLETE**

## 8.18.3 Contract and PAV boundary

`AssemblySvRequest`, `AssemblySvArtifact`, and `AssemblySvCollection` keep
assembly evidence separate from Phase 4 read evidence. Every artifact records
sample, caller, reference build, explicit haplotype FASTAs, raw VCF, finalized
BGZF/TBI, intermediates, commands, backend, and caller version.

PAV is modeled as a workflow adapter, not a fictitious single caller binary.
HiFiVar atomically writes PAV's official `config.json` and `assemblies.tsv`
inputs and invokes a configured PAV Snakefile through `CommandRunner`.
`HAP_h1` and `HAP_h2` remain independent inputs. Raw PAV output is retained.

The adapter runs PAV's default DAG instead of requesting HiFiVar's final VCF
as an upstream Snakemake target.  It validates PAV's native
`{sample}.vcf.gz` and TBI in the analysis directory, retains both files, then
atomically copies them to `{sample}.pav.assembly.sv.vcf.gz` and its index.
This contract addresses audit finding P8-PAV-001 without changing PAV itself.
The command deliberately has no explicit output target: PAV 2.4.6's default
`rule all` owns its internal target names. Real execution also requires an
explicit numeric `assembly_sv.pav.version`; the unresolved placeholder remains
valid for planning only and cannot be recorded as a completed caller version.
When `overwrite: true`, HiFiVar atomically replaces only its owned
`config.json`, `assemblies.tsv`, and captured outputs. It does not clear PAV's
work directory.

- PAV implementation: PASS
- PAV mock verification: PASS
- PAV upstream 2.4.6 tiny Linux/HPC validation: PASS (554/554 jobs)
- PAV adapter delta verification after P8-PAV-001: PASS (Linux/HPC, Apptainer)

### Production deployment contract

The real adapter validation fixes the supported production deployment at PAV
**2.4.6** through **Apptainer**. The deployment launcher preserves the
wrapper's Snakemake arguments:

```text
apptainer exec --bind <root>:<root> <image> snakemake "$@"
```

The validated image contains the upstream workflow at
`/opt/pav/Snakefile`. Because `PavWrapper` checks its configured Snakefile as a
host file before execution, deployment must copy the complete container
`/opt/pav` workflow tree to a host path below the bind root and configure that
host-side `Snakefile`. This keeps the Snakefile and its relative rule/support
files from the same PAV release together. The copy must be refreshed from the
container whenever the image is upgraded; mixing a stale host workflow copy
with a newer image is unsupported.

The validated local artifact was named `pav_latest.sif` and originated from
`library://becklab/pav/pav:latest`, but production identity is not derived from
that filename or floating tag. Deployment records must retain the source URI,
validated contained version 2.4.6, local versioned SIF path, and SIF SHA-256;
`latest` must not be re-pulled and assumed equivalent. The SIF, host workflow
copy, and executable launcher are site deployment artifacts—not files under
`src/` and not Python package contents. Site paths belong in effective
configuration. Native PAV and other container frameworks have not been
validated by this contract.

## 8.48.5 SVIM-asm and haplotypes

SVIM-asm uses the official alignment contract explicitly: each haplotype FASTA
is aligned separately with minimap2, sorted/indexed with samtools, and supplied
to `svim-asm haploid` or `svim-asm diploid`. FASTAs are never concatenated.
The native `variants.vcf` is retained before explicit bgzip/tabix finalization.

- SVIM-asm implementation: PASS
- SVIM-asm mock verification: PASS
- SVIM-asm official main CLI verification: PASS
- SVIM-asm Linux real verification: NOT_RUN

## 8.68.7 Validation and provenance

Validation checks the expected VCF/TBI, BGZF and tabix signatures, one exact
sample, INFO/SVTYPE, and exact reference-contig subset compatibility. It does
not rewrite records, contigs, breakpoints, SVTYPE, BND, or INS alleles. Run
reports are atomic JSON/YAML and do not hash large assemblies implicitly.

## 8.88.9 Workflow and integration

The config-disabled `assembly_sv.smk` module branches from Phase 7 hifiasm
outputs into independent `pav_assembly_sv` and `svim_asm_call` rules.
Unit, tiny fake integration, Unicode-path, haploid/diploid, failure, missing
input/output, overwrite, validation, provenance, and Snakemake dry-run tests
cover this boundary.

## Output policy

- `sample.pav.assembly.sv.vcf.gz`
- `sample.svim_asm.assembly.sv.vcf.gz`

These are independent evidence streams. Phase 8 creates no merged or final SV
callset and does not execute Jasmine or Truvari.

## External verification

PAV's official documentation confirms the analysis-directory
`config.json`/`assemblies.tsv` contract and independent haplotype columns. The
PAV 2.4.6 Apptainer adapter path, default DAG, native output discovery,
overwrite behavior, and version provenance have now passed Linux real-tool
delta validation.
SVIM-asm's official documentation confirms separate assembly-to-reference BAMs
and haploid/diploid modes. Installed SVIM-asm `--sample` placement remains a
separate tool-contract concern.

## Intentionally out of scope

Jasmine, Truvari harmonization, SURVIVOR, truth/confidence labeling, benchmark,
annotation, IGV/manual review, and Phase 10 are absent.
