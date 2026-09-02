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
as an upstream Snakemake target. It validates and retains PAV's native mixed
`{sample}.vcf.gz` and TBI in the analysis directory. It then streams a separate
SV-only plain VCF and finalizes it with configured bgzip/tabix executables as
`{sample}.pav.assembly.sv.vcf.gz` plus TBI.
This contract addresses audit finding P8-PAV-001 without changing PAV itself.
The command deliberately has no explicit output target: PAV 2.4.6's default
`rule all` owns its internal target names. Real execution also requires an
explicit numeric `assembly_sv.pav.version`; the unresolved placeholder remains
valid for planning only and cannot be recorded as a completed caller version.
When `overwrite: true`, HiFiVar replaces only its owned `config.json`,
`assemblies.tsv`, SV-only intermediate, and finalized output paths. It does not
clear PAV's work directory or rewrite PAV's native mixed VCF and index.

### Post-RC3 real-data handoff remediation

Whole-genome validation showed that PAV's native root `{sample}.vcf.gz` is a
mixed callset containing SNVs, short indels, and structural variants. PAV 2.4.6
has no official SV-only VCF target. Its own VCF writer nevertheless defines the
authoritative VARTYPE boundary: `sv_inv` records are structural, while
`svindel_ins` and `svindel_del` records are structural when `SVLEN >= 50`
before PAV negates deletion lengths in VCF output.

HiFiVar now mirrors that exact, version-locked rule by retaining INV records
and INS/DEL records with absolute `SVLEN >= 50`. Selection is streaming and
keeps PASS and non-PASS records, original IDs, alleles, FILTER, genotype,
haplotype, coordinates, and INFO fields. The mixed root VCF and its index remain
immutable provenance. Unknown SVTYPE, missing/scalar-invalid SVLEN, missing PAV
source version, or a non-2.4.6 PAV release fails explicitly. This is not an ad
hoc HiFiVar scientific threshold; it is the PAV 2.4.6 VARTYPE implementation.

The SV-only adapter delta passed Linux real-data validation against the
SCZ_BC2003 PAV 2.4.6 mixed root VCF: 29,962 eligible records were selected with
zero missing, extra, short-indel, SNV, or source-field differences. The derived
artifact is authorized as the sixth Phase 9 Jasmine source.

- PAV implementation: PASS
- PAV mock verification: PASS
- PAV upstream 2.4.6 tiny Linux/HPC validation: PASS (554/554 jobs)
- PAV adapter delta verification after P8-PAV-001: PASS (Linux/HPC, Apptainer)
- PAV mixed-to-SV-only adapter delta: LINUX REAL-DATA PASS

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
callset and does not execute Jasmine or Truvari. The PAV output is now a derived
SV-only artifact whose immutable raw evidence remains the native mixed root
VCF. Six-source downstream use is enabled only for this validated derived
artifact, never for the mixed PAV root VCF.

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
