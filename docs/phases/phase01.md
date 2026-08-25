# Phase 1: Core Bioinformatics Models

**Phase 1 status: COMPLETE**

## Phase 1.1: ReferenceGenome — COMPLETE

- Immutable `Contig` metadata with exact name and positive length
- Immutable `ReferenceGenome` metadata for FASTA, required FAI, explicit build,
  ordered contigs, and optional SHA256
- Lightweight FASTA validation plus FAI parsing reused from the Phase 0 API
- Exact contig compatibility with no automatic chromosome-name conversion
- JSON/YAML-friendly summary and optional detailed serialization
- Optional checksum calculation that is disabled during ordinary construction
- Nullable `reference.fasta` and `reference.build` configuration fields

Phase 1.1 does not download references, generate FAI files, manage sequence
dictionaries, parse BAM headers, or add biological workflow rules.

## Phase 1.2: Sample/InputDataset — COMPLETE

- `InputType` enum limited to FASTQ, BAM, and CRAM primary inputs
- Immutable `InputDataset` with suffix inference, exact type checking, ordered
  paths, duplicate rejection, and JSON/YAML-friendly serialization
- One or more long-read FASTQ files, or exactly one BAM/CRAM file
- Optional explicit BAM/CRAM index validation without index generation
- Immutable `Sample` with one primary dataset and strict ASCII-safe machine ID
- Unicode filesystem paths retained without allowing Unicode machine sample IDs
- Lightweight input validation reused from the Phase 0 API

Phase 1.2 does not parse BAM/CRAM headers, compute FASTQ QC or checksums, attach
derived artifacts, load sample sheets, or add workflow wildcards.

## Phase 1.3: SampleSheet/Pedigree/Cohort Metadata — COMPLETE

- Strict UTF-8/UTF-8-BOM TSV schema with required and optional canonical columns
- Ordered multi-sample loading through existing `Sample` and `InputDataset`
- Sheet-relative input paths, semicolon multi-FASTQ parsing, and input uniqueness
- Declared `Sex`, parent IDs, phenotype, and group metadata
- Parent existence, self/duplicate-parent, declared-sex, and cycle validation
- Partial and multi-generation pedigrees plus complete trio extraction
- JSON/YAML-friendly source and ordered record serialization

Phase 1.3 does not infer sex, parse alignment headers, validate references, or
perform trio/cohort calling. Config, CLI, and Snakemake are not integrated yet.

## Phase 1.4: AnalysisContext/RunManifest integration — COMPLETE

- `AnalysisContext` for single-sample and ordered SampleSheet modes
- Run-wide sample-ID and primary-input uniqueness plus obvious reference/config
  conflict detection
- Explicit FASTQ `not_applicable` and BAM/CRAM `not_checked` compatibility states
- Config-driven reference/sample-sheet construction with source-aware relative
  paths and nullable packaged defaults
- Versioned `RunManifest` with absolute non-resolved paths, UTC timestamp,
  effective config, file sizes, optional input checksums, and recursive secret
  redaction
- Atomic UTF-8 JSON/YAML output with explicit overwrite protection and readback
- Unicode-path end-to-end integration from config through both manifest formats

Phase 1.4 does not parse BAM/CRAM headers, require alignment indexes, compute
checksums by default, create run directories, launch Snakemake, add sample
wildcards, or implement QC/alignment/variant-calling modules. CLI and Snakemake
sample integration remain future work rather than Phase 1 requirements.
