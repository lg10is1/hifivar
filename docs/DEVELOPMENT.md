# Development

HiFiVar is developed with Python 3 using a `src/` package layout. Lightweight unit, configuration, CLI, and mock tests may run on Windows; workflow integration and external bioinformatics tools target Linux/HPC.

Install the current development package and run its tests with:

```bash
python -m pip install -e ".[dev,workflow]"
python -m pytest
```

Phase 0 release-candidate checks also use:

```bash
python -m pytest --cov=hifivar --cov-report=term-missing
python -m compileall src
python -m build
```

`build` creates local distribution artifacts only; publishing is a separate,
explicit release action.

Contributions must follow `AGENTS.md` and remain within the active development phase.

## Logging

Application modules use the shared HiFiVar logging namespace:

```python
from hifivar.logging_utils import configure_logging, get_logger

configure_logging(level="INFO")
logger = get_logger(__name__)
logger.info("HiFiVar task started")
```

## Configuration

Configuration is loaded in deterministic priority order: user overrides preset,
and preset overrides defaults. The merged configuration can be written for
reproducibility without changing any source YAML:

```python
from pathlib import Path

from hifivar.config import load_config, write_effective_config

config = load_config(
    Path("src/hifivar/resources/configs/default.yaml"),
    Path("src/hifivar/resources/configs/presets/standard.yaml"),
    Path("my_config.yaml"),
)
write_effective_config(config, Path("effective_config.yaml"))
```

Unknown keys and invalid types are errors. Environment-variable expansion is
not currently supported, and configuration contents should not be logged because
future files may contain secrets.

The optional reference config section records only the user-supplied FASTA path
and explicit build label. Phase 1.4 also provides the sample-sheet analysis
entry point:

```yaml
reference:
  fasta: /references/GRCh38.fa
  build: GRCh38
samples:
  sheet: /project/samples.tsv
```

All three values remain null in the packaged default, so foundation CLI commands
and the Phase 0 smoke DAG do not require analysis inputs. When an effective
`HiFiVarConfig` retains a user-config source, relative `reference.fasta` and
`samples.sheet` paths used by `AnalysisContext.from_config()` are interpreted
beside that user YAML. Relative paths without source provenance are rejected
rather than guessed from the process working directory. Paths inside a sample
sheet remain relative to the sheet itself.

## Reference genome

Phase 1.1 provides one immutable reference abstraction for future modules:

```python
from pathlib import Path

from hifivar.reference import ReferenceGenome

reference = ReferenceGenome.from_fasta(
    Path("/references/GRCh38.fa"),
    build="GRCh38",
)
```

The conventional `/references/GRCh38.fa.fai` must already exist. Contig names
and lengths come from that FAI in deterministic order; HiFiVar does not create
the index, infer the build, or convert names such as `1` and `chr1`. Primary
workflow references currently use uncompressed `.fa`, `.fasta`, or `.fna`.

Construction performs the lightweight FASTA check and reads the FAI, but does
not scan the entire FASTA. Request that only when provenance requires it:

```python
checksummed = reference.with_checksum()
# or ReferenceGenome.from_fasta(..., compute_checksum=True)
metadata = checksummed.to_dict(include_contigs=False)
```

The dictionary contains plain strings, integers, lists, and nulls suitable for
JSON/YAML. Full contig metadata is opt-in to avoid bloating ordinary summaries.

## Sample and primary input

Phase 1.2 represents one sample and its primary sequencing input through frozen
data models:

```python
from pathlib import Path

from hifivar.sample import InputDataset, Sample

dataset = InputDataset.from_files(
    [
        Path("/data/HG002.movie1.fastq.gz"),
        Path("/data/HG002.movie2.fastq.gz"),
    ]
)
sample = Sample(sample_id="HG002", input=dataset)
```

HiFi FASTQ input means one or more long-read files, such as separate movie or
delivery files. It is not modeled as Illumina-style R1/R2 paired-end data. BAM
and CRAM input each require exactly one file; mixed FASTQ/BAM/CRAM primary input
is invalid. BAM/CRAM indexes are optional at construction and can be required
later with `dataset.validate_index()`.

Machine `sample_id` values are deliberately ASCII-safe because they will be
used in filenames, workflow wildcards, result directories, and VCF metadata.
Allowed characters are letters, digits, underscore, hyphen, and dot; leading
dots, whitespace, separators, traversal, and Unicode IDs are rejected without
silent normalization. Filesystem paths remain Unicode-capable, so an ASCII
sample ID may refer to files inside directories such as `/data/测序数据/`.

`Sample.to_dict()` and `InputDataset.to_dict()` return only standard JSON/YAML
types. Paths are serialized as stored rather than resolved or converted to
absolute paths. Derived BAMs, VCFs, assemblies, QC metrics, checksums, pedigree,
and cohort metadata do not belong to these Phase 1.2 models.

Construction reuses the Phase 0 validation API. It checks only the first FASTQ
record and does not count reads, calculate N50/QV, compute coverage, checksum
large inputs, or scan an entire HiFi FASTQ. BAM/CRAM checks remain path-, suffix-,
non-empty-, and optional index-level only; headers and binary structure are not
parsed.

## Sample sheets and pedigree metadata

Phase 1.3 reads one strict UTF-8 or UTF-8-BOM TSV schema:

```text
sample_id	input	input_type	sex	father	mother	phenotype	group
HG002	data/HG002.fastq.gz	fastq	male	.	.	.	case
```

`sample_id` and `input` are required. `input_type`, `sex`, `father`, `mother`,
`phenotype`, and `group` are optional; column order is flexible, but canonical
lowercase names are required. Unknown and duplicate columns are errors. Blank
lines and lines beginning with `#` are ignored.

Load a sheet through the Python API:

```python
from hifivar.sample_sheet import SampleSheet

sheet = SampleSheet.from_tsv("project/samples.tsv")
record = sheet.get_record("HG002")
trios = sheet.get_trios()  # each tuple is child, father, mother
```

Relative input paths are interpreted relative to the sample-sheet directory,
not the process working directory. Absolute paths are retained. Multiple HiFi
FASTQ files use semicolons in input order:

```text
HG002	data/movie1.fastq.gz;data/movie2.fastq.gz
```

Empty semicolon components, mixed primary types, multiple BAM/CRAM files, reused
inputs across samples, duplicate sample IDs, and silently normalized sample IDs
are rejected. `Sample` and `InputDataset` remain the authoritative validation
models.

Declared sex accepts `M`, `F`, `male`, `female`, `unknown`, and `.`. Missing
metadata becomes `None`; no sex is inferred from sequence data. Parent IDs are
plain sample IDs resolved after all rows load. Parents must exist in the same
sheet, declared parent sex must be consistent when known, and pedigree cycles
are errors. Phenotype and group remain optional open text rather than a clinical
ontology.

Loading reads TSV metadata and performs the existing lightweight per-file checks
only. It does not scan complete FASTQs, parse BAM/CRAM headers, calculate input
checksums or coverage, validate against a reference, or run trio/cohort calling.

## Analysis context and run manifest

Phase 1.4 joins validated models without starting a workflow:

```python
from hifivar.context import AnalysisContext
from hifivar.manifest import RunManifest

context = AnalysisContext.from_config(config)
manifest = RunManifest.from_context(context)
manifest.write_json("provenance/run-manifest.json")
manifest.write_yaml("provenance/run-manifest.yaml")
```

The context requires one reference and at least one sample. Sample IDs and
primary input paths must be unique across the run, while different samples may
use different primary input types. Obvious `reference.fasta` and
`reference.build` conflicts are errors. Query contigs can be checked explicitly
with `context.validate_query_contigs(...)`; no chromosome names are rewritten.

FASTQ input is unaligned, so its reference status is `not_applicable` rather
than a compatibility claim. BAM/CRAM header/reference compatibility remains
`not_checked`, and the generic context does not require indexes. These checks
belong to later alignment-aware validation.

Manifest paths for reference, sample sheet, and primary inputs are absolute but
not symlink-resolved. The manifest retains standard JSON/YAML data only and
recursively redacts exact, case-insensitive `token`, `password`, `secret`,
`api_key`, and `apikey` keys. Input checksum calculation is disabled by default
for large HiFi files and is enabled only with
`RunManifest.from_context(..., compute_input_checksums=True)`. Reference
checksums are never calculated implicitly. Writers are atomic and require
`overwrite=True` to replace an existing manifest.

This layer does not implement a run-directory manager, workflow launcher,
Snakemake sample wildcards, QC, alignment, or variant analysis.

## Lightweight input QC

Phase 2.1 adds a Python-only QC framework on top of `AnalysisContext`:

```python
from hifivar.qc import run_input_qc

report = run_input_qc(context)
report.write_json("reports/input-qc.json")
report.write_yaml("reports/input-qc.yaml")
```

Validation and QC are intentionally separate. Missing, unreadable, empty, or
structurally unusable inputs raise the existing `InputValidationError` and stop
QC. Findings on inspectable data become `QCResult` statuses and issues. A
readable BAM/CRAM without a conventional BAI/CRAI is therefore WARN with
`ALIGNMENT_INDEX_MISSING`; it is not a generic input-validation failure and QC
does not create the index.

Current FASTQ metrics are limited to input type, ordered file paths, file count,
per-file/total byte size, and suffix-level `gzip`, `uncompressed`, or `mixed`
compression. Existing validation rechecks the first FASTQ record, but QC does
not scan all reads or calculate read count, N50, QV, yield, GC, or checksums.

Current BAM/CRAM metrics are limited to input type, path, byte size, and readable
conventional index presence. Headers, SM/RG values, contigs, binary integrity,
mapped/unmapped reads, duplicates, mapping quality, and coverage are not
inspected. No samtools, pbmm2, minimap2, pysam, or other external tool is used.

`run_input_dataset_qc()` supports one `InputDataset`; `run_input_qc()` preserves
the ordered single- or multi-sample context and returns a `RunQCReport`. Any
FAIL dominates WARN, WARN dominates PASS, and a report with no performed checks
is NOT_CHECKED. Reports contain version/UTC provenance and reuse the manifest's
shared atomic UTF-8 JSON/YAML writer. Existing files are protected unless
`overwrite=True` is explicit.

QC is read-only: it does not modify the context/input models, compute input or
reference checksums, create derived files, or start a workflow. Phase 2.1 has no
QC config section, CLI command, or Snakemake biological rule.

## Alignment planning interface

Phase 2.2 provides a side-effect-free handoff from a validated FASTQ context to
future alignment wrappers:

```python
from hifivar.alignment import AlignmentTool, build_alignment_requests

requests = build_alignment_requests(
    context,
    "results/alignment",
    tool=AlignmentTool.PBMM2,
    threads=24,
)
```

The returned frozen `AlignmentRequest` objects preserve sample and FASTQ order,
share the context reference, and use deterministic `{sample_id}.aligned.bam`
names by default. BAM/CRAM output is selected explicitly with
`AlignmentOutputFormat`. Paths are expanded but not resolved or translated, so
Linux/HPC path semantics remain under caller control.

All context samples must use FASTQ primary input. If any sample is already BAM
or CRAM, planning fails with `InputValidationError` and reports the incompatible
sample IDs and types; samples are never silently skipped. Output suffixes must
match their declared format, thread counts must be positive integers, and an
output cannot alias the input reference metadata paths.

`AlignmentBackend` is the minimal structural contract for a future wrapper: it
identifies its `AlignmentTool` and returns a shell-free `list[str]` from
`build_command()`. Phase 2.2 provides no concrete backend and does not create
directories, execute tools, inspect read groups, build indexes, or add CLI,
configuration, or workflow rules. Concrete pbmm2 work starts in Phase 2.3.

## Complete Phase 2 pipeline

Phase 2 now exposes the concrete Python orchestration boundary:

```python
from hifivar.phase2 import run_phase2

report = run_phase2(
    context,
    "results/alignment",
    dry_run=True,
)
report.write_json("results/phase2.json")
```

The `alignment` YAML section owns the selected tool/output format, threads,
memory, runtime, overwrite policy, pbmm2 preset/log level, index threads, and BAM
index format. Phase 2 execution currently permits only pbmm2 plus BAM. The
generic planner can represent minimap2 and CRAM output, but accepting those in
execution would falsely imply a wrapper exists.

FASTQ samples use `Pbmm2Wrapper`: validated FASTQ/reference input, CCS/HIFI
preset, stable RG/SM, on-the-fly coordinate sorting, and BAM output. Multiple
FASTQs use an ordered `fastq.fofn`. pbmm2's automatic BAM index is disabled so
the minimal `SamtoolsWrapper` owns one explicit `samtools index` step. Automatic
BAM index choice uses BAI unless a reference contig exceeds 2^29 bases, when CSI
is used. CRAM uses CRAI.

Existing BAM/CRAM inputs produce REUSE plans and keep their original paths.
They never invoke pbmm2. Their sort order is UNKNOWN because Phase 2 does not
parse headers; consequently, a missing index produces a warning and is not
silently rebuilt. Users must explicitly establish coordinate sorting before a
future indexing request.

`run_phase2()` performs input QC, planning, pbmm2 or reuse, output validation,
indexing where safe, alignment QC, and an atomic JSON/YAML `Phase2RunReport`.
Dry-run previews both pbmm2 and samtools commands without requiring either tool
or creating outputs. Real tool execution is a Linux/HPC verification task.

Alignment QC is intentionally lightweight. It records file size, BAM/CRAM
format, generated/existing source, declared sort order, and index presence. It
does not parse alignment headers or records, calculate coverage or mapping
statistics, call pysam, or run a variant caller.

## DeepVariant small-variant pipeline

Phase 3 consumes aligned, indexed BAM/CRAM—not raw FASTQ—and retains SNV/Indel
outputs separately:

```python
from hifivar.phase3 import collect_phase2_alignment_artifacts, run_phase3

phase3 = run_phase3(
    context,
    "results",
    alignment_artifacts=collect_phase2_alignment_artifacts(phase2),
    dry_run=True,
)
phase3.write_json("results/phase3.json")
```

Outputs are `results/small/{sample}.small.vcf.gz` and
`results/small/{sample}.g.vcf.gz`, each with a `.tbi`. The `small` YAML section
selects native, Docker, or Apptainer execution, an explicit image for container
modes, PACBIO model, threads, memory, runtime, and overwrite policy. Container
command construction remains isolated in `DeepVariantRuntime`; all execution
passes through `DeepVariantWrapper` and `CommandRunner`.

Output validation checks BGZF framing, VCF header, exact sample, reference
contig names, gVCF `NON_REF`, and TBI headers without loading variant records.
BAM/CRAM header compatibility is a declared Linux/external-tool assumption, not
silently inferred. Phase 3 does not normalize alleles, perform joint calling, or
mix SV/TR calls into the small-variant VCF.

## External commands

All external tools must be executed through `CommandRunner`; wrappers must not
implement their own subprocess handling:

```python
from hifivar.command import CommandRunner

runner = CommandRunner()
result = runner.run(["program", "--version"])
```

Execution uses an argument list with `shell=False`. Use `redact_values` for
sensitive arguments and `stdout_path`/`stderr_path` for output that should be
streamed to files instead of retained in memory. Captured output is decoded as
UTF-8 with invalid bytes replaced; redirected files receive the tool's original
bytes without text transcoding.

## CLI

The argparse CLI loads packaged defaults, then preset and user YAML, applies
explicit CLI overrides, and finally configures logging. Global options must be
placed before the subcommand:

```bash
hifivar --preset standard --log-level DEBUG config show
```

CLI modules dispatch application actions but must not contain tool-specific
execution logic. Future external commands remain behind wrappers and
`CommandRunner`.

The single source of default and preset YAML is
`src/hifivar/resources/configs/`; these files are installed as package data.

## Validation

Validation never silently fixes biological inputs. It reports hard failures
and leaves indexing, sorting, normalization, and contig-name changes to explicit
future workflow steps.

The Phase 0.8 API provides:

- readable file and directory checks, including optional non-empty checks;
- lightweight streaming checks for FASTA, FASTQ, VCF/VCF.GZ, and BED;
- basic FAI parsing with unique contig names;
- BAM/CRAM path, suffix, and optional index-presence checks;
- TBI/CSI presence checks for compressed VCF input;
- query-contig subset checks without automatic `chr` prefix conversion;
- non-empty expected-output checks; and
- chunked SHA256 checksums.

BAM/CRAM binary structure is not parsed in Phase 0.8. Reading a gzip VCF does
not prove that it is BGZF or tabix-compatible, and index files are checked only
for readable path presence. Deeper checks require future tool-aware validation.

## Snakemake workflow

Install workflow tooling separately from the lightweight package runtime:

```bash
python -m pip install -e ".[dev,workflow]"
```

HiFiVar resolves and validates default, preset, user, and CLI configuration
before Snakemake starts. Snakemake consumes only the resulting effective YAML;
it does not repeat deep merge or the full Python schema validation:

```text
HiFiVar Config API
    ↓
effective_config.yaml
    ↓
workflow/Snakefile
```

Create that handoff and exercise the default Phase 0 smoke target with:

```bash
hifivar --preset standard config dump-effective --output effective_config.yaml
snakemake \
    --snakefile workflow/Snakefile \
    --configfile effective_config.yaml \
    --cores 1 \
    --dry-run
snakemake \
    --snakefile workflow/Snakefile \
    --configfile effective_config.yaml \
    --cores 1
```

`workflow/Snakefile` defines the default target and includes shared conventions
from `workflow/rules/common.smk` plus executable smoke rules from
`workflow/rules/smoke.smk` and the disabled-by-default Phase 3 DeepVariant rule
from `workflow/rules/small.smk`. Null `paths.workdir` and `paths.outdir` values map to
relative `work/` and `results/`; configured paths override them. Rule logs use
`logs/<module>/<rule>.log`.

Rule names use `snake_case`. Foundation resources use `threads`, `mem_mb`, and
`runtime_min`; per-tool resource blocks will be introduced only with real tool
rules. Re-creatable intermediates may later use `temp()`, while important final
outputs may use `protected()` where appropriate.

Snakemake's `benchmark` directive records runtime/performance metadata and is
distinct from HiFiVar's future scientific variant-accuracy benchmark module.
DeepVariant runtime selection is config-driven through its wrapper rather than
duplicated in Snakemake. Scheduler profiles remain future Linux/HPC work. The
Phase 0 Python-only smoke DAG and Phase 3 dry-run are cross-platform, but Linux
remains the production workflow platform.
