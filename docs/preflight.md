# Calling preflight: check before a fresh batch run

[简体中文](zh_CN/preflight.md)

RC5 includes this preflight interface. It is absent from RC1–RC4
and the earlier RC5 candidates, including the review-only build. Its bounded
Linux gate has passed for the
[identified 15 October payload](rc5_candidate.md#validated-development-payload--已验收开发包),
which includes this extension.
Version alone does not identify the payload: compare its source manifest and
package hash, including after documentation revisions.

## Two levels, explicit execution

The existing `hifivar doctor` remains a foundation-environment summary. It does
not certify a caller deployment. The opt-in entry checks the **fresh calling
launcher** scope using its actual configuration/input/resource contracts:

```bash
hifivar --config analysis.yaml doctor --preflight --run-dir /your/project/run01 --json
```

The directory must not yet exist; it is inspected, **not created**. Relative
reference/sample/executable/catalog/image paths use the existing launcher's
documented resolution rules. No repair, index building, alignment, calling,
workflow execution, installation, download or automatic retry occurs. Unlike a
calling dry-run, metadata-only preflight does not run Snakemake at all. File
logging is not opened for this opt-in entry; redirect stdout yourself only to an
authorized project location if you want to retain its report.

Metadata checks include enabled-track scope, request-versus-configured scheduler
budgets, readable FASTA/FAI and explicit build, lightweight sample inputs,
conventional unique BAM index where required, workflow resources, enabled tool
lookup, and TR catalog structure/contigs/bounds/declared build/karyotype.
Disabled callers do not become mandatory dependencies. Alignment tools are
checked when alignment.enabled, even if some inputs could later be reused.
Opt-in within-sample read-SV harmonization additionally requires Jasmine,
bgzip/tabix and Truvari for the existing Phase 9 source comparisons. These
dependencies are not probed when harmonization is disabled.
With probes enabled, this launcher route requires tabix at least 1.23.1 for
Jasmine VCFv4.4 output. An older indexer is BLOCKED before analysis; preflight
does not upgrade it or rewrite a VCF header. This is a conservative deployment
floor for this route, not a change to independent VCFv4.3 workflow behavior.

To additionally allow **bounded native version and BAM metadata probes**, opt in:

```bash
hifivar --config analysis.yaml doctor --preflight --run-dir /your/project/run01 --probe-tools --timeout-seconds 30 --json
```

Bioinformatics version commands remain owned by existing wrappers and execute
through CommandRunner without a shell. Snakemake is probed in the isolated active
Python interpreter. Indexed BAM probes reuse the launcher validator: EOF check,
coordinate-sort declaration, every RG/SM, exact reference names/lengths and one
indexed region. They do **not** scan every alignment or certify read quality.
No probe log file or output directory is created by HiFiVar. Native programs may
manage their own caches; this is not an OS sandbox. Provision a trusted tool
environment and existing task-isolated temporary/cache locations before probing.
Credential/Python injection variables are not inherited by these probes.

Global options (`--config`, `--dry-run`) precede `doctor`; doctor-specific flags
follow it. `--dry-run` plus `--probe-tools` is refused rather than silently
executing tools. `--probe-tools`, `--json` and `--run-dir` require `--preflight`.

## Report states and exit codes

| State | Exit | Meaning |
| --- | --- | --- |
| CHECKED_SCOPE_PASS | 0 | Required checks in this bounded scope passed |
| BLOCKED | 1 | At least one failed check; fix the declared input/environment/configuration explicitly |
| INCOMPLETE | 2 | No failure, but at least one required check remains NOT_CHECKED |

Exit 2 is an intentional incomplete result, not a failed analysis. Metadata-only
checks normally remain INCOMPLETE because tool versions and BAM identities have
not been probed. One invalid launcher contract stops deeper probes; no empty
phase0 run is misrepresented as a valid batch analysis.

The JSON records each check, configured/resolved executables, scoped versions,
commands, exit codes and durations for probes. It never supplies an accuracy
metric. Resource observations are advisory: free filesystem bytes, inode counts
where supported, CPU count and OS process/file-descriptor limits where available.
They do not measure quota, per-UID remaining thread headroom, allocated RAM or
the storage required by a new WGS run. These remain explicitly NOT_CHECKED or
NOT_ESTIMATED. No caller/shard/resource setting is silently reduced.

## Deliberate incomplete/unsupported boundaries

- DeepVariant Apptainer images are checked for readable path metadata only;
  image digest/model/runtime provenance is not certified. Containers are never
  started or pulled, including Docker images. Their version check remains
  NOT_CHECKED even when --probe-tools is supplied.
- Native uBAM headers may be classified for explicit pbmm2 routing, but the
  expensive flag/RG record scan is deliberately not run here. uBAM remains
  INCOMPLETE; the actual raw launcher retains its full preflight contract.
- FASTQ validation reads a prefix, not all molecules or records. Assemblies are
  not produced, and no minimum-depth or sequencing-cost conclusion is made.
- CRAM/benchmark/review/low-depth routes remain subject to the launcher boundary.
  Later source supports explicit ready SV inputs and some assembly/annotation/cohort
  routes; this bounded probe registry does not certify every backend (see below).
  Harmonization stays per sample; no cross-sample SV clustering is introduced.
- This is fresh-run checking only, not resume identity, an HPC scheduler dry-run,
  Slurm resource validation or a replacement for actual tiny Linux acceptance.
- A scoped pass does not guarantee eventual workflow success, all-workflow
  readiness, raw-data immutability, accuracy, clinical suitability or human review.

Reports can contain sample names, paths and comments inherited from config errors.
Treat them as private runtime evidence; do not upload them automatically. Do not
embed credentials in any path/free text. Existing published packages and their
historical acceptance records remain unchanged.

## Later development-source corrections

Harmonization dependency checks follow `sv.harmonization.enabled` independently
of `sv.enabled`. Explicit ready VCFs therefore still require Jasmine, Truvari
and their configured BGZF/TBI tools even when every read caller is disabled.
A FASTQ used only as metadata context does not require an unused BAM checker;
actual alignment preparation or BAM metadata probes retain that requirement.
When source identity is enabled, bgzip also shares its wrapper's >=1.23.1 floor.
Presence without opt-in version probes remains INCOMPLETE; no analysis is run.

The later development-source registry additionally checks the enabled native
backends below. All commands remain owned by the existing wrappers; no analysis
command is constructed or run by preflight.

| Enabled branch | Version-only probe scope |
| --- | --- |
| assembly / hifiasm | hifiasm |
| assembly_sv / SVIM-asm | SVIM-asm, minimap2, samtools, bgzip, tabix |
| annotation / VEP | VEP program, not full cache contents |
| cohort / small variants | GLnexus and bcftools, not joint genotyping |

Parent and child switches both apply. A missing dependency blocks its group;
a version failure is reported separately and does not suppress observations
for other groups. Without --probe-tools, executable presence never substitutes
for version evidence. Successful observations record version_source=command
and a version-only validation scope, not database or result validation.

PAV and ANNOVAR still retain required NOT_CHECKED probe-scope observations.
PAV's configured source version is not the Snakemake adapter version; ANNOVAR's
declared source date is not an observed native version command. Neither is
silently promoted to a native PASS. No container starts, cache downloads,
assemblies, annotation jobs or GLnexus analyses occur. Actual configured input,
resource and launcher validation is still performed before probes. A scoped
PASS remains distinct from all-workflow readiness and accuracy.

These additions are in the GitHub RC5 release. The historical development builds used separate manifests. Historical package/native gates do not
certify a changed payload; compare its manifest and build hash.
