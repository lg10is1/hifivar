# Phase 9 - SV Harmonization / Integration

Status: **COMPLETE**

Phase 9 integrates independent read-derived and assembly-derived structural-variant
evidence. It does not determine truth, assign clinical meaning, or convert caller
support counts into confidence.

## Completed scope

| Subphase | Result |
| --- | --- |
| 9.1 Unified evidence model | Read/assembly sources, caller, original VCF and record ID, coordinates, SV type/length, native INFO, uncertainty, haplotype and run status are represented. |
| 9.2 Normalization boundary | VCF records are streamed; only unambiguous fields are extracted. Unknown complex types remain UNRESOLVED. |
| 9.3 Jasmine | Dedicated wrapper supports executable/version checks, deterministic arguments, input list, dry-run, logs, output validation and provenance. |
| 9.4 Truvari | Minimal bench comparison wrapper records summary and explicitly labels results comparison-only. |
| 9.5 SURVIVOR | Intentionally not introduced; Jasmine remains the primary harmonizer. |
| 9.6 Multi-caller integration | Sawfish, Sniffles2, pbsv, cuteSV, PAV and SVIM-asm artifacts remain intact and feed a per-sample harmonized VCF plus evidence TSV. |
| 9.7 Evidence classes | READ_ONLY, ASSEMBLY_ONLY, READ_AND_ASSEMBLY, and conservative UNRESOLVED are source categories only. |
| 9.8 Haplotype preservation | Assembly haplotype roles are retained in source provenance and the evidence table. |
| 9.9 BND/INS/complex boundary | Native BND and caller INFO are preserved; explicit INS length is recorded; ambiguous complex records are not forced into simple classes. |
| 9.10 Integration API | run_phase9() executes Jasmine, per-source Truvari comparison, and run provenance. |
| 9.11 Snakemake | Modular config-driven harmonize_sv rule consumes explicit enabled caller outputs. |
| 9.12 Partial evidence | NO_CALLS, NOT_RUN, FAILED, and DISABLED are distinct. Missing or failed callers are never interpreted as zero variants. |
| 9.13 Validation | BGZF/index/header/sample/reference validation and source traceability are required. |
| 9.14 Integration tests | Fake six-caller E2E and read-only Snakemake dry-run are covered. |
| 9.15 Scientific safeguards | Tests prohibit truth/confidence inference and protect original/native evidence. |

## Linux audit remediation

Jasmine 1.1.5 is now invoked through `bash` when the resolved executable is a
text launcher without a shebang. Source BGZF VCFs are stream-decompressed in
`work_directory` in the exact runnable-source order required by SUPP_VEC;
original caller VCFs are never changed. Jasmine's raw VCF is retained, missing
FORMAT definitions are copied only from actual source headers, and records are
external-sorted by declared contig order and POS before bgzip/tabix.

Truvari version detection now uses `truvari version`, with `--version` retained
only as an older-version fallback. These changes address P9-JASMINE-001/002/003
and P9-TRUVARI-001. Linux delta execution is still required before changing
the independent audit result to PASS.

Post-RC2 real-data validation identified that native IDs from Sawfish,
Sniffles2, pbsv, and cuteSV do not share a caller-prefix convention. Evidence
membership therefore uses Jasmine `SUPP_VEC` as the authoritative mapping to
the ordered runnable-source list whenever that field is present. `IDLIST`
remains unchanged as source-record provenance and is never used to override a
valid support vector. This corrects derived support counts without rewriting
the harmonized VCF or any caller VCF.

## Public API

- EvidenceRunStatus, EvidenceClass
- SVEvidenceSourceArtifact, SVEvidenceCollection, SVEvidenceRecord
- SVHarmonizationRequest, HarmonizedSvArtifact
- iter_sv_evidence(), write_evidence_table(), validate_harmonized_artifact()
- JasmineWrapper, JasmineCommandPlan, JasmineResult
- TruvariWrapper, TruvariRequest, TruvariResult
- Phase9Settings, Phase9RunReport, run_phase9()

## Deterministic outputs

- results/sv_harmonized/{sample}/{sample}.harmonized.sv.vcf.gz
- results/sv_harmonized/{sample}/{sample}.harmonized.sv.vcf.gz.tbi
- results/sv_harmonized/{sample}/{sample}.sv.evidence.tsv
- results/sv_harmonized/{sample}/truvari/{caller}/summary.json
- results/sv_harmonized/{sample}/{sample}.phase9.provenance.json

All raw caller VCFs and indexes remain separate and are not rewritten.

## Configuration

sv.harmonization contains only merge/comparison settings: backend, tool
executables, resources, maximum distance, distance mode and overwrite policy.
Caller-specific settings remain in their own read or assembly caller sections.
Unknown configuration keys continue to fail schema validation.

## Scientific boundary

- Caller count is metadata, not correctness.
- READ_AND_ASSEMBLY is an evidence-source class, not truth.
- Jasmine clustering does not create a benchmark truth set.
- Truvari output is concordance/comparison support only.
- No benchmark confidence model, annotation, review UI, or Phase 10 behavior is included.
- Full VCFs are processed with generators/streaming; pandas is not used.
- Harmonization remains per sample. Phase 12 cohort SV tables preserve
  source-native rows and do not invent cross-sample clustering; adding a
  cross-sample SV clustering contract requires a separate scientific design and
  benchmark, not an implicit extension of this phase.

## Verification status

Windows verification uses fake command runners, synthetic BGZF/TBI artifacts,
unit tests, fake integration tests, and Snakemake dry-run.

- Jasmine official command contract: reviewed against the upstream repository.
- Truvari bench contract: reviewed against official documentation.
- LINUX_REVALIDATION_REQUIRED: P9-JASMINE-001, P9-JASMINE-002,
  P9-JASMINE-003, P9-TRUVARI-001

The supported Linux/HPC deployment boundary is documented in
`docs/deployment.md`.

## Known limitations

- Jasmine SUPP_VEC is the authoritative deterministic source mapping when
  present. Native IDLIST values remain provenance because caller ID formats are
  heterogeneous; legacy records without either mapping are marked UNRESOLVED.
- Biological clustering thresholds require later benchmark calibration.
- PAV, Jasmine, and Truvari remediation passed independent Linux/HPC delta
  validation; HiPhase, hifiasm, SVIM-asm, and upstream-tool results remain
  independent validation evidence.
