# Explicit source review — development execution bridge

[简体中文](zh_CN/review_sources.md)

This RC5 source adds an opt-in final review consumer to 'hifivar run'.
It does not change historical release assets. The source-binding library remains
read-only; the separate execution bridge runs IGV through the existing wrapper,
seals actual command/output receipts and keeps human review NOT_REVIEWED.
This is orchestration, not a new caller, truth benchmark or automated adjudicator.

## 1. Select explicit sources and records

A dedicated schema-1 YAML chooses a source and nonempty candidate list:

    schema_version: 1
    review:
      - sample_id: HG002
        track: small
        origin: GENERATED
        producer: small.deepvariant
        output_role: vcf
        candidates:
          - review_id: example_snv
            variant_type: SNV
            site: {contig: chr21, position: 20000052, ref: G, alt: A}

This coordinate is syntax only, not a promised call. GENERATED resolves only the
selected owned attempt's enabled producer and explicitly selected VCF role.
Future filenames, another attempt, an old file, a version declaration or a dry-run
cannot substitute for a valid native receipt. Small review consumes a VCF, not gVCF.

EXTERNAL rows instead declare existing path, index_path, source_tool,
source_version and reference_build. Versions remain user-declared and unprobed.
Paths resolve beside the source plan. BGZF, conventional TBI, matching sample/
build/contig identity and exact selected records are required. External artifacts
and complete BND bundles are frozen in the launch input identity; they are not
relabelled as new calls. The general annotation/cohort plan parser still rejects
this separate review-only document.

## 2. Enable the final consumer

Add to the normal analysis config, retaining the required reference, sample sheet,
selected producers and adequate thread/memory budgets:

    review:
      enabled: true
      source_plan: review_sources.yaml
      selection_file: null
      igv_executable: /your/project/tools/igv.sh
      threads: 1
      memory_mb: 8000
      runtime_minutes: 240
      overwrite: false

The executable needs an appropriate graphical/headless deployment and Java.
HiFiVar does not install/start a display server. All caller and display software
must already be available in the configured Linux environment.

'source_plan' and 'selection_file' are mutually exclusive. The launcher supports
source_plan mode; selection_file remains the explicit standalone workflow/offline
review mode. Disabled defaults are unchanged. No automatic candidate ranking,
confidence threshold, chromosome renaming or VCF rewriting is introduced.

    hifivar --config analysis.yaml --dry-run run --run-dir /your/project/runs/run01
    hifivar --config analysis.yaml run --run-dir /your/project/runs/run02

A fresh run directory is required. Dry-run performs the existing read-only native
BAM preflight and previews producers; it does not run IGV or certify future source
records. Generated downstream and raw-alignment plan-only statuses retain their
own scope. Resume requires matching frozen controls, inputs and runtime identity;
there is no stale-lock removal, automatic retry or old-receipt inference.

## 3. Execution and alignment provenance

The launcher freezes the source intent before output creation, executes producers
with review disabled, then any existing annotation/cohort consumer. Only afterward
it creates a new attempt-local selection TSV, record-binding JSON and review-only
configuration. The final review DAG disables calling, assembly, phasing,
harmonization, annotation and cohort; it does not rerun those producers.

Review uses the actual prepared BAM context. For externally aligned BAMs it binds
this attempt's samtools header/EOF/index probes, declared RG/SM and exact reference
dictionary to frozen BAM/index identities. For raw-input preparation it verifies
the sealed same-attempt alignment handoff instead. Neither mode authenticates a
biological sample or proves full-record/alignment correctness. CRAM remains deferred.

Missing/failed producers or changed source/control/alignment bytes stop review.
A scheduler success code without the selected native review receipt is failure.
CommandRunner observes version/batch execution, return codes, duration and logs.
Required manifests, batch, provenance and every expected endpoint PNG are sealed
with their output hashes. PNG signature/dimensions are engineering checks, not
full rendering QA or human review. A BND with unknown/virtual endpoints keeps
those explicit; no image is invented for a nonrenderable endpoint.

## 4. Selection semantics

- variant_id selects a unique actual raw VCF ID. Missing or duplicated IDs fail.
- site selects exact CHROM/POS/full REF/full ALT, including raw ID '.'. The internal
  SITE label never rewrites the VCF. Multiallelic records retain all alleles;
  selection is record-level, not automatic allele normalization or allele-specific review.
- Choose ID OR SITE, never both. Ambiguous/duplicate source-record selection fails;
  no fuzzy positional fallback occurs.
- TR uses explicit TRID:<locus>, declared TRGT origin and catalog SHA256. Its repeat
  interval excludes the VCF left anchor. It remains an independent category.
- BND uses actual record ID and ALT anchors, not END as an event span. Single/
  virtual endpoints and unresolved MATEID remain explicit. Derived BND review
  requires its complete reversible native/ledger bundle; detached SITE selection
  cannot bypass that source binding.
- SNV/INDEL follows explicit sequence allele lengths. MNV-only, symbolic reference
  blocks and star alleles are not coerced. INS uses its anchor; DEL/DUP/INV requires
  a matching declared SVTYPE and nonempty END interval. These are not REF-sequence,
  genotype, biological-validity or accuracy checks.

## 5. Outputs and report

The selected attempt owns review/effective.yaml, selection.tsv,
source_bindings.json and materialization.json. Results under results/review contain
review.igv.batch, screenshots, JSON/YAML/TSV pending-review manifests and
review_provenance.json. The native step receipt binds the exact source/control/
alignment chain and all expected outputs. Failed attempts remain failed.

Reporting remains an explicit separate command, not an automatic part of run:

    hifivar report --run-dir /your/project/runs/run02 --attempt-id attempt_EXACT_ID --output-dir /your/project/reports/run02

Choose the real immediate-child attempt ID; do not guess one from mtimes. Reports
separate technical IGV completion from NOT_REVIEWED human status and
NOT_ASSESSED accuracy. Use the independent offline review/decision importer for
explicit reviewer records; no decision is inferred from screenshot availability,
caller support, annotation consequence or a pipeline exit code.

## 中文说明

开发源码新增“主入口→调用→可选注释/cohort→最终review consumer”的显式衔接，
不是把所有步骤塞进一个DAG，也没有改动历史发布包。配置 review.source_plan；
与 selection_file 互斥，后者仍用于独立工作流/离线复核。来源计划必须明确选定
producer/VCF角色或现成文件，以及原始ID、精确位点或TRID候选。新结果必须使用
同一attempt的有效执行回执，不能借用旧结果。多等位、TR和BND来源语义保持独立。

实际review使用本次准备好的BAM：既有BAM绑定原生头部/EOF/索引预检及输入哈希；
原始读段路线绑定同次alignment handoff，而不是把FASTQ当成IGV输入。review阶段
关闭所有producer及其他consumer，只有技术截图和独立证据输出。需要预装IGV、
Java和可用显示环境；软件不自动部署显示服务。退出码0但没有有效回执也会失败。

最终报告仍需单独调用report命令，并明确选择attempt。截图完成不代表人工判定完成，
默认NOT_REVIEWED，也不是新的准确率、临床意义或调用可信度结论。源码库本身的
resolve_review_targets仍只解析来源、不执行IGV。正式验收以精确候选包、运行和
测试范围为准，不根据共享rc版本号追认其他运行。已冻结benchmark不因此改变。
