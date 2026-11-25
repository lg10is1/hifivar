# Capabilities, execution and evidence / 功能、执行与证据

[简体中文](zh_CN/capabilities.md)

Applies to release source 0.1.0rc5. This page describes **RC5 release source**, not an upgrade of RC1–RC4 or every package reporting 0.1.0rc5. Select a
package by its trusted source/artifact manifest; see [candidate identity](rc5_candidate.md).
Implemented routing, successful execution and matched-truth accuracy are separate claims.

本页适用于 RC5 发行源码。功能存在、真实运行成功、匹配真值准确率是三个不同层次；
版本号相同不等于候选包内容相同。历史发行资产保持不变。

## What the launcher connects / 主入口能够衔接什么

The calling launcher accepts validated indexed BAM inputs. Opt-in alignment first
prepares BAMs from supported FASTQ/native uBAM. Explicit source plans then choose
new same-attempt outputs (GENERATED) or validated existing artifacts (EXTERNAL)
for downstream consumers.
It does not infer missing data, secretly rename samples/contigs or silently reuse
an older calling receipt. Requested branches keep separate artifact families.

| Layer / 层次 | Implementation and route / 实现与入口 | Evidence and boundary / 证据与边界 |
|---|---|---|
| Input/reference / 输入参考 | Sample sheets, indexes, exact RG/SM and contig/length checks, hashes and run identity | Header compatibility is not biological sample authentication or complete REF sequence identity |
| Alignment / 比对 | pbmm2/minimap2; opt-in sequential raw-input preparation before calling | Existing indexed BAM validation/reuse is not a new alignment; main-launcher CRAM remains deferred |
| SNV/Indel | DeepVariant; native VCF and gVCF kept separate | Regional and separate HG002 whole-genome evaluation exist; one bounded example is not a fresh WGS study |
| Read-SV | Sawfish, Sniffles2, pbsv, cuteSV; enable callers individually | Native caller evidence and matched-truth comparisons are distinct; retain historical identity qualifications |
| Tandem repeats / 串联重复 | TRGT with pinned reference-compatible catalog | TR is independent of SV/Indel; no general TR accuracy claim without independent matched truth |
| Assembly / 组装 | Explicit HiFi FASTQ to hifiasm; independent assembly context | Reused FASTAs do not acquire a fabricated new hifiasm receipt |
| Assembly-SV | PAV/SVIM-asm from generated or explicitly declared ready FASTAs | PAV mixed root remains immutable; derived SV-only selection is explicit. Latest ready-PAV native acceptance remains incomplete |
| Phasing / 定相 | HiPhase small-only route plus opt-in pbsv/Sawfish and/or TRGT joint sources; separate phased outputs | Identified installed-source regional joint cases passed independently, including a two-technical-label batch; not all-package/WGS or switch/flip accuracy |
| Cohort / 多样本 | GLnexus consumes real gVCFs; separate SV/TR matrices retain source and missing states | No cross-sample SV clustering or joint TR genotyping is implied; TR catalogs must match |
| SV harmonization / 整合 | Per-sample Jasmine with source evidence and Truvari comparisons; explicit ready/generated sources | Caller support is evidence metadata, not truth/confidence. BND representation is not an accuracy result |
| Annotation / 注释 | Independent ANNOVAR/VEP outputs and explicit SV/TR region overlap | Native scoped executions exist; a successful backend case does not clear every later route/deployment failure |
| Functional prediction / 功能预测 | Explicit selected-candidate AlphaGenome v1 CLI using the user's own API environment credential | Optional pinned SDK; local REF/context check and network consent. Offline contract acceptance is separate from authenticated service prediction; not call confidence or clinical truth |
| Review / 复核 | Offline review preparation/import, plus opt-in generated/external IGV consumer | Technical screenshots remain NOT_REVIEWED until explicit independent decisions; never overwrite caller VCFs |
| Benchmark / 基准评估 | Explicit query/truth plans; hap.py, Truvari, same-catalog TR comparison | One declared matched-truth sample per benchmark invocation; regional engineering evidence is not universal accuracy |
| Delivery / 交付 | Per-category outputs, logs, versions, config, input/source identities, attempt receipts and final report | Reporting is a separate explicit CLI; missing/failed evidence remains incomplete |

## One invocation, staged execution / 一次调用、分阶段执行

Configure only needed branches and provision their tools/resources beforehand.
The launcher uses the following enabled stages, not one monolithic all-tool DAG:

1. Optional raw-input alignment preparation.
2. Producer DAG: selected calling, assembly/assembly-SV, phasing and within-sample harmonization.
3. Optional annotation/cohort consumer DAG with explicit source plans.
4. Optional benchmark consumer DAG with matched truth/BED and native query roles.
5. Optional IGV consumer DAG with explicit candidate selection and actual BAM binding.
6. A **separate** report command after selecting the exact attempt.

Disabled stages are not run. This is not distributed Slurm execution, automatic
all-tool installation, automatic candidate interpretation or a claim that every
branch has been validated together at whole-genome scale. Small VCF, gVCF, SV,
TR, phased VCF, annotation, benchmark and review evidence remain distinct.

主入口按配置串接多个阶段，不自动打开所有工具。它不会把旧结果的存在当成本次运行
成功，也不把截图、来源支持数或功能影响当作真值。所有分轨原始结果保留。

## Choose the relevant guide / 选择运行说明

- [Indexed BAM batch execution / 已比对BAM批量运行](indexed_bam_run.md)
- [Complete small/SV/TR batch template / 完整批量多类型模板](batch_multivariant_run.md)
- [FASTQ/native-uBAM preparation / 原始输入准备](raw_input_run.md)
- [Ready haplotype FASTAs / 已有单倍型组装](ready_assemblies.md)
- [Within-sample SV integration / 同样本SV整合](sv_harmonization_run.md)
- [Optional joint small/SV/TR phasing / 可选多类别定相](joint_phasing.md)
- [Explicit annotation/cohort source plan / 显式下游来源](indexed_bam_run.md#explicit-downstream-handoff--显式下游衔接)
- [Owned benchmark execution / 基准评估衔接](benchmark_run.md)
- [Review execution / 复核衔接](review_sources.md) and [independent review decisions / 独立人工判定](review_cli.md)
- [Deployment / 依赖部署](deployment.md) and [preflight / 运行前检查](preflight.md)

Keep all runs under an explicitly chosen project root; the example paths below
are placeholders. Inspect a new plan, execute matching frozen controls, then
choose the actual immediate-child attempt ID recorded by that invocation:

~~~bash
hifivar --config analysis.yaml --dry-run run --run-dir /your/project/runs/run01
hifivar --config analysis.yaml run --run-dir /your/project/runs/run01 --resume
hifivar report --run-dir /your/project/runs/run01 --attempt-id attempt_ACTUAL_ID --output-dir /your/project/reports/run01
~~~

An ID placeholder must not be copied literally. Dry-run still reads inputs and
writes metadata; it does not validate future generated calls. Resume is not an
automatic repair/retry mechanism. Changed inputs/config/runtime need a new run;
retain failure logs and verify that old processes have terminated first.

## Evidence that must not be pooled / 不可混淆的证据

- An identified later installed-source case completed one analysis invocation
  with fresh small variants, all four read-SV callers, TRGT, within-sample
  Jasmine, joint small/SV/TR HiPhase and all three cohort tracks. Its two HG002
  labels are technical partitions (biological N=1), using regional BAMs and one
  catalog locus. The final report remained separate. This is not whole-genome,
  independent-cohort, TR/phasing accuracy or every-later-package acceptance.
- Historical two-sample regional calling and harmonization runs are separate
  executions, not one combined all-module experiment. Technical labels derived
  from one HG002 dataset are not independent biological samples.
- A later installed-candidate HG002 case completed fresh DeepVariant, Sniffles2,
  TRGT and HiPhase producers, then same-attempt small/SV benchmark and TR technical
  review. Query/index hashes and native producer receipts were independently
  checked. This is bounded orchestration evidence, not WGS throughput or a
  substitute for the fixed scientific benchmark protocol.
- Separate existing-query small/read-SV/assembly-SV comparisons establish the
  EXTERNAL route, not new calling. Generated TR review establishes source linkage,
  not independent TR truth accuracy.
- The small-variant WGS study, CMRG medically relevant gene evaluation and other
  draft-truth SV research retain their exact regions, filters, denominators and
  software identities. Do not merge these numbers or transfer old validation to
  a later same-version package. See [benchmarking](benchmarking.md).

已运行不等于准确率已证明。低深度标称值不是观测深度；嵌套技术抽样不是生物学重复；
没有普适最低深度、测序成本降低、临床效用或所有SV类型准确率结论。人工判读、作者
信息、数据许可与发布批准由相应责任人提供。未完成的后端/入口保持未完成状态，不能
用别的案例的绿色检查替代。

## New dedicated research entries / 新增独立研究入口

[SV depth-response study](sv_depth_study.md) runs per-depth/replicate read-SV calling and fixed-protocol comparison from existing BAMs. The historical small-variant depth curves are not SV results. [AlphaGenome](alphagenome.md) is an opt-in selected-variant API command, not an automatic main calling/annotation stage. Both belong to new RC5 source and require their exact new candidate identity; no RC1–RC4 asset is changed.
