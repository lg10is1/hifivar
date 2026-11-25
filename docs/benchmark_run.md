# Owned benchmark execution / 主入口中的显式 benchmark

[简体中文](zh_CN/benchmark_run.md)

## Availability / 可用范围

This is an **RC5 interface**, not a feature of RC1–RC4 or
every package labelled `0.1.0rc5`. Match the complete source and package SHA256
manifest. The previous standalone Phase 13 interface remains available.

这是 RC5 提供的接口，必须核对源码和安装包清单，不能仅凭版本字符串确认功能。
默认关闭，不会改变 caller、过滤或匹配参数，也不会改写历史 benchmark 结果。

`hifivar run` can now execute an opt-in benchmark-only consumer DAG after its
producer and annotation/cohort stages, before optional review. Reporting remains
an explicit, separate `hifivar report` invocation. The consumer runs the existing
hap.py (small variants), Truvari (read/assembly SV), or exact same-catalog TR
comparison. It does not add a caller or an accuracy-optimization algorithm.

主入口可选择衔接 small、read-SV、assembly-SV、TR 的既有比较器。四类产物保持分开；
一次 benchmark 明确选定样本表中的一个 truth 样本，而不是自动为全部样本匹配真值。
多个样本的 calling 与“每个样本均有匹配真值的 benchmark”是不同能力。

## Explicit source and truth / 显式指定查询与真值

Example source plan (`benchmark_sources.yaml`):

```yaml
schema_version: 1
benchmark:
  sv:
    sample_id: HG002
    track: sv
    origin: GENERATED
    producer: sv.sawfish
    output_role: vcf
```

Generated queries must be produced by the enabled producer in the **selected
attempt**, with its validated receipt. Old files or a successful scheduler exit
cannot substitute for that receipt. Generated assembly-SV queries require an
`assembly_sv.*` producer. To reuse an existing native VCF, replace that row with:

```yaml
  sv:
    sample_id: HG002
    track: sv
    origin: EXTERNAL
    path: /your/project/queries/HG002.sv.vcf.gz
    index_path: /your/project/queries/HG002.sv.vcf.gz.tbi
    source_tool: your_recorded_caller
    source_version: your_recorded_version
    reference_build: GRCh38
```

EXTERNAL means explicitly reused, **not re-called**; its tool/version is user
declared, not re-probed. Query VCF and conventional TBI must be real, readable
BGZF/indexed inputs. Native calls are required: no gVCF substitution or derived
BND review representation. There is no implicit contig renaming or liftover.
The guard streams ALT columns to reject `<NON_REF>` and `<*>` reference-block
alleles even under an ordinary VCF filename; legitimate bare `*` spanning
deletions remain allowed. This narrow record check is not complete VCF or REF
sequence validation, and full-file scanning/hashing adds input-size-dependent I/O.

EXTERNAL 的工具身份来自使用者记录，不能据此追认历史运行。GENERATED 则绑定本次
运行的 producer、路径、完整校验值和回执。参考构建不符必须失败；不能只加 `chr`
前缀就把 GRCh37 数据当作 GRCh38。源 VCF、truth、BED 和 reference 均保持不变。

Example configuration fragment (merge into a valid launcher configuration):

```yaml
benchmark:
  enabled: true
  source_plan: benchmark_sources.yaml
  benchmark_id: HG002_explicit_scope
  sample_id: HG002
  overwrite: false
  sv:
    enabled: true
    truth_vcf: /your/resources/HG002.GRCh38.truth.vcf.gz
    truth_version: your_pinned_truth_version
    truth_source: your_recorded_official_source
    confident_bed: /your/resources/HG002.GRCh38.confident.bed
    confident_bed_version: your_pinned_region_version
    truvari_executable: /your/envs/truvari/bin/truvari
    threads: 1
    memory_mb: 8000
```

Use independently matched, legally available truth resources. Names above are
placeholders, not downloadable resources or recommended matching thresholds.
Keep all configured native thresholds explicit and consistent with the study.
Enabled source-plan categories must match exactly: `small_variants` uses track
`small`; `sv` and `assembly_sv` use track `sv`; `tr` uses track `tr` and the same
explicit `catalog_id: sha256:<catalog digest>` as `benchmark.tr.catalog_id`.
Small variants require a pinned confident BED. Owned SV/assembly comparisons also
require an explicit pinned BED; `latest` is rejected. Configuring both source_plan
and legacy query_vcf is rejected. Legacy query_vcf remains standalone-only.

Resources are frozen before execution. Header sample, contig names/lengths,
explicit reference build, index structure and BED bounds are checked. This is
**not authentication of a provider's truth resource, a full REF-allele sequence
check, or proof that two biological samples are identical**. Independently
validate reference sequence/resource identity and the evaluation scope before
interpreting metrics. No benchmark is appropriate for unmatched SCZ truth.

## Execution and evidence / 执行与证据

```bash
hifivar --config analysis.yaml --dry-run run --run-dir /your/project/runs/bench01
hifivar --config analysis.yaml run --run-dir /your/project/runs/bench01 --resume
```

With source_plan, dry-run reports `BENCHMARK_SOURCE_PLAN_ONLY`: the downstream
query paths and comparison are not claimed to have run. A real invocation
materializes exact query paths only after source receipts are available, then
runs comparisons serially within the benchmark consumer job. Truth/query/BED
hashes, tools, native commands, metric artifacts and source-receipt links are
sealed in the selected-attempt benchmark receipt. TR-only exact comparison is
marked IN_PROCESS, not native third-party execution. Enabled native job resource
reservations must fit the launcher's total budget.

输出位于 `results/benchmark/<benchmark_id>/`，包含分轨结果、原生比较器输出、汇总
manifest 与 metrics TSV；`attempts/<attempt_id>/benchmark/` 保存实际配置、来源和
衔接记录。人工复核、病理解释与 benchmark 仍是独立语义。

Report completion means verified execution/output bytes. Inspect per-comparison
PASS/PARTIAL, the original metrics, filter (`PASS` versus `ALL`), variant class,
region, truth version and denominator. A native comparison is scoped evidence,
not universal accuracy. Overlapping descriptive strata must not be pooled as
independent observations. Reusing a TR VCF as its own truth is forbidden by file
identity, but matching bytes in a separate file is **not** independent truth.

工程小区域验收不能当作全基因组性能；单个 HG002 数据集也不能证明普遍最低深度、
临床有效性或测序成本下降。失败不自动重试；修复配置或软件后使用新的可追溯目录。

## Valid zero-refinement result / 合法的无需细化结果

With Truvari 5.4.0, a successful bench --refine may report zero eligible
refinement regions and exit without writing refined summaries. The new
development contract recognizes only the narrowly verified perfect-match case:
the current command exits zero; native parameters match the request; native
summary counts and rates reconcile with all four readable, indexed assignment
VCFs; both TP sides are positive and FP/FN are zero; the native logs establish the zero-region
exit; and no partial refined or phab outputs exist.

Such a result records refinement_status: NOT_NEEDED and uses the actual
summary.json as **bench** metrics. refined_metrics stays empty. No refined
summary is manufactured or copied from the bench summary. Receipts retain the
proof and its file identities; reuse rechecks them rather than trusting a status
flag. Ordinary completed refinement continues to use its native refined summary.

该分支只处理有完整证据的 Truvari 5.4.0 完全匹配情形，不把一般的缺文件视为成功。
非零退出、部分 refined 输出、日志或参数不符、分区计数不符等仍失败；其他原生
“无需细化”情形暂不放宽。这里的完全匹配仅指该比较口径的计数，不代表基因型完全一致。
比较不同单元时保留每行实际指标来源，不能把 bench
行统一标成 refined，也不能据此增加“算法提高准确率”的结论。

This is an RC5 output-contract fix. Previously failed attempts remain
failed in their historical reports; validate a newly identified candidate in a
separate attempt. A matching package version string alone is insufficient.
