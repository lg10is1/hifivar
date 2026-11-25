# Joint small/SV/TR phasing / 小变异、SV与TR联合定相

[简体中文](zh_CN/joint_phasing.md)

This guide describes a **RC5 release source**, not RC1–RC4 or
every wheel reporting `0.1.0rc5`. Small-only HiPhase has separate native evidence.
Identified installed-source regional joint cases have independently passed,
including a later two-technical-label batch with three cohort tracks. They are
not evidence for every same-version/later package or whole-genome scope;
local fake-tool and DAG tests alone are not native acceptance.

本说明对应 GitHub RC5 发行包。已有小变异定相运行不等于多类别定相已经验收，
相同版本号也不证明包内容相同。请以所用包的来源清单和实际验收范围为准。

后续指定安装包已有独立的区域联合定相与两技术标签批量验收，二者各保留自己的
attempt和包身份；并非所有新包都已原生验收，也不等于定相准确率。完整配置及
范围见[批量多类型指南](batch_multivariant_run.md)。

## Inputs and configuration / 输入与配置

HiPhase receives the same indexed single-sample BAM/reference, a DeepVariant VCF,
and optionally **one** pbsv or Sawfish VCF and/or a TRGT VCF. These are distinct
inputs to a single phasing invocation, not a merged VCF. Sniffles2, cuteSV, PAV,
SVIM-asm and harmonized outputs are not silently treated as supported HiPhase inputs.

Add the following to an otherwise complete, validated analysis configuration:

```yaml
small:
  enabled: true
sv:
  enabled: true
  pbsv:
    enabled: true
tr:
  enabled: true
  catalog: /your/project/resources/repeats.bed
  catalog_reference_build: GRCh38
  karyotype: XX  # Example declaration only; set from actual sample metadata.
phasing:
  enabled: true
  sv_caller: pbsv  # null (default), pbsv, or sawfish; explicit choice, not auto-selection.
  include_tr: true  # false by default.
  global_realignment: auto  # auto, enabled, or disabled.
```

Enable only intended SV callers in the complete config. Selecting `sv_caller`
requires that caller's producer to be enabled; `include_tr` requires TRGT.
The launcher and direct DAG reject disabled/missing dependencies rather than
reenabling producers. Small calling remains required for this workflow route.
The independent `run_phase6` API can instead accept explicit validated artifacts
via `sv_variant_artifacts` and `tandem_repeat_artifacts`, with the same selections.

中文：上例只是需要合并进完整配置的片段。参考、样本表、工具路径、目录和资源需另行准备；
`XX`不是性别推断，必须按真实样本声明设置。选择的SV caller或TRGT被关闭时直接报错，
不会为了定相自动开启其他caller。可以只加SV、只加TR，或二者都加；默认仍只接小变异。

## Separate outputs / 分类型输出

| Role / 角色 | Output / 输出 |
| --- | --- |
| Small variants / 小变异 | `phasing/{sample}.phased.vcf.gz` and `.tbi` |
| Selected SV caller / 所选SV来源 | `phasing/{sample}.sv.phased.vcf.gz` and `.tbi` |
| TRGT / 串联重复 | `phasing/{sample}.tr.phased.vcf.gz` and `.tbi` |

Only selected extra outputs are declared. The wrapper sends paired repeated
`--vcf` / `--output-vcf` flags in small, SV, TR order, runs tabix independently
for each output, and records all outputs in the same phasing receipt. Original
caller VCFs/indexes are not overwritten; the bridge rechecks their identities.

中文：三个类别始终保留独立文件。同一次HiPhase可以使用共享的相位信息，但输出VCF
不是新的calling结果，也不能把不同类别或phased/unphased文件重复相加为变异总数。
报告和后续注释不能在未声明的情况下以phased文件替换原始caller输出。

## Realignment and validation / 重比对与校验

- `global_realignment: auto` preserves the previous small-only
  `--disable-global-realignment` behavior. Adding SV/TR enables the upstream
  global-realignment default. Explicit enabled/disabled choices are recorded.
- Joint inputs and outputs require the exact sample and reference contig names
  and matching declared lengths. SV/TR category headers, BGZF and TBI are checked;
  these checks do not prove biological identity or complete REF sequence equivalence.
- Joint output validation streams each ordered CHROM/POS/ID/REF/ALT key and
  unordered GT alleles against its own source. It permits phase orientation,
  PS and PF changes; it does not promise byte-identical non-phase annotations.
- Existing output/index or source aliases fail before commands unless overwrite
  is explicitly permitted; source aliases remain forbidden even with overwrite.
- Failure or a missing category output does not create a successful phasing receipt.

HiPhase may leave a variant unphased or mark `PF=TR_OVERLAP` when a TR overlaps a
small/SV call. That is upstream phase-selection behavior, not evidence that the
variant disappeared or was independently rejected. GT phase, phase sets, caller
support, functional impact and manual review are different evidence types.

中文：`auto`在小变异单轨时维持旧参数；加SV/TR后使用上游推荐的全局重比对默认行为，
可能增加资源需求。所有资源仍须按数据设置，不能把小型测试配置当作通用生产配置。
记录保真只核对上述身份字段和GT等位集合，不等于相位正确率，也不代替switch/flip
benchmark、独立相位真值或人工复核。没有相位真值时必须保持准确率未评估。

## Evidence and scope / 验收边界

This increment does not change calling algorithms, truth sets, matching policies,
historical metrics or releases. It adds neither haplotagged BAM output nor pedigree
phasing. Exact-package native joint acceptance must identify its
inputs, configured categories, native command, indexes and actual output content;
it cannot be borrowed from earlier small-only execution.

Official interface: [HiPhase user guide](https://github.com/PacificBiosciences/HiPhase/blob/main/docs/user_guide.md).
See also [capabilities](capabilities.md) and [historical Phase 6](phases/phase06.md).
