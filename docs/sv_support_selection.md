# Caller-support selection / 按 caller 支持数筛选

[简体中文](zh_CN/sv_support_selection.md)

This option belongs to RC5. Historical releases and
earlier packages with the same version string do not acquire it automatically.
Use the exact accepted source manifest and package checksum.

本选项属于当前开发源码。历史 Release 与先前同版本候选包保持不变，不能仅凭
版本号判断是否包含该功能。原生完整 Jasmine 输出始终保留。

## Configure a run / 运行配置

Add this setting to a complete [within-sample harmonization configuration](sv_harmonization_run.md):

```yaml
sv:
  harmonization:
    enabled: true
    minimum_caller_support: 2
```

```bash
hifivar --config analysis.yaml --dry-run run --run-dir /your/project/runs/support2
hifivar --config analysis.yaml run --run-dir /your/project/runs/support2 --resume
```

The first command freezes and previews a new run; execute its matching identity
with the second command. Use an unused run directory for a different threshold.
Once any support bundle has been created, the workflow refuses to reuse that
output root. An interrupted batch with existing support bundles therefore needs
a new attempt/output root; unrestricted in-place resume is not supported.
The setting is a strict integer from 1 to 6, additionally bounded by the actual
selected input callers for **each** sample. Booleans, fractions and numeric
strings are rejected. A batch may select different source sets per sample, but
the configured threshold must be feasible for all of them.

门槛为“至少 k 个不同 caller 来源支持”，不是恰好 k 个，也不是 read 数、样本数
或单倍型数。仅参与本次 Jasmine 的 COMPLETED/NO_CALLS 来源占据 SUPP_VEC 的
位置；禁用、未运行或失败来源不算可用输入。若某样本只有两个输入，不能设置为3。
改变门槛需新运行目录，不覆盖旧分析。

## Retained outputs / 保留的输出

The default `1` preserves existing behavior and creates no extra selection bundle.
For `k > 1`, HiFiVar first completes the native Jasmine VCF/evidence, then creates:

```text
results/sv_harmonized/<sample>/
  <sample>.harmonized.sv.vcf.gz[.tbi]     # complete native result
  <sample>.sv.evidence.tsv               # complete native evidence
  support-ge-<k>/
    <sample>.support-ge-<k>.sv.vcf.gz
    <sample>.support-ge-<k>.sv.vcf.gz.tbi
    <sample>.support-ge-<k>.sv.evidence.tsv
    manifest.json
```

Membership uses the popcount of `SUPP_VEC` in the exact actual Jasmine input-list
order. If `SUPP` is present, it must agree. The complete evidence table must
agree row by row before a subset can be published. Retained VCF header/records
and evidence rows are byte-preserved; no alleles, genotypes or filters are
rewritten. Original caller and native harmonized files are hashed before/after.
The manifest records the threshold, source order, counts, versions, commands,
input/output checksums and validation. It is written last, only on success.

筛选结果为独立派生包，目录必须全新，即使其他设置允许 overwrite 也不能覆盖。
失败文件保留且不签发成功 manifest。空集合允许，但仍需合法 header、BGZF 和索引。
完整原生 VCF、evidence 及来源记录保留，便于重新选择门槛和追溯。

## Consumer boundary / 下游边界

Existing native output roles, BND representation, per-source Truvari concordance
and downstream consumer routing continue to use the complete native result.
They are not silently redirected to a subset. A downstream analysis must
explicitly select and identify the derived artifact. The current selection is
recordwise: BND mate completeness is **NOT_ASSESSED**. Do not substitute this VCF
for a paired-breakend/junction bundle without a separate pairing validation.
Benchmarks must be rerun on each selected query using matched truth, regions and
fixed comparison settings; do not reuse TP/FP labels from the unfiltered query.

现有自动下游仍绑定完整原生结果；使用筛选集合时需明确指定该派生文件并单独验收。
本功能不会自动运行金标准 benchmark，也不改变 Jasmine 聚类距离或配对规则。

## Choosing a threshold / 如何选择

- Discovery: retain all candidates and inspect the evidence; singleton support
  does not prove a call is false.
- Balanced selection: compare precision, recall and F1 on matched benchmark
  data, then lock a threshold for independent evaluation.
- Precision-prioritized selection: test a higher threshold and report its recall
  cost. More support is not a calibrated probability of correctness.

研究目的不同可以选择不同门槛。开发者的探索性阈值曲线不应变成所有数据的默认
推荐：深度、caller 组合、组装质量和变异类型都会改变结果。多个 caller 共用相同
测序数据与比对，不能假定其证据独立。支持数不是置信度、致病性或临床验证结果。
应在独立样本、区域或预先划分的数据上确认所选策略，再作普遍结论。
