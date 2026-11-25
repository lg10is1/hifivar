# Explicit BND representation bundles / 显式 BND 派生表示

[简体中文](zh_CN/bnd_representation.md)

This development-source increment provides a separate Python API and an
explicitly configured Phase 9 workflow option, disabled by default. It is not a
new caller. Existing native Jasmine VCF output roles and historical release
artifacts are unchanged; source support does not imply a published package.

本开发源码增量提供独立 Python 接口，以及默认关闭的 Phase 9 显式配置选项。
不新增 caller，不自动改变现有流程；原生 VCF 角色和历史 Release 资产不变。

## Phase 9 opt-in / 显式启用

```yaml
sv:
  harmonization:
    enabled: true
    bnd_representation:
      enabled: true
      policy_id: HIFIVAR_JASMINE_BND_REPRESENTATION_V1
      policy_version: 1
      policy_sha256: 9fa679df41dfadc07d7ba06ffea2e0217e3663746b72674f603d1499a61f2f32
      mode: OMIT_BND_END_AND_SVLEN_V1
      timeout_seconds: 300
```

This is part of the existing explicitly enabled within-sample harmonization
workflow, not a new standalone CLI subcommand. Ordinary generated read-SV
selection and its prerequisites remain required by `hifivar run`. External and
assembly sources keep their independent workflow interface.

The per-sample bundle is written beside native outputs at
`results/sv_harmonized/SAMPLE/bndrepr-v1/`. The rule declares that directory only
when enabled. Phase 9 provenance records the exact policy and bundle result;
attempt receipts record all ten mandatory bundle files separately from native
roles. Dry-run records only planned paths, never successful derivation.
An existing bundle is rejected before running Jasmine for any sample, even if
native overwrite was explicitly allowed. No-overwrite is not disabled by a
workflow resume flag. The Snakemake interface also rejects an existing bundle
while constructing the DAG, before scheduler output cleanup can occur. Reusing
an already completed bundle is not implemented by rerunning this producer.

默认关闭时不生成该目录或相关回执文件。启用后，每个样本各自保存完整包；策略身份
必须明确且完全匹配。已有包会在执行任何样本的工具前拒绝复用，native 的 overwrite
设置不能授权覆盖此包。dry-run 只显示计划，不能声称已派生。

Existing per-source Truvari concordance still consumes **native** Jasmine
outputs; this option does not silently change those comparisons or a truth
benchmark. Annotation, review, cohort and benchmark consumers are not
automatically redirected. Explicit derived-role/junction consumers remain a
separate integration gate.

## Contract / 数据契约

The supported policy is `HIFIVAR_JASMINE_BND_REPRESENTATION_V1`, version 1,
mode `OMIT_BND_END_AND_SVLEN_V1`. Its canonical semantic JSON SHA256 is:

`9fa679df41dfadc07d7ba06ffea2e0217e3663746b72674f603d1499a61f2f32`

Its normative grammar is pinned to the [VCFv4.4 specification](https://github.com/samtools/hts-specs/blob/c101c797106d8a1006b9b06f75e8bcc6062cf81e/VCFv4.4.tex).
Omitting BND END/SVLEN is this software's explicit compatibility policy, not a
claim that the specification mandates this edit for every BND producer.

Input is a completed, consistently bound `JasmineResult` from Jasmine 1.1.5,
with its native VCF/TBI, raw VCF, evidence table, source order and reference.
A fresh absolute output directory is required; overwrite is unsupported.
This result object is an API binding, not a cryptographic attestation of an
arbitrary object constructed by an untrusted caller.

输入必须是已完成且身份一致的 `JasmineResult`，而不是预期将来出现的文件路径。
输出目录必须为全新绝对路径；不支持覆盖。结果对象的字段一致性检查不等于对任意
外部构造对象的密码学认证。

Only BND INFO `END` and `SVLEN` tokens are removed in the derived VCF. Two
explicit provenance INFO fields identify the policy and mandatory ledger row.
Record count/order, CHROM/POS/ID/REF/ALT/QUAL/FILTER/FORMAT/sample fields,
all other INFO tokens, SUPP_VEC and IDLIST are preserved. Non-BND records,
including existing POS0 insertion representations, retain their exact bytes.
No filtering, coordinate clamping, allele normalization or synthetic mate link
is performed. Unsupported or malformed grammar fails the whole derivation.

派生文件只移除 BND 的 `END`、`SVLEN` 原始 INFO 项，并追加两个溯源字段。
不删记录、不重排、不改变断点、等位基因、基因型或 caller 成员关系。非 BND 记录逐字节
不变。原有 POS0 插入不会被平移到 1。非法或未支持表示会明确失败，而非静默跳过。

## Bundle and validation / 产物与验证

- Binary copies of native VCF/TBI, raw VCF and native evidence.
- A separate `SAMPLE.bndrepr-v1.sv.vcf.gz` and TBI (bgzip/tabix >=1.23.1).
- `records.jsonl`: every record's ordinal, hashes, action and reversible INFO ledger.
- `junctions.tsv`: BND local/mate anchors, orientations and separate mate-link QC.
- `native_result.json` and a success `manifest.json` published only after validation.
- `FAILED.json`, without a success manifest, on an ordinary execution failure.

The native-plus-ledger package is reversible: derived records plus the ledger
must reproduce every decompressed native data-record byte, including LF/CRLF.
The standalone derived VCF is **not lossless**. Native compressed bytes are
preserved by binary copy, not promised to be reproduced by re-compression.
Sequential BGZF readback and bounded local index queries are different checks.
Virtual boundaries and contigs with nonordinary region-query syntax are recorded
as not tested by ordinary interval queries, never silently converted.

“可逆”指整个原生文件＋派生文件＋账本的数据包；单独派生 VCF 并非无损副本。
BGZF/TBI 成功不证明事件正确。局部染色体索引也不会按 ALT 中另一染色体的坐标检索。

## Interpretation and consumer boundary / 解释与下游边界

MATEID resolution, reciprocity, orientation, exact coordinate agreement and
explicit CIPOS-aware consistency are recorded separately. Unresolved original
IDs are not silently renamed. Coordinate inequality alone does not invalidate
a mate relation in the presence of microhomology/uncertainty.

The junction table contains two VCF anchors, not an event-spanning interval.
Existing interval annotation/review consumers must not read native BND END as
a local event span or assume a virtual 0/N+1 anchor is an ordinary interval.
This API does not automatically redirect those consumers. Their explicit
junction-aware integration remains a separate gate. INS/DEL benchmarks remain
INS/DEL-only; this representation cannot establish BND accuracy.

双断点表不等于跨区间注释。现有注释、人工复核和 benchmark 不会自动改用派生文件；
下游必须显式选择产物并通过其自身的双端坐标契约。caller 支持数、断点表示兼容性、
人工证据和真实准确率不能混为一谈。

## Explicit selected consumers / 显式选择的消费者

Development annotation/review consumers now distinguish native VCF and an explicitly verified BND bundle, query anchors independently, and retain unknown/virtual endpoints. No automatic consumer redirection or BND accuracy claim is added. See [junction-aware consumers](bnd_junction_consumers.md). Earlier separate-gate descriptions above preserve the original representation-only boundary.
