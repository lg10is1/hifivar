# Explicit SV source identity / 显式 SV 来源身份

[简体中文](zh_CN/sv_source_identity.md)

RC5 includes this opt-in interface. Consult the validation record for native execution on a specific build.
Historical release assets and canonical caller VCFs are not rewritten.

本功能由 RC5 提供，默认关闭。单元测试不代表新的已安装原生工作流
已经验收。历史 Release 和原始 caller VCF 不会被回写。

## Why it is separate / 为什么单独处理

VCF IDs need not uniquely identify records. Two records with the same ID can
have different alleles, genotypes, filters or coordinates. Dropping either is
not an identity fix. Conversely, a caller name plus a non-unique ID is not a
complete provenance key. This feature creates an explicit analysis-only view
and records the declared caller namespace, exact source SHA256 and one-based
record ordinal. Distinct callers may legitimately produce identical file bytes;
their provenance keys remain separate without rejecting such concordance.

重复 ID 不等于重复变异。本功能不去重、不丢弃记录，也不修改坐标、等位基因、
FILTER、INFO、FORMAT、GT、phase 或科学阈值。稳定来源键由原始压缩文件的
SHA256、声明的 caller 命名空间与记录序号组成；不同 caller 即使输出完全相同的
文件字节，也保留独立来源键，不因此拒绝合法一致结果。分析视图和原始调用文件
的角色明确分开。

## Configuration / 配置

The strict configuration key is `sv.harmonization.input_identity`. Enabled
runs require the exact packaged policy ID, version, semantic SHA256 and mode.
Use the installed package to print the complete supported mapping:

启用时必须显式声明安装包内策略的 ID、版本、语义 SHA256 和模式。以下命令仅
打印受支持配置，不执行分析；将输出加入自己的 YAML，不能省略或猜测策略哈希。

```bash
python -c 'from hifivar.sv_input_identity import SourceIdentityPolicy,POLICY,POLICY_SHA256; import yaml; p=SourceIdentityPolicy(POLICY["policy_id"],POLICY["policy_version"],POLICY_SHA256,POLICY["mode"]); print(yaml.safe_dump({"sv":{"harmonization":{"enabled":True,"input_identity":p.to_dict()}}},sort_keys=False))'
```

Default:

```yaml
sv:
  harmonization:
    input_identity:
      enabled: false
```

This option belongs to the existing Phase 9 workflow or the supported launcher
harmonization path; it does not enable a caller or a cross-sample SV cluster.
Select sources and executable paths using the ordinary harmonization settings.
The policy changes no Jasmine distance or Truvari benchmark threshold.

该开关进入既有 Phase 9 工作流和受支持的整合入口，不会自动开启 caller，也不
支持跨样本 SV 聚类。来源、程序路径及资源仍由原有配置指定。操作必须使用新的
结果目录；已有 `source-identity-v1` 目录在调度前被拒绝，overwrite 也不能绕过。

## Identity policy / 身份策略

- Every occurrence of a globally duplicated ID, and each missing `.` ID, gets a
  deterministic `HFVID1_<caller>_<source_sha256>_<ordinal>` analysis ID.
- Unique IDs remain unchanged. Existing `HFVID1_` IDs and ambiguous/multi-token
  IDs are rejected, not silently interpreted. Repeated IDs across callers are
  handled using the same explicit policy.
- MATEID/PARID references to duplicated original IDs are rejected. No missing
  mate, parent or biological event relationship is guessed.
- Source and reference identities, sample and contigs, exact non-ID record
  bytes, order and original line endings are verified. Derived views are BGZF
  compressed and indexed, with native index listing and a bounded ordinary
  coordinate query when applicable. Indexability is not biological validation.
- A fresh owned bundle is retained on failure with `FAILED.json`. Only complete
  preparation publishes `manifest.json`. A completed preparation manifest is
  not proof that subsequent Jasmine or a whole sample analysis succeeded.

原始文件始终保持不变。策略只改分析视图的 ID 字段，额外添加一条策略头；每条
记录都保留原 ID、分析 ID、两侧记录哈希及原始记录键。非唯一 mate/parent 引用
明确失败，不通过删 BND 或猜配对来放行。源文件变化、非法坐标或参考不符也失败，
不做 chr 命名转换或坐标钳制。

## Retained outputs / 保留结果

The per-sample `source-identity-v1/` directory contains:

- Ordered `NN.<caller>.analysis.vcf.gz` plus `.tbi` views.
- One `NN.<caller>.analysis.records.jsonl` ledger per source.
- A streaming SQLite identity lookup and a policy/tool/input/output manifest.
- Preparation logs and owned work files. The directory is part of the Phase 9
  output contract; successful launcher receipts bind its retained files.

For the opt-in path, the evidence TSV retains the existing columns and adds:

| Field | Meaning / 含义 |
|---|---|
| `source_record_ids` | Original caller IDs, in actual member order; duplicates and `.` are allowed and are not unique keys / 原始 ID，不承诺唯一。 |
| `analysis_record_ids` | Actual unique IDs in Jasmine IDLIST / 实际工具成员 ID。 |
| `source_record_keys` | JSON list of caller + exact source SHA256 + record ordinal keys / 唯一来源键。 |
| `source_record_callers` | Caller for each member, aligned to the keys / 每个成员对应的 caller。 |
| `source_files` | Canonical files for supporting callers / 原始 caller 文件。 |
| `representation_files` | Actual ordered analysis views for supporting callers / 实际分析视图。 |
| `source_identity_manifest` | Mandatory identity manifest / 必须保留的身份清单。 |

The IDLIST members must match SUPP_VEC by caller name. Every source record must
be assigned exactly once; missing, repeated or unknown members fail evidence
publication rather than becoming plausible-looking provenance. Default-off TSV
columns and interpretation remain unchanged.

SUPP_VEC 仍按真实输入列表解释。证据表按名称和唯一来源键连接，不依赖列位置。
视图单独不是无损交付，必须同时保存原始调用文件、账本和清单。该策略不证明
调用准确率、GT/定相准确率、临床可信度或支持数置信度。BND 表示与后续消费者
仍遵循自己的独立 opt-in 合同。
