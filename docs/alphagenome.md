# AlphaGenome: explicit user-key functional prediction / 用户自有凭据功能预测

[简体中文](zh_CN/alphagenome.md)

This page describes **RC5 release source**, not a capability retroactively added to RC1–RC4 or an older same-version RC5 wheel. Match the source manifest and candidate artifact SHA256 before installation.

## Scope / 范围

The new command uses the official AlphaGenome **0.10.0** v1 SDK to score an explicit candidate TSV. It preserves scorer, modality, tissue/ontology, gene/strand and raw/quantile score rows rather than collapsing them into a calling-confidence score. It is not a new caller or a clinical pathogenicity classifier.

当前支持人类GRCh38、chr1–chr22/chrX/chrY，以及明确的ACGT等位序列。small候选可以包括SNV/Indel；标为SV的候选仅接受有共同锚点、长度差≥50bp且序列已解析的INS/DEL。单个allele长度上限为10,000bp。符号ALT、BND、INV、CPX及未经等位序列重建的TR不支持；不自动重命名contig、liftover、重建插入序列或更改原始调用。

The local FASTA/FAI must match every selected REF and contain the entire model context. This does not prove that every base of the local FASTA is identical to the provider's reference asset. The model checkpoint is explicitly named; no server-weight digest is inferred.

## Install / 安装

From the matching RC5 source checkout:
```bash
python -m pip install '.[alphagenome]'
hifivar functional alphagenome --help
```

The optional extra pins alphagenome==0.10.0. Core HiFiVar remains usable without it. SDK installation/import/offline contract tests and an authenticated service prediction are separate gates; developer tests use fake transport and do not establish biological accuracy.

Use the [plan](../examples/alphagenome/plan.template.yaml) and [selection header](../examples/alphagenome/selected_variants.template.tsv) as templates; they require actual user inputs before execution.

## Plan / 配置

```yaml
schema_version: 1
selection_file: /your/project/selected_variants.tsv
reference:
  fasta: /your/project/reference/GRCh38.fa
  build: GRCh38
alphagenome:
  model_version: ALL_FOLDS
  requested_modalities: [RNA_SEQ, ATAC]
  ontology_terms: [UBERON:0000955]
  sequence_length: 1048576
  max_variants: 20
  max_rows_per_variant: 100000
  max_allele_length: 10000
  credential_env: HIFIVAR_ALPHAGENOME_API_KEY
  connect_timeout_seconds: 15
  rpc_timeout_seconds: 60
data_sharing_consent: false
```

Supported checkpoint names: ALL_FOLDS, FOLD_0, FOLD_1, FOLD_2, FOLD_3. Sequence lengths in this pinned SDK are 16,384, 131,072, 524,288 or 1,048,576 bases. Use documented uppercase modality names from the SDK; the adapter exposes ATAC, CONTACT_MAPS, DNASE, CHIP_TF, CHIP_HISTONE, CAGE, PROCAP, RNA_SEQ, SPLICE_SITES, SPLICE_SITE_USAGE and SPLICE_JUNCTIONS.

The adapter retains separate stranded-gene tracks and does not restrict tidy
rows to matching gene/track strands. These recorded combinations are not an
additional biological strand-matching filter; assess the retained metadata for
the intended downstream interpretation.

当前适配器保留独立链轨道及基因/轨道链组合，未额外进行链匹配生物学筛选。

**Ontology selection is a post-scoring result filter.** The SDK score_variant call has no ontology argument; choosing a tissue does not reduce remote inference or limit which model tracks the service computes. An empty matching result is reported explicitly, not converted into a zero effect.

候选TSV必须恰好包含以下列：
```text
sample source_variant_id variant_category contig position reference_bases alternate_bases source_annotation selection_reason
```

Column names must be unique, with exactly these nine fields per data record.
Reordering the columns and valid TSV quoting are supported. Duplicate headers,
missing or extra data fields and malformed quoted records are rejected before
credential access or scoring. Errors do not echo record contents; field-count
and value errors identify the data-record index.

列名必须唯一，每条数据记录必须恰好有九个字段；可重排列顺序并使用合法TSV引号。
重复列、缺失或多余字段以及不完整引号记录会在读取凭据或评分前失败。

Fields above are tab-separated, not space-separated. Position is 1-based VCF POS; category is small or sv. source_annotation is an existing provenance file, resolved relative to the TSV when not absolute. Use actual selected records from your accepted annotation/callset, not illustrative alleles. Each sample/source_variant_id pair must be unique. The input TSV and caller outputs are not rewritten.

## Offline first / 先离线核对

```bash
hifivar --dry-run functional alphagenome --plan alphagenome.yaml --output-dir /your/project/functional/run01
```

Dry-run validates the plan, explicit variant forms, existing annotation provenance, REF and full context against FASTA/FAI. It does not inspect an API key, import the optional SDK, create the output directory or send requests.

## Supply your own key / 使用自己的API授权

Acquire service authorization through the official AlphaGenome service. This adapter currently supports its official API-key authentication; OAuth/ADC and arbitrary third-party endpoints are not claimed.

For an interactive Linux Bash session, avoid placing the literal key in a command argument or a repository file:
```bash
read -r -s -p "AlphaGenome API key: " HIFIVAR_ALPHAGENOME_API_KEY
export HIFIVAR_ALPHAGENOME_API_KEY
```

Set data_sharing_consent: true only after deciding that the selected coordinates/alleles may be sent to the official service. Then explicitly enable network execution:
```bash
hifivar functional alphagenome --plan alphagenome.yaml --output-dir /your/project/functional/run01 --allow-network
unset HIFIVAR_ALPHAGENOME_API_KEY
```

The adapter reads only the explicitly named variable after opt-in; it does not fall back to legacy AlphaGenome variables. Do not put a literal API key in YAML/TSV, CLI arguments, reports or Git. User service access/quota remains separate from software installation.

## Outputs, timing and failures / 输出和失败

A new output directory is required. Dry-run returns PLANNED without creating files. For real execution, alphagenome.result.json records RUNNING/COMPLETE/FAILED as appropriate, explicit selection/settings, local REF/context checks, SDK/model identifiers and structured scores. A failed request retains a sanitized failure report and makes no successful-prediction claim. Provider exception text and credential values are not persisted.

The connect timeout is separate from the per-RPC deadline. The SDK can make up to five transient attempts; HiFiVar adds no retry and does not claim a hard deadline for the entire variant batch. At most20 selected variants are accepted. Actual quota/latency depends on the user's service account and context size.

This is a separate explicit post-annotation command. annotation.functional_enabled is not an automatic hifivar run branch; the launcher points users to this route. Functional impact, call confidence, truth matching and manual review remain independent.

## Primary references

- [Official pinned SDK v0.10.0](https://github.com/google-deepmind/alphagenome/tree/v0.10.0)
- [Official DNA client and scorers](https://www.alphagenomedocs.com/api/models.html)
- [Official variant-scoring explanation](https://www.alphagenomedocs.com/variant_scoring.html)
