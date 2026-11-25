# 基准测试

[English](benchmarking.md) | [简体中文](benchmarking.zh-CN.md)

本文档描述 HiFiVar 发布候选版本 0.1.0rc4 的已验证 HG002 基准测试结果。所有指标均源自已冻结、已审计的工件，且特定于以下文档中记录的确切工具版本、参数、参考基因组构建和真相数据集。

## 重要范围声明

- 基准测试性能是**数据集特定、版本特定、配置特定和区域特定**的。不同的配置、样本或技术可能产生不同的结果。
- CMRG 读段 SV 基准测试仅覆盖**具有挑战性的医学相关基因**，**不是**全基因组 SV 真相基准测试。
- 没有匹配真值集的样本结果，不作为准确率基准。
- 调用器支持计数是描述性证据，不是生物学真相。
- 注释或功能影响不是调用置信度。
- 不提供临床解释。

## 1. 全基因组小变异基准测试

### 配置

| 项目 | 值 |
|------|-----|
| 样本 | GIAB HG002 |
| 数据 | PacBio HiFi CCS 11kb |
| 参考基因组 | GRCh38 (human_GRCh38_no_alt_analysis_set) |
| 评估范围 | chr1–chr22 置信区域 |
| 调用器 | DeepVariant 1.10.0 (PACBIO 模型) |
| 基准测试工具 | hap.py 0.3.15 |
| hap.py 运行时 | Python 2.7.15 传统验证环境 |
| DeepVariant 分片数 | 4（从 8 减少，因资源限制；无科学影响） |

### 规范指标

| 变异类别 | 精确度 | 召回率 | F1 |
|----------|--------|--------|-----|
| SNP | 0.995444 | 0.977927 | 0.986608 |
| INDEL | 0.916605 | 0.893150 | 0.904726 |

组合指标（精确度=0.9843，召回率=0.9665，F1=0.9753）是**非规范的**，不得用作标题性能。它是从 SNP 和 INDEL 行派生的，仅用于内部比较。

### 分层小变异基准测试

结果按 GIAB/GA4GH v3.1 基因组难度类别分层。仅**所有困难区域**和**不在所有困难区域**可作为互补类别呈现。不要重叠分层区域。

#### SNP 按分层

| 分层 | 精确度 | 召回率 | F1 |
|------|--------|--------|-----|
| all_difficult_regions | 0.988451 | 0.971033 | 0.979665 |
| not_in_all_difficult_regions | 0.997079 | 0.979560 | 0.988242 |
| low_complexity | 0.983696 | 0.968712 | 0.976146 |
| homopolymer | 0.977721 | 0.963155 | 0.970384 |
| tandem_repeats | 0.991771 | 0.974935 | 0.983281 |
| segmental_duplications | 0.973572 | 0.956138 | 0.964776 |
| mappability | 0.929612 | 0.846068 | 0.885875 |
| GC_content | 0.995844 | 0.979658 | 0.987685 |
| functional_regions | 0.994127 | 0.973027 | 0.983464 |

#### INDEL 按分层

| 分层 | 精确度 | 召回率 | F1 |
|------|--------|--------|-----|
| all_difficult_regions | 0.888783 | 0.860543 | 0.874435 |
| not_in_all_difficult_regions | 0.984558 | 0.970714 | 0.977587 |
| low_complexity | 0.881102 | 0.851154 | 0.865869 |
| homopolymer | 0.853250 | 0.823097 | 0.837902 |
| tandem_repeats | 0.958893 | 0.927645 | 0.943010 |
| segmental_duplications | 0.950951 | 0.938938 | 0.944906 |
| mappability | 0.866345 | 0.804769 | 0.834422 |
| GC_content | 0.923517 | 0.901508 | 0.912380 |
| functional_regions | 0.813596 | 0.937028 | 0.870961 |

除非有单独分析支持，否则不要断言基因组属性导致了较低的准确性。仅描述观察到的差异。

## 2. 基于读段的 SV 基准测试 (CMRG GRCh38)

### 配置

| 项目 | 值 |
|------|-----|
| 样本 | HG002 |
| 真相 | GIAB CMRG GRCh38 v1.00 StructuralVariant |
| 范围 | 仅具有挑战性的医学相关基因 |
| 基准测试工具 | Truvari v5.4.0 |
| 最小 SV 大小 | 50 bp |
| 原始 CMRG 记录数 | 250 |
| 有资格的基准测试基础计数 | 216 |
| 排除 | 34 条低于 50 bp 的记录 |

### 为什么有 250 条记录但基础计数是 216？

CMRG 真相 VCF 包含 250 条总 SV 记录。Truvari 应用 `--sizemin 50`，按 SV 大小过滤记录。只有 SV 大小 >= 50bp 的记录有资格进行基准测试，将基准测试范围减少到 216 条记录。排除的 34 条记录是小于 50bp 的小变异，不属于标准 SV 基准测试范围。

### 前 NIST SV v0.6 基准测试

使用 NIST SV Tier1 v0.6 的先前基准测试是 **INVALID_REFERENCE_BUILD_MISMATCH**：

- 真相：GRCh37
- 查询调用/参考：GRCh38
- 染色体重命名未执行坐标转换

所有先前针对 NIST v0.6 的每调用器 Truvari 指标对于 GRCh38 性能解释无效，不得出现在性能表或图中。

### 调用器指标

| 调用器 | 版本 | TP | FP | FN | 精确度 | 召回率 | F1 | 基因型一致性 |
|--------|------|-----|-----|-----|--------|--------|------|--------------|
| Sawfish | 2.2.1 | 196 | 0 | 20 | 1.000 | 0.907 | 0.951 | 0.944 |
| Sniffles2 | 2.8.0 | 160 | 7 | 56 | 0.958 | 0.741 | 0.836 | 0.875 |
| pbsv | 2.11.0 | 142 | 32 | 74 | 0.816 | 0.657 | 0.728 | 0.873 |
| cuteSV | 2.1.4 | 34 | 3 | 182 | 0.919 | 0.157 | 0.269 | 0.882 |

**历史 Sawfish 样本身份限制。** 其原始及导出 VCF 的物理样本名为
`UnknownSampleName`，从原始 reads 到派生 BAM 的完整身份链尚未闭合。表中数值作为
历史结果保留，不计作严格通过样本身份验证的比较点。制作比较图或调用器排名时，应
排除该点或明确标记限制；后来的区域 Sawfish 运行不能补证这个旧产物的身份。

每个公开声明都必须说明这些是 **CMRG 区域结果**，**不是全基因组 SV 性能**。不要声称一个调用器普遍优于另一个。

### Truvari 参数

```
--sizemin 50
--sizemax 100000
--refdist 500
--pctseq 0.5
--pctsize 0.5
--includebed HG002_GRCh38_CMRG_SV_v1.00.bed
```

### 调用器输入 BAM

| 调用器 | BAM | 状态 |
|--------|-----|------|
| Sawfish | HG002.nohardclip.sorted.bam | 已移除硬剪切（工具兼容性工件） |
| Sniffles2 | HG002_ccs11kb.sorted.bam | 原始（规范原始比对） |
| pbsv | HG002.nohardclip.rg.sorted.bam | 已移除硬剪切 + RG 标签（工具兼容性工件） |
| cuteSV | HG002_ccs11kb.sorted.bam | 原始（规范原始比对） |

## 3. 可重复性

所有基准测试工件都有校验和。冻结的证据文件是：

- `HG002_BENCHMARK_FROZEN_SUMMARY.tsv` — 机器可读摘要
- `HG002_BENCHMARK_FROZEN_METADATA.yaml` — 包含校验和的完整元数据
- `HG002_BENCHMARK_FROZEN_CHECKSUMS.sha256` — 所有工件校验和
- `HG002_BENCHMARK_INTERPRETATION.md` — 解释指南

详见[可重复性](reproducibility.md)文档。

## 4. 局限性

详见[已知局限性](known_limitations.md)，包括：

- 单一 HG002 数据集
- 一种 PacBio CCS 11kb 配置
- CMRG 受限区域 SV 范围
- 缺少全基因组 GRCh38 SV 真相基准测试
- 传统 hap.py 环境
- 不提供临床解释

## 5. 引用

引用基准测试结果时，请引用各个工具：

- **DeepVariant**: Poplin et al., Nature Biotechnology 36, 983–987 (2018). doi:10.1038/nbt.4235
- **hap.py**: Li et al., Nature Methods 14, 587–590 (2017). doi:10.1038/nmeth.4326
- **Truvari**: Alonge et al., Genome Biology 24, 168 (2023). doi:10.1186/s13059-023-02999-w
- **Sawfish**: Huang et al., Genome Research 33, 1–17 (2023). doi:10.1101/gr.277943.123
- **Sniffles2**: Sedlazeck et al., bioRxiv (2023). doi:10.1101/2023.12.27.573420
- **pbsv**: PacBio, GitHub https://github.com/PacificBiosciences/devbio/wiki/pbsv
- **cuteSV**: Jiang et al., Nature Communications 11, 6251 (2020). doi:10.1038/s41467-020-19985-w

## 引用验证

CITATION_VERIFICATION_REQUIRED: 正式发表前验证所有 DOI 和 URL。
# Phase 15 低深度 HiFi 实验

可选的 Phase 15 使用“参数 profile × 目标深度 × 重复次数”矩阵，并由用户
预先登记 SNP/INDEL 等类别的准确度阈值。它支持已比对 BAM 工程模式和原始
FASTQ 分子级嵌套抽样模式；后者用于独立重新比对、calling 和 truth benchmark。
它区分目标深度与实际观测深度，不会把一次实验或单个样本结果宣称为通用最低
测序深度，也不会自动声明某个参数 profile 最优。具体契约见
[Phase 15 文档](phases/phase15.md)。
