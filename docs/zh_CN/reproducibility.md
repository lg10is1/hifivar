# 可重复性

本文档描述如何重现 HG002 基准测试结果并验证所有基准测试工件的完整性。

## 冻结的证据文件

| 文件 | 描述 |
|------|------|
| `HG002_BENCHMARK_FROZEN_SUMMARY.tsv` | 机器可读摘要表 |
| `HG002_BENCHMARK_FROZEN_METADATA.yaml` | 包含校验和的完整元数据 |
| `HG002_BENCHMARK_FROZEN_CHECKSUMS.sha256` | 所有工件校验和 |
| `HG002_BENCHMARK_INTERPRETATION.md` | 解释指南 |

## 验证步骤

### 1. 验证校验和

```bash
sha256sum -c HG002_BENCHMARK_FROZEN_CHECKSUMS.sha256
```

### 2. 验证小变异 VCF

```bash
bcftools view -h HG002_ccs11kb.wgs.small.vcf.gz | head -20
bcftools query -l HG002_ccs11kb.wgs.small.vcf.gz
```

预期：样本 HG002，GRCh38 染色体（chr1-chr22）。

### 3. 验证 CMRG 真相 VCF

```bash
bcftools view -h HG002_GRCh38_CMRG_SV_v1.00.vcf.gz | head -20
bcftools view -H HG002_GRCh38_CMRG_SV_v1.00.vcf.gz | wc -l
```

预期：250 条记录，GRCh38 染色体，样本 HG002。

### 4. 验证调用器 VCF

```bash
for caller in sawfish sniffles2 pbsv cutesv; do
  echo "--- ${caller} ---"
  bcftools view -H HG002.${caller}.sv.vcf.gz | wc -l
  bcftools query -l HG002.${caller}.sv.vcf.gz
done
```

### 5. 验证 Truvari 结果

```bash
for caller in sawfish sniffles2 pbsv cutesv; do
  echo "--- ${caller} ---"
  python3 -c "
import json
d = json.load(open('truvari_cmrG/${caller}/summary.json'))
tp = d['TP-base']
fn = d['FN']
bc = d['base cnt']
print(f'TP={tp} FN={fn} TP+FN={tp+fn} base_cnt={bc} match={tp+fn==bc}')
"
done
```

所有调用器应显示 TP+FN = 216 = base_cnt。

### 6. 验证精确度分母

```bash
for caller in sawfish sniffles2 pbsv cutesv; do
  echo "--- ${caller} ---"
  python3 -c "
import json
d = json.load(open('truvari_cmrG/${caller}/summary.json'))
tp = d['TP-comp']
fp = d['FP']
prec = d['precision']
calc = tp / (tp + fp) if (tp + fp) > 0 else 0
print(f'TP={tp} FP={fp} TP+FP={tp+fp} recorded={prec} calculated={calc} match={abs(prec-calc)<1e-10}')
"
done
```

## 参考信息

| 项目 | 值 |
|------|-----|
| 参考 FASTA | human_GRCh38_no_alt_analysis_set.fasta |
| 参考构建 | GRCh38 |
| 参考 SHA256 | 9cce8b926416dd96b152deea85188495b75f7ac8d634cc723a017067be8702b7 |

## 工具版本

| 工具 | 版本 |
|------|------|
| DeepVariant | 1.10.0 |
| hap.py | 0.3.15 |
| Truvari | v5.4.0 |
| Sawfish | 2.2.1 |
| Sniffles2 | 2.8.0 |
| pbsv | 2.11.0 |
| cuteSV | 2.1.4 |

## 可重复性的局限性

- hap.py 环境使用 Python 2.7.15，这是传统环境。未来版本的 hap.py 可能产生不同的结果。
- 外部工具版本、容器镜像、参考基因组、数据库和缓存不与 HiFiVar 捆绑。不同的版本可能产生不同的结果。
- DeepVariant 运行使用 4 个分片，因资源限制。不同的分片数可能影响运行时间但不影响科学结果。
- CMRG 结果仅覆盖具有挑战性的医学相关基因，不是全基因组 SV 性能。
