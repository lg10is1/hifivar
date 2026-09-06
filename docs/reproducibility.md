# Reproducibility

This document describes how to reproduce the HG002 benchmark results
and verify the integrity of all benchmark artifacts.

## Frozen evidence files

| File | Description |
|------|-------------|
| `HG002_BENCHMARK_FROZEN_SUMMARY.tsv` | Machine-readable summary table |
| `HG002_BENCHMARK_FROZEN_METADATA.yaml` | Full metadata with checksums |
| `HG002_BENCHMARK_FROZEN_CHECKSUMS.sha256` | All artifact checksums |
| `HG002_BENCHMARK_INTERPRETATION.md` | Interpretation guidelines |

## Verification steps

### 1. Verify checksums

```bash
sha256sum -c HG002_BENCHMARK_FROZEN_CHECKSUMS.sha256
```

### 2. Verify small-variant VCF

```bash
bcftools view -h HG002_ccs11kb.wgs.small.vcf.gz | head -20
bcftools query -l HG002_ccs11kb.wgs.small.vcf.gz
```

Expected: sample HG002, GRCh38 contigs (chr1-chr22).

### 3. Verify CMRG truth VCF

```bash
bcftools view -h HG002_GRCh38_CMRG_SV_v1.00.vcf.gz | head -20
bcftools view -H HG002_GRCh38_CMRG_SV_v1.00.vcf.gz | wc -l
```

Expected: 250 records, GRCh38 contigs, sample HG002.

### 4. Verify caller VCFs

```bash
for caller in sawfish sniffles2 pbsv cutesv; do
  echo "--- ${caller} ---"
  bcftools view -H HG002.${caller}.sv.vcf.gz | wc -l
  bcftools query -l HG002.${caller}.sv.vcf.gz
done
```

### 5. Verify Truvari results

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

All callers should show TP+FN = 216 = base_cnt.

### 6. Verify precision denominators

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

## Reference information

| Item | Value |
|------|-------|
| Reference FASTA | human_GRCh38_no_alt_analysis_set.fasta |
| Reference build | GRCh38 |
| Reference SHA256 | 9cce8b926416dd96b152deea85188495b75f7ac8d634cc723a017067be8702b7 |

## Tool versions

| Tool | Version |
|------|---------|
| DeepVariant | 1.10.0 |
| hap.py | 0.3.15 |
| Truvari | v5.4.0 |
| Sawfish | 2.2.1 |
| Sniffles2 | 2.8.0 |
| pbsv | 2.11.0 |
| cuteSV | 2.1.4 |

## Limitations of reproducibility

- The hap.py environment uses Python 2.7.15, which is a legacy
  environment. Future versions of hap.py may produce different results.
- External tool versions, container images, references, databases, and
  caches are not bundled with HiFiVar. Different versions may yield
  different results.
- The DeepVariant run used 4 shards due to resource limits. Different
  shard counts may affect runtime but not scientific results.
- CMRG results cover challenging medically relevant genes only and are
  not whole-genome SV performance.
