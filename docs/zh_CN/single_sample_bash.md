# HiFiVar 0.1.0rc2 单样本 Bash 运行指南

本文给出一个可复制的 Linux 单样本运行流程，使用一个未比对的 PacBio
HiFi uBAM 和 GRCh38 reference。它适用于：

- 已经获得计算节点的交互作业；
- 独占的 Linux 计算服务器；
- 用 `sbatch` 申请到一个节点后，在该节点内部以普通 Bash 运行。

不要在共享 HPC 登录节点直接执行全基因组比对或变异检测。

## 1. 当前 RC2 的输入边界

公开版 `0.1.0rc2` 的 Snakemake calling DAG 接受已经排序并建立索引的
BAM/CRAM，不直接把 uBAM 接入 calling DAG。因此，单样本完整路径是：

```text
PacBio HiFi uBAM
  -> pbmm2 align --preset CCS --sort
  -> coordinate-sorted BAM + BAI
  -> HiFiVar sample sheet
  -> config validation + Snakemake dry-run
  -> 单样本 calling tracks
  -> 可选 review / annotation / benchmark
```

uBAM 是未比对 BAM，不能因为扩展名为 `.bam` 就当成已有比对结果直接交给
DeepVariant、SV caller 或 TRGT。

## 2. 本例使用的数据

本例先运行 `SAMPLE_A`：

```bash
export SAMPLE_ID="SAMPLE_A"
export UBAM="/data/project/pacbio_smrtlink/data/example_data/batch_1/hifi_reads/sample_a.hifi_reads.bam"
export REFERENCE_FASTA="/data/project/references/GRCh38.fa"
export RUN_ROOT="/work/hifivar/realdata_validation/single_${SAMPLE_ID}"
```

`RUN_ROOT` 必须是新的、专用于本次运行的目录。不要把原始数据放入可清理的
`work/` 目录。

## 3. 需要由用户填写的部署参数

以下值不能由 HiFiVar 猜测：

```bash
export CONDA_ROOT="<CONDA或MINIFORGE安装目录>"
export PBMM2_ENV="<已经真实验证的pbmm2环境名>"
export HIFIVAR_ENV="hifivar-rc2-public-test"
export DEEPVARIANT_IMAGE="<DeepVariant-1.10.0的已验证SIF绝对路径>"
export THREADS=16
export SAMTOOLS_INDEX_THREADS=4
```

示例中的尖括号必须替换。不要把 `latest` 当作正式版本，也不要在未验证的环境中
一次性开启全部 caller。

## 4. 创建运行目录

```bash
set -euo pipefail

mkdir -p "$RUN_ROOT"/{configs,alignments,logs/alignment,logs/workflow,work,tmp,results,validation}
cd "$RUN_ROOT"
```

推荐结构：

```text
single_SAMPLE_A/
├── configs/
├── alignments/
├── logs/
│   ├── alignment/
│   └── workflow/
├── results/
├── tmp/
├── validation/
└── work/
```

## 5. 输入与 reference 预检查

```bash
test -r "$UBAM"
test -s "$UBAM"
test -r "$REFERENCE_FASTA"
test -s "$REFERENCE_FASTA"
test -s "${REFERENCE_FASTA}.fai"

source "$CONDA_ROOT/etc/profile.d/conda.sh"
conda activate "$PBMM2_ENV"

command -v pbmm2
command -v samtools
pbmm2 --version
samtools --version | head -1
samtools quickcheck -v "$UBAM"
pbmm2 align --help >"$RUN_ROOT/validation/pbmm2_align_help.txt" 2>&1
```

必须根据服务器上真实 `pbmm2 align --help` 核对参数。reference 的 contig 命名、
BAM header 和所有下游资源必须兼容；HiFiVar 不会静默执行 `chr1`/`1` 转换。

## 6. 将 uBAM 比对为排序 BAM

```bash
export ALIGNED_BAM="$RUN_ROOT/alignments/${SAMPLE_ID}.aligned.bam"
export ALIGNED_BAI="${ALIGNED_BAM}.bai"
export ALIGNMENT_DONE="$RUN_ROOT/validation/${SAMPLE_ID}.alignment.complete.tsv"
```

首次运行前拒绝覆盖已有产物：

```bash
if [[ -e "$ALIGNED_BAM" || -e "$ALIGNED_BAI" || -e "$ALIGNMENT_DONE" ]]; then
  echo "发现已有比对产物；请先核验，不要静默覆盖：$ALIGNED_BAM" >&2
  exit 13
fi
```

执行 pbmm2：

```bash
READ_GROUP=$'@RG\tID:'"$SAMPLE_ID"$'\tSM:'"$SAMPLE_ID"$'\tPL:PACBIO'

pbmm2 align \
  "$REFERENCE_FASTA" \
  "$UBAM" \
  "$ALIGNED_BAM" \
  --preset CCS \
  --sort \
  --bam-index NONE \
  --rg "$READ_GROUP" \
  -j "$THREADS" \
  --log-level INFO \
  >"$RUN_ROOT/logs/alignment/${SAMPLE_ID}.pbmm2.stdout.log" \
  2>"$RUN_ROOT/logs/alignment/${SAMPLE_ID}.pbmm2.stderr.log"

samtools quickcheck -v "$ALIGNED_BAM"
samtools index -@ "$SAMTOOLS_INDEX_THREADS" "$ALIGNED_BAM"
test -s "$ALIGNED_BAI"
```

验证 header、排序状态、样本名和至少一个 mapped read：

```bash
samtools view -H "$ALIGNED_BAM" >"$RUN_ROOT/validation/${SAMPLE_ID}.header.sam"
grep -q '^@SQ' "$RUN_ROOT/validation/${SAMPLE_ID}.header.sam"
grep -q 'SO:coordinate' "$RUN_ROOT/validation/${SAMPLE_ID}.header.sam"
grep -q $'SM:'"$SAMPLE_ID" "$RUN_ROOT/validation/${SAMPLE_ID}.header.sam"

MAPPED_READS="$(samtools view -c -F 4 "$ALIGNED_BAM")"
if [[ "$MAPPED_READS" -le 0 ]]; then
  echo "比对结果中没有 mapped reads" >&2
  exit 14
fi

printf 'sample_id\tbam\tindex\tmapped_reads\n%s\t%s\t%s\t%s\n' \
  "$SAMPLE_ID" "$ALIGNED_BAM" "$ALIGNED_BAI" "$MAPPED_READS" \
  >"$ALIGNMENT_DONE"
```

只有 `ALIGNMENT_DONE` 存在且上述检查全部成功，才能进入 calling。

## 7. 创建单样本表

创建 `$RUN_ROOT/configs/samples.tsv`。这是一个真正的 Tab 分隔文件：

```bash
printf 'sample_id\tinput\tinput_type\n%s\t%s\tbam\n' \
  "$SAMPLE_ID" "$ALIGNED_BAM" \
  >"$RUN_ROOT/configs/samples.tsv"
```

检查：

```bash
column -t -s $'\t' "$RUN_ROOT/configs/samples.tsv"
```

## 8. 创建最小单样本 calling 配置

先只开启真实部署已经验证的 DeepVariant。创建
`$RUN_ROOT/configs/calling.yaml`：

```yaml
project:
  name: example_sample_a_rc2

reference:
  fasta: /data/project/references/GRCh38.fa
  build: GRCh38

samples:
  sheet: <RUN_ROOT>/configs/samples.tsv

runtime:
  threads: 16
  tmpdir: <RUN_ROOT>/tmp

paths:
  workdir: <RUN_ROOT>/work/calling
  outdir: <RUN_ROOT>/results/calling

logging:
  level: INFO
  file: <RUN_ROOT>/logs/workflow/hifivar.log

small:
  enabled: true
  execution_mode: apptainer
  deepvariant_image: <VALIDATED_DEEPVARIANT_1.10_IMAGE>
  model_type: PACBIO
  threads: 16
  memory_mb: 64000
  runtime_minutes: 2880
  overwrite: false

sv:
  enabled: false
tr:
  enabled: false
phasing:
  enabled: false
assembly:
  enabled: false
assembly_sv:
  enabled: false
review:
  enabled: false
annotation:
  enabled: false
cohort:
  enabled: false
benchmark:
  enabled: false
```

把所有 `<RUN_ROOT>` 和 image 占位符替换成绝对路径。

## 9. 验证配置并生成 effective config

```bash
conda activate "$HIFIVAR_ENV"
unset PYTHONPATH PYTHONHOME

hifivar --version
hifivar doctor
snakemake --version
apptainer --version

hifivar --config "$RUN_ROOT/configs/calling.yaml" config validate
hifivar --config "$RUN_ROOT/configs/calling.yaml" config dump-effective \
  --output "$RUN_ROOT/configs/calling.effective.yaml"

WORKFLOW_ROOT="$(python -c 'from hifivar.package_resources import installed_workflow_root; print(installed_workflow_root())')"
test -f "$WORKFLOW_ROOT/Snakefile"
```

DeepVariant 分片运行前检查文件描述符：

```bash
ulimit -n
if [[ "$(ulimit -n)" -lt 4096 ]]; then
  echo '文件描述符上限不足；DeepVariant至少需要4096，建议65536。' >&2
  exit 20
fi
```

## 10. 必须先 dry-run

```bash
snakemake \
  --snakefile "$WORKFLOW_ROOT/Snakefile" \
  --configfile "$RUN_ROOT/configs/calling.effective.yaml" \
  --cores 1 \
  --dry-run \
  --printshellcmds \
  2>&1 | tee "$RUN_ROOT/logs/workflow/snakemake.dry-run.log"
```

检查 DAG 中只有预期样本和已启用分支。出现路径、reference、image、sample ID
问题时应先修配置并重新生成 effective config，不要直接正式运行。

## 11. 在一个计算节点上正式执行

```bash
cd "$RUN_ROOT"

snakemake \
  --snakefile "$WORKFLOW_ROOT/Snakefile" \
  --configfile "$RUN_ROOT/configs/calling.effective.yaml" \
  --cores "$THREADS" \
  --resources mem_mb=64000 deepvariant_slots=1 \
  --rerun-incomplete \
  --printshellcmds \
  2>&1 | tee "$RUN_ROOT/logs/workflow/snakemake.controller.log"
```

`--cores 16` 表示当前这一台节点最多使用 16 核，不表示向 Slurm 提交 16 个节点。
首次真实运行不建议加 `--keep-going`。

如果在普通 Linux 服务器后台运行，可在站点允许时使用：

```bash
nohup snakemake \
  --snakefile "$WORKFLOW_ROOT/Snakefile" \
  --configfile "$RUN_ROOT/configs/calling.effective.yaml" \
  --cores "$THREADS" \
  --resources mem_mb=64000 deepvariant_slots=1 \
  --rerun-incomplete \
  --printshellcmds \
  >"$RUN_ROOT/logs/workflow/nohup.stdout.log" \
  2>"$RUN_ROOT/logs/workflow/nohup.stderr.log" &
```

在 HPC 上更推荐把前一个正式命令放进 `sbatch`，或使用多样本 Slurm 指南，
而不是依靠登录会话中的 `nohup`。

## 12. 运行后检查

```bash
snakemake \
  --snakefile "$WORKFLOW_ROOT/Snakefile" \
  --configfile "$RUN_ROOT/configs/calling.effective.yaml" \
  --summary

find "$RUN_ROOT/results" -type f -size +0c -print | sort
find "$RUN_ROOT/logs" -type f -size +0c -print | sort
```

单样本 small-variant 结果应保持独立：

```text
<sample>.small.vcf.gz
<sample>.small.vcf.gz.tbi
<sample>.small.g.vcf.gz
<sample>.small.g.vcf.gz.tbi
```

具体目录以 dry-run 和 effective config 为准。不得用文件“存在”代替 VCF header、
index、sample、contig 和 provenance 校验。

## 13. 逐步增加其他分支

推荐顺序：

1. 单样本 pbmm2 + DeepVariant；
2. 只增加一个已经真实验证的 read-SV caller；
3. 再逐个增加其他 read-SV caller；
4. 有兼容 GRCh38 catalog 时才开启 TRGT；
5. 再执行 phasing、assembly、assembly-SV；
6. review、annotation、benchmark 作为下游可选分支。

每增加一个分支都必须：修改 user config、重新 `config validate`、重新生成
effective config、重新 dry-run，然后才执行。

## 14. 失败与恢复

- 保留原始 uBAM、reference、BAM、VCF/gVCF、日志和 quarantine 产物。
- 不要用删除完好结果的方式“从头再来”。
- 外部工具成功但输出校验失败时，应先检查 quarantine 和 stdout/stderr。
- 修复运行目录中的配置、环境或资源后，再使用 `--rerun-incomplete`。
- `FAILED`、`NOT_RUN`、`MISSING_INPUT` 不能解释为阴性生物学结果。
- 不要让两个 Snakemake controller 同时写同一个 `RUN_ROOT`。

## 15. 完成标准

单样本运行只有在以下条件全部满足后才算完成：

- uBAM 可读取且通过 `samtools quickcheck`；
- pbmm2 输出为 coordinate-sorted BAM；
- BAM index、`@SQ`、sample read group 和 mapped reads 检查通过；
- HiFiVar config validation 和 dry-run 通过；
- 所有已启用规则成功；
- 预期 VCF/gVCF 及 index 通过校验；
- stdout、stderr、effective config、工具版本和 provenance 已保留；
- 未启用的分支明确为 `DISABLED` 或 `NOT_RUN`。
