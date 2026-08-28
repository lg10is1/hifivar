# HiFiVar 0.1.0rc2 多样本 Slurm 运行指南

本文说明如何在 Slurm 集群上批量处理多个 PacBio HiFi uBAM。详细英文部署稿见
[`../slurm_multi_sample.md`](../slurm_multi_sample.md)。本文是可执行的中文流程，
不是把三个单样本脚本手工复制三遍。

## 1. 推荐调度模型

```text
ubam_samples.tsv
  -> Slurm array：每个 uBAM 一个 pbmm2 比对任务
  -> aligned_samples.tsv
  -> Snakemake + Slurm executor
       -> 每个样本的 DeepVariant
       -> 每个样本的已启用 read-SV caller
       -> 每个样本的 TRGT（需要统一 catalog）
  -> 所有单样本结果通过验证
  -> GLnexus / SV cohort / TR cohort
```

不要为每个样本手写一套 calling 命令。sample sheet 决定样本集合，Snakemake
负责展开 `sample × rule`，Slurm 负责资源分配。

## 2. 三个真实 uBAM 与 reference

```text
SAMPLE_A  /data/project/pacbio_smrtlink/data/example_data/batch_1/hifi_reads/sample_a.hifi_reads.bam
SAMPLE_B  /data/project/pacbio_smrtlink/data/example_data/batch_1/hifi_reads/sample_b.hifi_reads.bam
SAMPLE_C  /data/project/pacbio_smrtlink/data/example_data/batch_2/hifi_reads/sample_c.hifi_reads.bam
```

Reference：

```text
/data/project/references/GRCh38.fa
```

这些 BAM 是未比对 uBAM，必须先比对，不能直接作为 calling 输入。

## 3. 设置运行根目录

```bash
export RUN_ROOT="/work/hifivar/realdata_validation/example_3sample_rc2"
export REFERENCE_FASTA="/data/project/references/GRCh38.fa"

mkdir -p "$RUN_ROOT"/{configs/slurm,slurm,alignments,logs/alignment,logs/workflow,work,tmp,results,validation}
cd "$RUN_ROOT"
```

需要从集群管理员或现有作业脚本确认：

```text
<SLURM_PARTITION>
<SLURM_ACCOUNT>
<CONDA_ROOT>
<VALIDATED_PBMM2_ENV>
<VALIDATED_DEEPVARIANT_1.10_IMAGE>
```

## 4. 创建 uBAM manifest

创建 `$RUN_ROOT/configs/ubam_samples.tsv`，使用真正的 Tab：

```text
sample_id	ubam
SAMPLE_A	/data/project/pacbio_smrtlink/data/example_data/batch_1/hifi_reads/sample_a.hifi_reads.bam
SAMPLE_B	/data/project/pacbio_smrtlink/data/example_data/batch_1/hifi_reads/sample_b.hifi_reads.bam
SAMPLE_C	/data/project/pacbio_smrtlink/data/example_data/batch_2/hifi_reads/sample_c.hifi_reads.bam
```

检查行数和文件：

```bash
awk -F '\t' 'NR>1 {print NR-1, $1, $2}' "$RUN_ROOT/configs/ubam_samples.tsv"
while IFS=$'\t' read -r sample ubam; do
  [[ "$sample" == "sample_id" ]] && continue
  test -s "$ubam"
done <"$RUN_ROOT/configs/ubam_samples.tsv"

test -s "$REFERENCE_FASTA"
test -s "${REFERENCE_FASTA}.fai"
```

## 5. 创建 pbmm2 Slurm array 脚本

保存为 `$RUN_ROOT/slurm/align_ubam_array.sbatch`：

```bash
#!/usr/bin/env bash
#SBATCH --job-name=hifivar-align
#SBATCH --partition=<SLURM_PARTITION>
#SBATCH --account=<SLURM_ACCOUNT>
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/alignment/%A_%a.stdout.log
#SBATCH --error=logs/alignment/%A_%a.stderr.log

set -euo pipefail

: "${RUN_ROOT:?RUN_ROOT未传入}"
: "${REFERENCE_FASTA:?REFERENCE_FASTA未传入}"

source <CONDA_ROOT>/etc/profile.d/conda.sh
conda activate <VALIDATED_PBMM2_ENV>

command -v pbmm2 >/dev/null
command -v samtools >/dev/null
pbmm2 --version
samtools --version | head -1

manifest="$RUN_ROOT/configs/ubam_samples.tsv"
row="$(awk -F '\t' -v task="$SLURM_ARRAY_TASK_ID" 'NR == task + 1 {print; exit}' "$manifest")"
if [[ -z "$row" ]]; then
  echo "array task没有对应manifest行：$SLURM_ARRAY_TASK_ID" >&2
  exit 10
fi

IFS=$'\t' read -r sample_id ubam <<<"$row"
if [[ -z "$sample_id" || -z "$ubam" ]]; then
  echo "manifest行格式错误：$row" >&2
  exit 11
fi
if [[ ! -r "$ubam" || ! -s "$ubam" ]]; then
  echo "uBAM不可读或为空：$ubam" >&2
  exit 12
fi

output="$RUN_ROOT/alignments/${sample_id}.aligned.bam"
index="${output}.bai"
done_file="$RUN_ROOT/validation/${sample_id}.alignment.complete.tsv"
if [[ -e "$output" || -e "$index" || -e "$done_file" ]]; then
  echo "拒绝覆盖已有输出：$output" >&2
  exit 13
fi

samtools quickcheck -v "$ubam"
read_group=$'@RG\tID:'"$sample_id"$'\tSM:'"$sample_id"$'\tPL:PACBIO'

pbmm2 align \
  "$REFERENCE_FASTA" \
  "$ubam" \
  "$output" \
  --preset CCS \
  --sort \
  --bam-index NONE \
  --rg "$read_group" \
  -j "$SLURM_CPUS_PER_TASK" \
  --log-level INFO

samtools quickcheck -v "$output"
samtools index -@ 4 "$output"
test -s "$index"

header="$(samtools view -H "$output")"
grep -q '^@SQ' <<<"$header"
grep -q 'SO:coordinate' <<<"$header"
grep -q $'SM:'"$sample_id" <<<"$header"

mapped="$(samtools view -c -F 4 "$output")"
if [[ "$mapped" -le 0 ]]; then
  echo "样本没有mapped reads：$sample_id" >&2
  exit 14
fi

printf 'sample_id\tbam\tindex\tmapped_reads\n%s\t%s\t%s\t%s\n' \
  "$sample_id" "$output" "$index" "$mapped" >"$done_file"
```

提交前先核对真实版本和 CLI：

```bash
source <CONDA_ROOT>/etc/profile.d/conda.sh
conda activate <VALIDATED_PBMM2_ENV>
pbmm2 --version
pbmm2 align --help
samtools --version | head -1
```

## 6. 先跑一个 array task

从 `RUN_ROOT` 提交，因为日志路径是相对提交目录：

```bash
cd "$RUN_ROOT"

TEST_JOB_ID="$(sbatch \
  --parsable \
  --array=1-1%1 \
  --export=ALL,RUN_ROOT="$RUN_ROOT",REFERENCE_FASTA="$REFERENCE_FASTA" \
  slurm/align_ubam_array.sbatch)"

echo "$TEST_JOB_ID"
squeue -j "$TEST_JOB_ID"
```

`1-1` 只运行 `SAMPLE_A`。确认完成记录、BAM、BAI、日志都通过后，再提交剩余
样本。因为脚本拒绝覆盖，已经完成的第一项不要再次提交；可运行：

```bash
REST_JOB_ID="$(sbatch \
  --parsable \
  --array=2-3%2 \
  --export=ALL,RUN_ROOT="$RUN_ROOT",REFERENCE_FASTA="$REFERENCE_FASTA" \
  slurm/align_ubam_array.sbatch)"
```

如果在全新目录中一次提交三个样本：

```bash
ALIGN_JOB_ID="$(sbatch \
  --parsable \
  --array=1-3%2 \
  --export=ALL,RUN_ROOT="$RUN_ROOT",REFERENCE_FASTA="$REFERENCE_FASTA" \
  slurm/align_ubam_array.sbatch)"
```

`1-3%2` 表示共三个 task，同时最多运行两个。不要只为了“更快”就提高并发；pbmm2
会同时消耗 CPU、内存和共享文件系统带宽。

## 7. 监控和验收比对

```bash
JOB_ID="<本次实际提交返回的job ID>"
squeue -j "$JOB_ID" -o '%.18i %.9P %.28j %.8T %.10M %.9l %.6D %R'
sacct -j "$JOB_ID" \
  --format=JobID,JobName%28,State,ExitCode,Elapsed,AllocCPUS,ReqMem,MaxRSS,NodeList
```

三个任务必须全部 `COMPLETED` 且 `ExitCode=0:0`：

```bash
for sample in SAMPLE_A SAMPLE_B SAMPLE_C; do
  test -s "$RUN_ROOT/validation/${sample}.alignment.complete.tsv"
  test -s "$RUN_ROOT/alignments/${sample}.aligned.bam"
  test -s "$RUN_ROOT/alignments/${sample}.aligned.bam.bai"
  samtools quickcheck -v "$RUN_ROOT/alignments/${sample}.aligned.bam"
done
```

任何 task 为 `FAILED`、`CANCELLED`、`TIMEOUT` 或 `OUT_OF_MEMORY` 时，都不要启动
calling。

## 8. 创建 HiFiVar aligned sample sheet

创建 `$RUN_ROOT/configs/aligned_samples.tsv`：

```text
sample_id	input	input_type
SAMPLE_A	<RUN_ROOT>/alignments/SAMPLE_A.aligned.bam	bam
SAMPLE_B	<RUN_ROOT>/alignments/SAMPLE_B.aligned.bam	bam
SAMPLE_C	<RUN_ROOT>/alignments/SAMPLE_C.aligned.bam	bam
```

把 `<RUN_ROOT>` 替换为绝对路径。

## 9. 创建三样本 calling 配置

首先只开启 DeepVariant：

```yaml
project:
  name: example_3sample_rc2

reference:
  fasta: /data/project/references/GRCh38.fa
  build: GRCh38

samples:
  sheet: <RUN_ROOT>/configs/aligned_samples.tsv

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

## 10. 配置验证和 DAG dry-run

```bash
conda activate hifivar-rc2-public-test
unset PYTHONPATH PYTHONHOME

hifivar --config "$RUN_ROOT/configs/calling.yaml" config validate
hifivar --config "$RUN_ROOT/configs/calling.yaml" config dump-effective \
  --output "$RUN_ROOT/configs/calling.effective.yaml"

WORKFLOW_ROOT="$(python -c 'from hifivar.package_resources import installed_workflow_root; print(installed_workflow_root())')"

snakemake \
  --snakefile "$WORKFLOW_ROOT/Snakefile" \
  --configfile "$RUN_ROOT/configs/calling.effective.yaml" \
  --cores 1 \
  --dry-run \
  --printshellcmds
```

三行 sample sheet 会展开为三个 `deepvariant_small` 样本任务。以后启用四个 read-SV
caller 时，会产生最多十二个样本/caller 任务，但实际依赖关系由 DAG 决定。

## 11. 准备 Snakemake Slurm executor

不要直接污染已经验证的 RC2 环境；建议克隆：

```bash
mamba create -n hifivar-rc2-slurm --clone hifivar-rc2-public-test
conda activate hifivar-rc2-slurm
python -m pip install 'snakemake-executor-plugin-slurm==2.7.1'

snakemake --version
python -m pip show snakemake-executor-plugin-slurm
snakemake --help | grep -A2 -B2 slurm
```

`2.7.1` 是本部署文档固定的候选版本，不表示已经在你的 Slurm 站点完成验证。

## 12. 创建 Slurm profile

保存为 `$RUN_ROOT/configs/slurm/profile.v9+.yaml`：

```yaml
executor: slurm
jobs: 20
latency-wait: 120
rerun-incomplete: true
keep-going: false
printshellcmds: true
show-failed-logs: true

default-resources:
  slurm_account: "<SLURM_ACCOUNT>"
  slurm_partition: "<SLURM_PARTITION>"
  mem_mb: 4000
  runtime: 60

set-resources:
  deepvariant_small:
    runtime: 2880
  read_based_sv:
    runtime: 1440
  tandem_repeat:
    runtime: 720
  cohort_small_variants:
    runtime: 1440
  cohort_sv:
    runtime: 240
  cohort_tr:
    runtime: 240
```

关键说明：

- `jobs: 20` 是同时处于提交/运行状态的 Snakemake 作业上限，不是 CPU 数；
- rule 中的 `threads` 对应每个作业的 CPU；
- `mem_mb` 对应作业总内存；
- HiFiVar 根据 `small.max_concurrent_samples: 1` 注册全局
  `deepvariant_slots`，默认整个 DAG 一次只允许一个 DeepVariant 样本运行；只有在
  tiny DAG 证实每个样本使用独立 TMPDIR 且节点资源足够后，才调高该配置值；
- 当前 HiFiVar config 记录 `runtime_minutes`，rule 记录 `runtime_min`；
- Slurm executor 使用标准资源名 `runtime` 映射 `sbatch --time`，所以 profile 中需要
  显式设置 `runtime`，不能假设它会自动翻译 `runtime_min`。

## 13. 使用 Slurm executor 正式运行

先 dry-run：

```bash
conda activate hifivar-rc2-slurm
unset PYTHONPATH PYTHONHOME

snakemake \
  --snakefile "$WORKFLOW_ROOT/Snakefile" \
  --configfile "$RUN_ROOT/configs/calling.effective.yaml" \
  --profile "$RUN_ROOT/configs/slurm" \
  --dry-run \
  --printshellcmds
```

从登录/提交节点启动 controller：

```bash
snakemake \
  --snakefile "$WORKFLOW_ROOT/Snakefile" \
  --configfile "$RUN_ROOT/configs/calling.effective.yaml" \
  --profile "$RUN_ROOT/configs/slurm" \
  --printshellcmds \
  2>&1 | tee "$RUN_ROOT/logs/workflow/snakemake.controller.log"
```

controller 必须保持存活以监控作业。站点允许时可使用 `tmux`/`screen`。不要在同一
输出目录启动第二个 controller。

## 14. 单节点 Slurm 回退方案

如果暂时无法安装 Slurm executor plugin，可以申请一个大节点，在节点内部运行本地
Snakemake。保存为 `$RUN_ROOT/slurm/run_calling_single_node.sbatch`：

```bash
#!/usr/bin/env bash
#SBATCH --job-name=hifivar-calling
#SBATCH --partition=<SLURM_PARTITION>
#SBATCH --account=<SLURM_ACCOUNT>
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=3-00:00:00
#SBATCH --output=logs/workflow/%j.stdout.log
#SBATCH --error=logs/workflow/%j.stderr.log

set -euo pipefail
source <CONDA_ROOT>/etc/profile.d/conda.sh
conda activate hifivar-rc2-public-test
unset PYTHONPATH PYTHONHOME

: "${RUN_ROOT:?RUN_ROOT未传入}"

ulimit -n 65536 || true
if [[ "$(ulimit -n)" -lt 4096 ]]; then
  echo 'DeepVariant文件描述符上限不足，至少4096，建议65536。' >&2
  exit 20
fi

WORKFLOW_ROOT="$(python -c 'from hifivar.package_resources import installed_workflow_root; print(installed_workflow_root())')"

snakemake \
  --snakefile "$WORKFLOW_ROOT/Snakefile" \
  --configfile "$RUN_ROOT/configs/calling.effective.yaml" \
  --cores "$SLURM_CPUS_PER_TASK" \
  --resources mem_mb=120000 deepvariant_slots=1 \
  --rerun-incomplete \
  --printshellcmds
```

提交：

```bash
cd "$RUN_ROOT"
sbatch --export=ALL,RUN_ROOT="$RUN_ROOT" slurm/run_calling_single_node.sbatch
```

此模式只使用一个节点。`--cores 32` 不等于 32 个 Slurm 作业。

## 15. 初始资源建议

以下只是 RC2 起点，应根据 `sacct`、覆盖度、日志和集群限制调整：

| 阶段 | 三样本任务数 | 每任务CPU | 每任务内存 | 初始时限 |
|---|---:|---:|---:|---:|
| pbmm2 | 3 | 16 | 48 GB | 48 h |
| DeepVariant | 3 | 16 | 64 GB | 48 h |
| Sawfish | 3 | 16 | 32 GB | 24 h |
| Sniffles2 | 3 | 8 | 16 GB | 12 h |
| pbsv | 3 | 8 | 32 GB | 24 h |
| cuteSV | 3 | 8 | 16 GB | 12 h |
| TRGT | 3 | 8 | 16 GB | 12 h |
| GLnexus cohort | 1 | 8 | 已验证三样本 WGS 工作负载建议 192–200 GB | 24 h |

## 16. 多样本 cohort handoff

只有三个样本的对应单样本 track 全部通过验证后，才能创建 cohort manifest。表头必须
为：

```text
sample	track	state	source_path	index_path	source_tool	source_version	reference_build	catalog_id
```

- small-variant 行使用每个样本的 DeepVariant gVCF 和 TBI；
- `state` 使用明确的 `CALLED`/`FAILED`/`NOT_RUN` 等状态；
- SV 优先使用每个样本 harmonized SV VCF；
- TR 必须使用相同 catalog ID；
- 缺失结果不能伪装成阴性样本。

单独创建 cohort config：

```yaml
cohort:
  enabled: true
  cohort_id: EXAMPLE_3SAMPLE_RC2
  input_manifest: <RUN_ROOT>/configs/cohort_inputs.tsv
  overwrite: false
  small_variants:
    enabled: true
    glnexus_executable: glnexus_cli
    bcftools_executable: bcftools
    preset: DeepVariantWGS
    threads: 8
    # 必须显式设置。已验证三样本 WGS 峰值约 153 GB，此处 192 GB 留出余量。
    memory_gb: 192
    runtime_minutes: 1440
  sv:
    enabled: false
  tr:
    enabled: false
```

GLnexus 应使用已经验证的 `hifivar-glnexus-1.4.1` 环境或等价的固定版本部署。
192 GB 不是普适默认值；扩大 cohort 前应根据样本数、gVCF 密度、站点限制和峰值 RSS
重新估算。

## 17. 推荐推进顺序

1. 对 `SAMPLE_A` 做小子集 smoke；
2. 完整运行 `SAMPLE_A` 的 pbmm2 + DeepVariant；
3. 为该样本逐个增加 read-SV caller；
4. 有兼容 catalog 后增加 TRGT；
5. 批量比对并分析 `SAMPLE_B`、`SAMPLE_C`；
6. 所有单样本 gVCF 通过校验后运行 GLnexus；
7. 上游契约完整后才运行 SV/TR cohort；
8. 根据 `sacct` 调整并发、内存和时限，再扩展更多样本。

## 18. 失败处理

- 查看 rule 日志和 `sacct` 的 State、ExitCode、MaxRSS；
- 只修复运行目录内的配置、路径、环境或资源问题；
- 不删除原始 uBAM、reference、BAM、VCF/gVCF 或 quarantine；
- 理解失败原因后再使用 `--rerun-incomplete`；
- 首次验收不要启用 `--keep-going`；
- 不让多个 controller 写同一结果树；
- 工具退出成功但产物校验失败时，仍属于失败，必须保存现场。

## 19. 完成标准

多样本运行完成必须同时满足：

- 所有 array task 为 `COMPLETED`、`ExitCode=0:0`；
- 三个 aligned BAM/BAI 和 completion record 全部存在并通过检查；
- effective config 和 dry-run 已保存；
- 每个启用的样本/track 都有明确状态、日志、工具版本和 provenance；
- cohort 只消费验证通过且 reference/catalog 兼容的输入；
- 未运行或失败的 track 没有被解释为阴性生物学结果。
