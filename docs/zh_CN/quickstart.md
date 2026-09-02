# HiFiVar 0.1.0rc4 完整中文快速入门

[English](../quickstart.md) | [简体中文](quickstart.md)

本指南说明 Linux/HPC 上从安装、输入准备、配置、dry-run 到单样本和多样本运行
的完整路径。真实 WGS 前必须先用 tiny 数据验证。

## 1. 先理解公开版本边界

HiFiVar 采用“core + 独立外部工具”模式。wheel 不包含 caller、reference、
database、catalog 或 WGS 数据。打包的 small/SV/TR 工作流以已经比对并建立索引
的 BAM/CRAM 为输入。`0.1.0rc4` 没有统一的 `hifivar run`，也没有打包的
Snakemake alignment rule；原始 FASTQ/uBAM 需要先显式比对。hifiasm assembly
分支则使用 HiFi FASTQ。

## 2. 安装与验证

按照[中文安装指南](installation.md)下载 Release 三个文件、执行 SHA256 校验并
安装 wheel。随后运行：

```bash
hifivar --version
hifivar --help
hifivar config validate
hifivar doctor
```

预期输出版本为 `hifivar 0.1.0rc4`。

## 3. 建立项目目录

```text
analysis/
├── config.yaml
├── effective_config.yaml
├── samples.tsv
├── references/
├── containers/
├── databases/
├── work/
├── results/
└── logs/
```

原始 FASTQ/uBAM、正式 BAM/CRAM 和 reference 不要放进可随时清理的 `work/`。

## 4. 准备 reference

使用未压缩 FASTA，并准备 `.fai`：

```bash
samtools faidx /data/reference/GRCh38.fa
```

必须明确 reference build。检查 BAM/CRAM、VCF、TR catalog 和数据库是否采用
相同 build 与 contig 风格。HiFiVar 不会静默执行 `chr1` 与 `1` 的互换。

## 5. 准备样本表

创建 UTF-8 TSV：

```text
sample_id	input	input_type	sex
SAMPLE01	/data/alignments/SAMPLE01.bam	bam	unknown
SAMPLE02	/data/alignments/SAMPLE02.cram	cram	unknown
```

BAM 应非空、按 caller 要求进行 coordinate sort，并有 BAI；CRAM 需要匹配的
reference 与 CRAI。开始昂贵计算前检查 sample name、read group、index、contig
和 reference compatibility。

PacBio uBAM 是未比对数据，不能伪装为可复用的 aligned BAM。应先通过显式
pbmm2 步骤生成 sorted/indexed BAM。assembly 分支需显式提供或转换为 HiFi FASTQ。

## 6. 生成配置

源码目录可复制 minimal 示例：

```bash
cp examples/minimal/config.yaml config.yaml
cp examples/minimal/samples.tsv samples.tsv
```

替换全部 `/data/...` 与 `/work/...` 示例路径。默认 minimal 配置不启用生物学
caller，适合先检查配置和 DAG。全局 CLI 参数必须放在子命令前：

```bash
hifivar --config config.yaml config validate
hifivar --config config.yaml config dump-effective --output effective_config.yaml
```

逐项检查并保存 `effective_config.yaml`，不要在配置中保存密码或 token。

## 7. 查找已安装 workflow 并 dry-run

```bash
WORKFLOW_ROOT="$(python -c 'from hifivar.package_resources import installed_workflow_root; print(installed_workflow_root())')"
test -f "$WORKFLOW_ROOT/Snakefile"

snakemake \
  --snakefile "$WORKFLOW_ROOT/Snakefile" \
  --configfile effective_config.yaml \
  --cores 1 \
  --dry-run \
  --printshellcmds
```

基础 dry-run 成功只说明 core 与 DAG 可加载，不表示真实 caller 已安装或科学参数
已经正确。

## 8. 一次只启用一个分支

按[工具部署矩阵](../deployment.md)准备固定版本 executable/container/database。
不要使用 `latest`。每启用一个分支，都重新生成 effective config 并再次 dry-run。

主要结果保持分离：

- small variant：DeepVariant VCF 与 gVCF；
- read-SV：四个 caller 的独立 VCF 与 harmonized evidence；
- TR：TRGT 独立 VCF；
- assembly：hifiasm 后接 PAV 与 SVIM-asm；
- cohort：GLnexus small VCF，SV/TR 使用各自 cohort artifact；
- review、annotation、benchmark、report：可选下游分支。

support count、harmonization、人工 review 和 annotation impact 都只是 evidence，
不能自动解释为 truth、call confidence 或 clinical pathogenicity。

## 9. 单样本 Bash 运行

dry-run 通过后：

```bash
snakemake \
  --snakefile "$WORKFLOW_ROOT/Snakefile" \
  --configfile effective_config.yaml \
  --cores 8 \
  --printshellcmds \
  --rerun-incomplete
```

`runtime.tmpdir` 应指向容量足够的数据盘。DeepVariant 临时目录按样本隔离为
`<tmp-root>/deepvariant/<sample>/tmp`。文件描述符上限至少 4096，条件允许时建议
65536。完整模板见[单样本 Bash 运行指南](single_sample_bash.md)。

## 10. 多样本 Slurm 运行

先在登录/交互节点完成 dry-run，再使用本集群批准的 Snakemake profile：

```bash
snakemake \
  --snakefile "$WORKFLOW_ROOT/Snakefile" \
  --configfile effective_config.yaml \
  --profile /path/to/site-profile \
  --rerun-incomplete \
  --printshellcmds
```

profile 必须映射 `threads`、`mem_mb` 和 `runtime_min`。DeepVariant 并发由
`small.max_concurrent_samples` 控制，tiny 多样本测试确认 TMPDIR 隔离和资源充足
后才能提高。GLnexus 应显式申请与 cohort 规模匹配的内存。详见
[多样本 Slurm 运行指南](slurm_multi_sample.md)。

## 11. 失败恢复与结果检查

- 保存 stdout、stderr、effective config、tool version 和 container digest；
- 诊断失败后使用 `--rerun-incomplete`；
- 不删除原始 FASTQ/BAM/CRAM/VCF；
- 未理解 partial-track 行为前不要使用 `--keep-going`；
- 校验 BGZF/TBI、sample、contig、header 和 expected outputs；
- 记录 HiFiVar/Git SHA、reference checksum、database/catalog/truth version。

输出约定见[Outputs](../outputs.md)，故障诊断见[Troubleshooting](../troubleshooting.md)。

## 12. 验证范围与限制

`0.1.0rc4` 的 package 安装、资源访问、Snakemake regression 及部分真实 Linux/HPC
工具路径已通过 Release 验证。但新的集群、reference、container、database 或
scheduler profile 仍需自己的 tiny real-tool validation，之后才能运行 WGS。
HiFiVar 仅用于科研，不提供临床解释。
