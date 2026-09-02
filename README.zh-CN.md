# HiFiVar

[English](README.md) | [简体中文](README.zh-CN.md)

HiFiVar 是面向 PacBio HiFi 全基因组变异分析的模块化、可复现工作流框架。
它负责安全地编排成熟生物信息工具，不替代这些工具的 calling 算法。

当前公开候选版本为 **0.1.0rc4**，正式运行目标是 Linux/HPC，仅供科研使用，
不是临床诊断系统。

## 已包含的能力

- pbmm2、DeepVariant、Sawfish、Sniffles2、pbsv、cuteSV、TRGT、HiPhase、
  hifiasm、PAV、SVIM-asm、Jasmine、Truvari、GLnexus、ANNOVAR、VEP、
  IGV 和 hap.py 的安全 wrapper、边界与 provenance；
- 相互独立的 SNV/Indel、SV 和 TR 结果体系；
- 可选下游分支的模块化 Snakemake rules；
- QC、manifest、benchmark、review 和 final report 数据模型；
- wheel、sdist、本地 Conda recipe、Docker core 和 Apptainer core 定义。

Python 包和 core container **不包含**外部 caller、容器镜像、参考基因组、
数据库、cache 或 WGS 数据，启用相应分支前必须独立部署并验证。

## 当前执行边界

公开 Snakemake 的 small/SV/TR 分支以已经比对并建立索引的 BAM/CRAM 为输入。
Phase 2 Python API 可规划和运行 pbmm2 比对，但 `0.1.0rc4` 尚无统一的
`hifivar run` 命令，也没有打包的 Snakemake alignment rule，因此不能宣称
“FASTQ 一条命令完成全部变异分析”。原始 FASTQ/uBAM 必须先显式完成比对；
hifiasm assembly 分支则以 HiFi FASTQ 为输入。

## 系统要求

- 生产环境：Linux/HPC；
- Python 3.10–3.12；
- Snakemake 8 或 9；
- 已建立索引的参考 FASTA；
- 各启用分支所需的固定版本 executable/container/database；
- 与数据规模匹配的 CPU、内存、存储和调度器配置。

Windows 可用于 Python、配置、打包、mock 测试和 Snakemake dry-run，不作为
正式 external-tool 生产运行平台。

## 从 GitHub Release 安装

```bash
curl -fLO https://github.com/lg10is1/hifivar/releases/download/v0.1.0-rc4/hifivar-0.1.0rc4-py3-none-any.whl
curl -fLO https://github.com/lg10is1/hifivar/releases/download/v0.1.0-rc4/hifivar-0.1.0rc4.tar.gz
curl -fLO https://github.com/lg10is1/hifivar/releases/download/v0.1.0-rc4/SHA256SUMS
sha256sum -c SHA256SUMS

python3 -m venv hifivar-env
source hifivar-env/bin/activate
python -m pip install --upgrade pip
python -m pip install './hifivar-0.1.0rc4-py3-none-any.whl[workflow]'
```

当前尚未发布到 PyPI、Bioconda 或 conda-forge，请使用 GitHub Release 文件或
固定 tag 的源码安装。

## 验证安装

```bash
hifivar --version
hifivar --help
hifivar config validate
hifivar doctor
```

预期版本为 `hifivar 0.1.0rc4`。

## 五分钟 smoke test

```bash
WORKFLOW_ROOT="$(python -c 'from hifivar.package_resources import installed_workflow_root; print(installed_workflow_root())')"
hifivar --preset standard config dump-effective --output effective_config.yaml
snakemake \
  --snakefile "$WORKFLOW_ROOT/Snakefile" \
  --configfile effective_config.yaml \
  --cores 1 \
  --dry-run
```

这一步只验证核心安装与 DAG，不调用真实生物信息工具。

## 开始真实分析

1. 阅读[完整中文快速入门](docs/zh_CN/quickstart.md)；
2. 按[中文安装指南](docs/zh_CN/installation.md)完成核心安装；
3. 复制并修改[minimal 示例](examples/minimal/README.md)；
4. 按[部署矩阵](docs/deployment.md)部署需要的外部工具；
5. 先验证 reference、index、sample、contig 和 dry-run；
6. 单样本使用[中文 Bash 指南](docs/zh_CN/single_sample_bash.md)；
7. 多样本集群使用[中文 Slurm 指南](docs/zh_CN/slurm_multi_sample.md)。

HiFiVar 不会静默进行 `chr1`/`1` 转换，不会覆盖 raw caller VCF，也不会把
caller 支持数、人工 review、注释影响或 benchmark 状态等同于生物学真值或
致病性。

## 更多文档

- [完整中文快速入门](docs/zh_CN/quickstart.md)
- [中文安装指南](docs/zh_CN/installation.md)
- [单样本 Bash 运行](docs/zh_CN/single_sample_bash.md)
- [多样本 Slurm 运行](docs/zh_CN/slurm_multi_sample.md)
- [输出说明（英文）](docs/outputs.md)
- [常见问题（英文）](docs/troubleshooting.md)
- [外部工具部署矩阵（英文）](docs/deployment.md)

## License

HiFiVar 使用 [Apache License 2.0](LICENSE)。第三方工具、容器、参考数据和
数据库遵循各自许可证，HiFiVar 不重新分发它们。
