# HiFiVar

[English](README.md) · [中文文档](docs/zh_CN/index.md)

HiFiVar 是面向 PacBio HiFi 的科研分析框架，可按配置衔接小变异、结构变异和串联重复分析，并保留来源和运行记录。**0.1.0rc5** 面向 Linux 和 HPC 环境。请按启用的分支另行准备调用器、参考与索引。[GitHub 发布页](https://github.com/lg10is1/hifivar/releases/tag/v0.1.0-rc5)

1. 下载安装固定版本。[安装](docs/zh_CN/installation.md)
2. 从小样本和一个已配置的调用器开始，检查计划后再执行。[快速入门](docs/zh_CN/quickstart.md)
3. 多样本参考[批量运行](docs/batch_multivariant_run.md)与[中文 Slurm 指南](docs/zh_CN/slurm_multi_sample.md)。

[RC5 主入口](docs/zh_CN/indexed_bam_run.md)使用已比对的索引 BAM。[可选输入准备](docs/raw_input_run.md)可将支持的 FASTQ 或 native uBAM 比对成 BAM。主入口暂不接入 CRAM。支持 CRAM 的独立 Snakemake 分支按自身条件运行。历史版 0.1.0rc4 仍可由旧 tag 安装。[功能与证据](docs/capabilities.md)

原生结果、派生结果与人工判定分别记录。更换工具和参考后先做小样本真实试跑。[输出](docs/zh_CN/outputs.md) · [基准](docs/benchmarking.zh-CN.md) · [已知局限](docs/zh_CN/known_limitations.md)

[Apache-2.0](LICENSE)。第三方工具和资源保留各自许可。
