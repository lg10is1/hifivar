# HiFiVar RC5 中文文档

[English](../index.md)

Start with installation and quick start. / 首次使用先看安装与快速入门。

| English | 简体中文 |
| --- | --- |
| [Install HiFiVar 0.1.0rc5](../installation.md) | [HiFiVar 0.1.0rc5 安装](installation.md) |
| [HiFiVar 0.1.0rc5 quick start](../quickstart.md) | [HiFiVar 0.1.0rc5 中文快速入门](quickstart.md) |
| [Development](../DEVELOPMENT.md) | [开发与测试](DEVELOPMENT.md) |
| [AlphaGenome: explicit user-key functional prediction](../alphagenome.md) | [AlphaGenome 可选功能预测](alphagenome.md) |
| [One-invocation batch small/SV/TR analysis](../batch_multivariant_run.md) | [一次配置的批量多类型分析](batch_multivariant_run.md) |
| [Owned benchmark execution](../benchmark_run.md) | [由主入口衔接的基准评估](benchmark_run.md) |
| [Benchmarking](../benchmarking.md) | [基准测试](../benchmarking.zh-CN.md) |
| [BND junction-aware annotation and review](../bnd_junction_consumers.md) | [双端断点的注释与复核](bnd_junction_consumers.md) |
| [Explicit BND representation bundles](../bnd_representation.md) | [显式 BND 表示](bnd_representation.md) |
| [Capabilities, execution and evidence](../capabilities.md) | [功能、执行和证据](capabilities.md) |
| [Conda, Docker and Apptainer](../containers.md) | [Conda、Docker 和 Apptainer](containers.md) |
| [Validated deployment matrix](../deployment.md) | [依赖部署与真实工具边界](deployment.md) |
| [Public distribution validation](../distribution_validation.md) | [发行包验证与历史记录](distribution_validation.md) |
| [IGV screenshot integrity and review status](../igv_snapshot_quality.md) | [IGV 截图完整性与复核状态](igv_snapshot_quality.md) |
| [Indexed BAM batch launcher](../indexed_bam_run.md) | [用索引 BAM 运行 RC5](indexed_bam_run.md) |
| [Joint small/SV/TR phasing](../joint_phasing.md) | [联合定相](joint_phasing.md) |
| [Known Limitations](../known_limitations.md) | [已知局限性](known_limitations.md) |
| [Linux/HPC execution](../linux_hpc.md) | [Linux/HPC 运行](linux_hpc.md) |
| [Output path safety](../output_path_safety.md) | [输出路径保护](output_path_safety.md) |
| [Output layout](../outputs.md) | [输出布局与完成状态](outputs.md) |
| [Calling preflight: check before a fresh batch run](../preflight.md) | [批量运行前检查：先知道缺什么，再开始计算](preflight.md) |
| [Candidate PCR](../primer_design.md) | [候选位点的本地引物设计](primer_design.md) |
| [HiFiVar Project Status](../project_status.md) | [项目状态与适用范围](project_status.md) |
| [FASTQ and native uBAM batch execution](../raw_input_run.md) | [由 FASTQ 或 native uBAM 开始](raw_input_run.md) |
| [HiFiVar 0.1.0rc5 release identity](../rc5_candidate.md) | [RC5 发行身份](rc5_candidate.md) |
| [Existing haplotype FASTAs: explicit reuse contract](../ready_assemblies.md) | [复用已有单倍型组装](ready_assemblies.md) |
| [Reproducibility](../reproducibility.md) | [可重复性](reproducibility.md) |
| [Offline review preparation and decision import](../review_cli.md) | [IGV 复核包与人工判定导入](review_cli.md) |
| [Execute candidate IGV snapshots](../review_snapshots.md) | [执行候选 IGV 截图](review_snapshots.md) |
| [Explicit source review — development execution bridge](../review_sources.md) | [显式复核来源](review_sources.md) |
| [One-sample Linux Bash workflow](../single_sample_bash.md) | [HiFiVar 0.1.0rc5 单样本 Bash 运行指南](single_sample_bash.md) |
| [Multi-sample execution with Slurm](../slurm_multi_sample.md) | [HiFiVar 0.1.0rc5 多样本 Slurm 运行指南](slurm_multi_sample.md) |
| [SV depth-response study](../sv_depth_study.md) | [SV 深度响应与支持数研究](sv_depth_study.md) |
| [Within-sample SV integration](../sv_harmonization_run.md) | [同一供者的 SV 整合](sv_harmonization_run.md) |
| [Explicit SV source identity](../sv_source_identity.md) | [显式 SV 来源身份](sv_source_identity.md) |
| [Caller-support selection](../sv_support_selection.md) | [按调用器支持数筛选](sv_support_selection.md) |
| [Plot one existing TRGT locus](../trgt_plot.md) | [绘制一个已有 TRGT 位点](trgt_plot.md) |
| [Troubleshooting](../troubleshooting.md) | [故障排查](troubleshooting.md) |
| [Candidate review and validation preparation](../validation_preparation.md) | [候选变异复核与PCR验证准备](validation_preparation.md) |
| [VEP output safety](../vep_output_safety.md) | [VEP 输出保护](vep_output_safety.md) |

## Development phases / 开发阶段

- [00 English](../phases/phase00.md) · [基础设施](phases/phase00.md)
- [01 English](../phases/phase01.md) · [生物信息数据模型](phases/phase01.md)
- [02 English](../phases/phase02.md) · [质控与比对](phases/phase02.md)
- [03 English](../phases/phase03.md) · [DeepVariant 小变异](phases/phase03.md)
- [04 English](../phases/phase04.md) · [读段 SV](phases/phase04.md)
- [05 English](../phases/phase05.md) · [TRGT 串联重复](phases/phase05.md)
- [06 English](../phases/phase06.md) · [HiPhase 定相](phases/phase06.md)
- [07 English](../phases/phase07.md) · [单样本组装](phases/phase07.md)
- [08 English](../phases/phase08.md) · [组装型 SV](phases/phase08.md)
- [09 English](../phases/phase09.md) · [SV 整合](phases/phase09.md)
- [10 English](../phases/phase10.md) · [IGV 与人工复核](phases/phase10.md)
- [11 English](../phases/phase11.md) · [注释与候选评分](phases/phase11.md)
- [12 English](../phases/phase12.md) · [队列分析](phases/phase12.md)
- [13 English](../phases/phase13.md) · [真值基准](phases/phase13.md)
- [14 English](../phases/phase14.md) · [报告与打包](phases/phase14.md)
- [15 English](../phases/phase15.md) · [低深度研究](phases/phase15.md)
