# Batch small/SV/TR example / 批量多类型变异模板

[简体中文](README.zh-CN.md)

For the **RC5 release source** with joint-phasing and generated
consumer support. Not a runnable example for RC1-RC4 or arbitrary same-version
RC5 packages. See the complete bilingual [guide](../../docs/batch_multivariant_run.md).

Copy config.yaml, samples.tsv and source_plan.yaml into a new analysis directory
inside your chosen project root. Replace all /work paths and sample metadata;
provision compatible pinned tools, reference/index and TR catalog. BAM @RG SM
must exactly equal the sheet and plan IDs. Add/remove all three plan rows whenever
you add/remove a sample. Generated TR rows bind the configured catalog's actual
SHA256; no invented hash or filename identity is needed.

Scope: existing indexed BAMs -> small + four read-SV + TR -> within-sample
harmonization and small/SV/TR phasing -> GLnexus + separate SV/TR summaries.
Final report is a separate command. No assembly, annotation, benchmark or review
is enabled by this example. SV/TR summaries are not cross-sample clustering or
joint TR genotyping. Resource values are scheduler examples, not guarantees.

中文：本模板对应 GitHub RC5 发行包，不适用于旧RC1-RC4，也不能只凭0.1.0rc5版本号判断
可用性。先复制三文件到项目内的分析目录，再替换路径、真实样本ID/sex及工具入口。
样本表使用真实制表符；BAM中所有RG的SM必须与样本表/来源计划精确一致。
每增删一个样本，同步增删来源计划中的small_variants、sv、tr三行。
GENERATED TR自动绑定实际catalog内容SHA256，不手填虚假hash。

模板不安装工具、不自动创建原始BAM索引、不重命名样本/染色体、不默认覆盖或重试。
资源需按节点和数据检查。先用独立小型数据验证，不能把区域结果当全基因组结果。
