# 绘制一个已有 TRGT 位点

[English](../trgt_plot.md)

这是 RC5 中独立的绘图入口，不在 RC1–RC4 或更早的同版本 RC5 包内。
请核对准确的源码清单、包 SHA256 和验收范围，不能仅凭版本号判断功能存在。
它不运行 TR genotype，不改 calling DAG，不执行 IGV，也不替代人工复核。

## 输入与支持范围

准备**同一次原始 TRGT 分析**的 VCF、TBI、spanning-read BAM、BAI，以及相应
reference FASTA/FAI 和 repeat catalog。当前只支持原生 **TRGT 5.1.0**，版本探测必须
精确报告 `5.1.0` 或已观测部署构建 `5.1.0-ec66463`，以及一个明确的
`TRID:<catalog ID>`、allele/motifs 图及 SVG。未开放 methylation、waterfall、
deepdive、任意 native 参数或 threads 参数。

计划必须精确绑定七个文件的 SHA256：`fasta/fai/catalog/vcf/tbi/spanning_bam/bai`。
哈希证明输入字节与声明一致，不证明 producer 身份，也不能单凭哈希确认 spanning
BAM 确实由所给 VCF 的同一次分析产生。初版来源模式仅支持 `EXTERNAL_FILES`。
原始 producer 记录应另外保留，不能用手填的 attempt 字符串代替。

程序检查单个声明样本、唯一 TRID、catalog/VCF 坐标、MOTIFS、STRUC、完整 GT
和明确的序列等位基因；从实际 FASTA 随机读取包含前置 anchor 的完整 REF 并核对。
FAI 必须是 `FASTA.fai`；VCF/TBI、BAM/BAI 使用显式常规索引名。reference 与 catalog
build 必须一致，不静默改染色体名、等位表示或基因型。这些检查不等于完整 BAM
来源审计或匹配真值准确率 benchmark。

## 准备严格 YAML 计划

复制[通用示例](../../examples/review/trgt_plot.yaml)到自己的项目目录，替换所有路径、
样本、位点和七个 SHA 占位符；未编辑模板不能直接运行。例如对实际文件计算：

```bash
sha256sum reference.fa reference.fa.fai repeats.bed sample.tr.vcf.gz sample.tr.vcf.gz.tbi sample.tr.spanning.bam sample.tr.spanning.bam.bai
```

所有以下字段都必须明确提供：

- `schema_version: 1`、`sample_id`、`selector: TRID:<原始 catalog ID>`；
- `reference: {fasta, fai, build}`、`catalog: {path, build}`；
- `source: {mode, vcf, tbi, spanning_bam, bai}`；
- `expected_sha256`：精确七个角色、64 位小写十六进制 SHA；
- `runtime: {executable, supported_version}`：版本字符串为 `"5.1.0"`；
- `plot: {type, show, format, max_allele_reads}`：allele/motifs/svg，加正整数显示上限；
- `execution: {timeout_seconds}`：有限正数，不超过 1200 秒；
- `output: {directory}`：一个尚不存在的新目录，父目录必须已存在。

相对路径以计划文件所在目录为基准。未知、缺失、重复、合并、别名或非字符串 YAML
key 都拒绝。`runtime.executable` 是一个可执行文件名或路径，不是拼接 shell 参数的
字符串。loader 会单独绑定原始计划字节 SHA，不把它与规范化后的 request 混淆。
timeout 对版本探测和 plot 命令分别生效，不是整个 consumer 的总墙钟上限；
输入验证和重复哈希也需要时间。

## 先预览，再执行

```bash
hifivar --dry-run review tr-plot --plan trgt_plot.yaml
hifivar review tr-plot --plan trgt_plot.yaml
```

dry-run 完整验证文件、SHA 和选中位点，但不探测版本、不运行 TRGT、不创建输出
目录。即使只画一个位点，流式核对整份 reference 的 SHA 仍会产生读取时间。
实际绘图前需另行部署 TRGT 5.1.0，并正确填写其 executable。

执行通过统一 CommandRunner，在独占新目录内记录版本及 plot 各自的 stdout/stderr、
有效输入、实际 argv、工具身份和文件 hashes。版本日志保留原始字节；receipt 中
`tool_version` 规范记录为 `5.1.0`，`tool_version_reported` 完整保留实际观测 token。
未知后缀、版本行中的额外尾注及 stdout/stderr 中相互冲突的版本报告会被拒绝，
不静默去除后缀放行。
计划的 `runtime.supported_version` 仍填 `"5.1.0"`。支持这个明确部署 token 不等于
新包已通过原生绘图验收，新包仍需自己的执行证据。主图为 `plot.svg`。正常成功 CLI 状态
为 `COMPLETE_TECHNICAL_SVG`，`PLANNED` 只代表预览；两者均不代表准确率或人工
复核完成。

## 完成与失败的判断

固定证据包括 `effective_input.json`、`plot.svg`、各自的 `version.*.log` 和
`command.*.log`、保留的 `.completion.pending.json`，最后发布 `completion.json`。
程序再次核对来源、工具、计划与输出身份；异常还可能留下 `failure.json`。

必须同时满足**进程 exit 0、completion 与实际文件匹配、不存在 `failure.json`**，
才能接受本次技术结果。只有文件存在不算成功；若恰在 receipt link 后中断，可能留下
receipt 但没有正常返回。不能据此默认为成功或在该目录自动续跑。

不提供 overwrite/resume，不自动重试，不删除或递归清理。已有或晚到的目标会明确
失败。失败目录保留，确认原因后才自行选择另一个新目录。该目录不是抵御恶意并发
改写的操作系统沙箱，receipt 发布也不是多文件物理原子事务或断电恢复保证。

SVG 使用保守静态白名单检查，拒绝脚本、动态结构和外部资源。但结构检查不是像素
QA：用于论文前仍需在可信 viewer 中查看字体、读段和布局。排版衍生图与原生 SVG
分开保存，不改写原始证据。

## 与人工 review 分开

receipt 固定保留 `manual_review: NOT_REVIEWED`、`accuracy: NOT_ESTABLISHED`。
原生等位 panel 的先后不代表父源/母源，也不证明 phasing 准确率；图中长度、motif
不是临床判定。

如需纳入后续离线 review 包，按[复核指南](review_cli.md)在 selection 的
`trgt_visualization` 列显式填写 SVG 路径。这不会执行 IGV、改旧复核包或替用户完成
人工判读。图和有效计划可能含样本名、绝对路径，未获得共享许可前作为私有证据保存。
