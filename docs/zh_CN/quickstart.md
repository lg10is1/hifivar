# HiFiVar 0.1.0rc5 中文快速入门

[English](../quickstart.md)

1. 按[安装说明](installation.md)确认版本。
2. 复制[配置](../../examples/minimal/config.yaml)和[样本表](../../examples/minimal/samples.tsv)。替换 `/data/...` 和 `/work/...` 路径。准备 FASTA 和 FAI、坐标排序的 BAM 和索引。样本 ID 与 BAM 的 RG/SM 一致，参考、BAM 和目录的染色体名称及长度应匹配。
3. 示例默认关闭调用器。配置一个已部署的调用器、镜像或可执行程序和对应资源。[部署](deployment.md)
4. 先用前两行预演并检查计划。确认输入和分析分支后，在计算节点运行最后一行，沿用同一配置、输入和目录。

~~~bash
hifivar --config config.yaml config validate
hifivar --config config.yaml --dry-run run --run-dir ./runs/sample01
hifivar --config config.yaml run --run-dir ./runs/sample01 --resume
~~~

主入口使用索引 BAM。[可选原始输入准备](raw_input_run.md)支持指定 FASTQ 和 native uBAM，预演返回 `RAW_INPUT_PLAN_ONLY`。比对完成后检查完整调用计划。主入口暂不接入 CRAM。

运行后检查各分支的原生 VCF、索引和回执。[输出](outputs.md) · [单样本 Bash](single_sample_bash.md) · [批量](batch_multivariant_run.md) · [Slurm](slurm_multi_sample.md)
