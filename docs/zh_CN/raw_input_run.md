# 由 FASTQ 或 native uBAM 开始

[English](../raw_input_run.md)

在配置中显式启用 alignment 和匹配的参考、后端、线程及临时目录。检查样本表中的输入类型。dry-run 只验证输入路线并返回 RAW_INPUT_PLAN_ONLY。实际过程先生成索引 BAM，随后检查并运行调用 DAG。要验证索引 BAM 路线可从已生成比对单独试跑。

## 运行命令

下面的示例使用项目路径占位符。完整配置字段和参数以英文参考为准。

~~~bash
hifivar --config analysis.yaml --dry-run run --run-dir /your/project/runs/raw_batch01
hifivar --config analysis.yaml run --run-dir /your/project/runs/raw_batch01 --resume
~~~

[完整参数、配置和边界](../raw_input_run.md)
