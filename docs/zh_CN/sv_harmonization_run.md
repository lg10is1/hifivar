# 同一供者的 SV 整合

[English](../sv_harmonization_run.md)

配置同一样本的多个 SV 来源，并明确生成或外部来源及顺序。Jasmine 整合前检查参考和 VCF 索引，整合后保留逐条来源向量和证据。配置比较时另看 Truvari 的匹配范围。跨样本聚类属于不同任务。

## 运行命令

下面的示例使用项目路径占位符。完整配置字段和参数以英文参考为准。

~~~bash
hifivar --config analysis.yaml doctor --preflight --run-dir /your/project/runs/sv_batch01 --probe-tools --json
hifivar --config analysis.yaml --dry-run run --run-dir /your/project/runs/sv_batch01
hifivar --config analysis.yaml run --run-dir /your/project/runs/sv_batch01 --resume
~~~

[完整参数、配置和边界](../sv_harmonization_run.md)
