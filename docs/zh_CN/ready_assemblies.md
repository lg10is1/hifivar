# 复用已有单倍型组装

[English](../ready_assemblies.md)

为每个样本显式指定两个已有单倍型 FASTA、来源与参考。验证配对、索引及文件身份，在配置里选择需要的组装 SV 调用器并分别部署依赖。PAV 混合原生 VCF 与 SV 派生视图各自保存。已有组装保留外部来源。

## 运行命令

下面的示例使用项目路径占位符。完整配置字段和参数以英文参考为准。

~~~bash
hifivar --config analysis.yaml run --run-dir /your/project/runs/ready01
hifivar report --run-dir /your/project/runs/ready01 --output-dir /your/project/reports/ready01
~~~

[完整参数、配置和边界](../ready_assemblies.md)
