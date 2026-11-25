# 按调用器支持数筛选

[English](../sv_support_selection.md)

设置最低 caller 来源数，范围不得超过实际来源数。筛选产物需保留原生整合记录到逐记录来源证据的映射。若要评价不同门槛，分别在同一参照和规则下比较 P、R、F1 及错误计数。

## 运行命令

下面的示例使用项目路径占位符。完整配置字段和参数以英文参考为准。

~~~bash
hifivar --config analysis.yaml --dry-run run --run-dir /your/project/runs/support2
hifivar --config analysis.yaml run --run-dir /your/project/runs/support2 --resume
~~~

[完整参数、配置和边界](../sv_support_selection.md)
