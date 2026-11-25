# 由主入口衔接的基准评估

[English](../benchmark_run.md)

显式计划绑定待评价产物、真值和可信区域。GENERATED 只指本次实际产生并验收的结果，EXTERNAL 需提供外部来源身份。hap.py、Truvari 与同目录 TR 比较保留独立任务和结果。查看两侧分母、工具退出状态以及同一次回执，失败尝试原样保留。

## 运行命令

下面的示例使用项目路径占位符。完整配置字段和参数以英文参考为准。

~~~bash
hifivar --config analysis.yaml --dry-run run --run-dir /your/project/runs/bench01
hifivar --config analysis.yaml run --run-dir /your/project/runs/bench01 --resume
~~~

[完整参数、配置和边界](../benchmark_run.md)
