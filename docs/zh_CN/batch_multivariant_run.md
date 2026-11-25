# 一次配置的批量多类型分析

[English](../batch_multivariant_run.md)

样本表绑定每个索引 BAM 与自己的样本 ID。配置指定小变异、SV、TR 和可选定相；为消费阶段分别声明 GENERATED 或 EXTERNAL 来源。先 dry-run，确认要启动的 producer 与后续任务。执行同一配置后记录真实 attempt ID，再使用独立 report 命令。示例中的路径、镜像、工具及目录均需替换。

## 运行命令

下面的示例使用项目路径占位符。完整配置字段和参数以英文参考为准。

~~~bash
hifivar --config config.yaml config validate
hifivar --config config.yaml doctor --preflight --run-dir /your/project/runs/batch01 --json
hifivar --config config.yaml doctor --preflight --run-dir /your/project/runs/batch01 --probe-tools --timeout-seconds 30 --json
~~~

~~~bash
# Reads real inputs and writes plan metadata; no calling.
hifivar --config config.yaml --dry-run run --run-dir /your/project/runs/batch01

# After reviewing the exact plan, execute the same frozen controls.
hifivar --config config.yaml run --run-dir /your/project/runs/batch01 --resume

# Use the ACTUAL immediate-child attempt ID recorded by that invocation.
hifivar report --run-dir /your/project/runs/batch01 --attempt-id attempt_ACTUAL_ID --output-dir /your/project/reports/batch01
~~~

[完整参数、配置和边界](../batch_multivariant_run.md)
