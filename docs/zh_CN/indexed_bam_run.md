# 用索引 BAM 运行 RC5

[English](../indexed_bam_run.md)

先填写样本表和配置。要求非空的坐标排序 BAM 与索引，表中样本 ID 对应 RG/SM，参考 contig 和长度兼容。选择所需分支。检查实际计划后，以同一目录和配置续跑。支持范围外的 CRAM 应使用相应独立路线。对显式下游来源，预演可返回 GENERATED_SOURCE_PLAN_ONLY；待产生实际产物后再构建消费计划。

## 运行命令

下面的示例使用项目路径占位符。完整配置字段和参数以英文参考为准。

~~~bash
# Plan: performs real read-only BAM preflight and writes launch metadata/logs.
hifivar --config analysis.yaml --dry-run run --run-dir /your/project/runs/batch01

# Execute the same frozen configuration and input bytes after reviewing the plan.
hifivar --config analysis.yaml run --run-dir /your/project/runs/batch01 --resume

# Only after diagnosing a failed job and checking no scheduler/children remain:
hifivar --config analysis.yaml run --run-dir /your/project/runs/batch01 --resume --rerun-incomplete
~~~

[完整参数、配置和边界](../indexed_bam_run.md)
