# 执行候选 IGV 截图

[English](../review_snapshots.md)

提交候选选择及已核对的 VCF、BAM、索引和参考，配置 IGV、Java 与所需显示环境。先预演路径与候选范围，再生成图像和会话。检查会话和产物状态。人工复核需导入明确的决策表。

## 运行命令

下面的示例使用项目路径占位符。完整配置字段和参数以英文参考为准。

~~~bash
hifivar --config config.yaml --dry-run review snapshots \
  --settings capture_settings.yaml --output-dir /work/review_candidate_01
hifivar --config config.yaml review snapshots \
  --settings capture_settings.yaml --output-dir /work/review_candidate_01
~~~

~~~bash
hifivar review import-decisions --bundle-dir /work/review_candidate_01 \
  --decisions reviewed.tsv --output-dir /work/review_decisions_01
~~~

[完整参数、配置和边界](../review_snapshots.md)
