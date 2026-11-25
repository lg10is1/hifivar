# 故障排查

[English](../troubleshooting.md)

先确认实际解释器与 HiFiVar 版本。检查参考、样本表、索引、可执行程序和配置，根据第一条失败日志定位前置环节。在新配置下先预演并确认任务状态，按回执恢复。保留已有输入、成功产物和失败证据。

## 运行命令

下面的示例使用项目路径占位符。完整配置字段和参数以英文参考为准。

~~~bash
hifivar --config config.yaml config validate
~~~

[完整参数、配置和边界](../troubleshooting.md)
