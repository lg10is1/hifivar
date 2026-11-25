# SV 深度响应与支持数研究

[English](../sv_depth_study.md)

根据模板填写各技术单元的索引 BAM、实测覆盖、参考、真值与可信区域。冻结 caller 和比较参数。先验证计划，再运行所选单元。支持数研究单独使用对应的成功原生调用。覆盖标签与实测覆盖分别记录，同供者技术重复不增加供者数。

## 运行命令

下面的示例使用项目路径占位符。完整配置字段和参数以英文参考为准。

~~~bash
python -m pip install '.[sv-benchmark]'
hifivar benchmark sv-depth --help
~~~

~~~bash
hifivar --dry-run benchmark sv-depth --manifest sv_depth.yaml --output-dir /your/project/sv_depth/study01
hifivar benchmark sv-depth --manifest sv_depth.yaml --output-dir /your/project/sv_depth/study01 --unit-id 5x-r1 --caller sniffles2
~~~

~~~bash
hifivar benchmark sv-depth --manifest sv_depth.yaml --output-dir /your/project/sv_depth/study01 --resume
~~~

~~~bash
hifivar --dry-run benchmark sv-depth-support --manifest support.yaml --output-dir /your/project/sv_depth/support_attempt01
hifivar benchmark sv-depth-support --manifest support.yaml --output-dir /your/project/sv_depth/support_attempt01 --unit-id 5x-r1
~~~

[完整参数、配置和边界](../sv_depth_study.md)
