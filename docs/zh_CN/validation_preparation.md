# 候选变异复核与PCR验证准备

[English](../validation_preparation.md)

本页所列复核准备入口随 RC5 提供。
功能用于准备可检查材料，人工判读、实验成功和GIAB准确率分别记录。

## IGV自动截图

`review prepare`仍只准备文件。新增`review snapshots`会执行IGV，要求明确的
`review.source_plan`外部来源、兼容的参考及样本BAM/索引。GENERATED来源继续通过
原有launcher工作流，不把预期输出当作已完成调用结果。

```bash
# 从本候选源码目录安装；IGV、Java和samtools另行配置。
python -m pip install '.[review]'
hifivar --config config.yaml --dry-run review snapshots \
  --settings capture_settings.yaml --output-dir /work/review_01
hifivar --config config.yaml review snapshots \
  --settings capture_settings.yaml --output-dir /work/review_01
```

窗口可通过YAML配置。SNV/小INDEL采用近景，长DEL/INV等采用整体加两端点；BND
按选中ALT显示可解析端点，不猜测未知位置。截图前执行samtools样本、参考、排序、
EOF与索引访问检查；截图后检查PNG完整性、解码、尺寸和数量。
直接截图入口还为每张图保存IGV session XML，核对实际参考、BAM/VCF资源及窗口，
防止把错误的全基因组空视图当作成功。IGV 2.19.8的genome命令保留字面引号；
普通参考路径使用未加引号的batch参数；原生2.19.8部署中的含空格参考尚未可靠通过，
此版在执行前明确拒绝，不声称已经支持。含空格的输出目录可采用安全相对session
文件名及明确工作目录；不自动移动或替换参考。

打开输出中的`index.html`批量查看。默认状态为`NOT_REVIEWED`。有效图片不证明
轨道载入、变异或基因型正确；其他IGV用户偏好尚未全部固定，需原生截图验收。
同一BAM上的人工观察也不能替代PCR或独立真值。

实际复核后填写`manual_decisions.tsv`中的判定、reviewer、带时区时间和说明，另存
TSV并导入独立目录。原图、VCF与BAM保持不变。

```bash
hifivar review import-decisions --bundle-dir /work/review_01 \
  --decisions reviewed.tsv --output-dir /work/review_decisions_01
```

完整配置和限制见[截图执行说明](../review_snapshots.md)。

## 引物候选设计

首版接收同一参考上的精确双等位SNV及序列明确、简单锚定DEL/INS。输入计划明确
VCF、样本、contig、POS、REF和ALT，检查该样本携带ALT、参考REF一致、contig长度
以及记录身份。复杂替换、多ALT、符号插入、INV/BND/DUP和不完整GT暂不进入设计。

```bash
# primer3_core需作为外部程序在PATH中，或在计划中填写明确路径。
hifivar --dry-run primers design --plan examples/primer_design/plan.yaml \
  --output-dir /work/primer_candidates_01
hifivar primers design --plan examples/primer_design/plan.yaml \
  --output-dir /work/primer_candidates_01
```

示例是人工构造的小序列，不是人类数据、真实PCR或准确率结果。用户需换成自己的
匹配参考与候选来源。所有输出使用新目录，失败后保留证据，不覆盖原VCF。

默认目标是三对位置不同的候选；产物范围、Tm、GC及引物长度可在计划`settings`
调整。合格不足时返回实际0/1/2对，不能凑数。Primer3通过统一CommandRunner运行，
保留版本、参数、日志、序列方向和REF/ALT产物重构。

输出：

- `primer_sequences.tsv` 提供可直接复制的正向和反向引物序列，均按 5′→3′ 输出，并列出引物对编号、预测 REF/ALT 长度及验证状态。CLI 完成时返回该文件路径。
- `primers.tsv` 提供完整候选详情，包括序列、坐标、Tm/GC、产物长度和计算指标。
- `predicted_amplicons.fasta`：REF与ALT预测扩增序列。
- `predicted_bands.tsv`：按原GT列出的可能产物长度类别；实际条带未测量。
- `primer_design_summary.json`及逐位点目录：来源、候选不足原因和原生调用记录。

**SNV默认是PCR后Sanger测序方案。** 普通PCR的REF/ALT产物通常等长，常规跑胶
不能确认碱基替换。DEL/INS即便长度不同，也未保证胶分辨率、扩增效率、等位脱落或
非特异产物。首版未执行全基因组特异性、邻近样本变异检查，状态保留`NOT_ASSESSED`；
实验保留`NOT_RUN`。这些候选不能直接称为“已验证可靠引物”。ARMS/RFLP及复杂SV
连接PCR需另行实现和验收。

Primer3默认热力学参数依赖部署；自定义参数目录路径会记录，但当前尚未封存整个
参数数据库。进入正式实验前应完成固定版本的特异性检查、参数归档和实验室审核。

详细方法见[引物核心说明](../primer_design.md)。
