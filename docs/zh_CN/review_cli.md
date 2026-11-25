# IGV 复核包与人工判定导入

RC5还提供独立的 [`review tr-plot --plan` 入口](trgt_plot.md)，只绘制一个
已有 TRGT 5.1.0 位点。它不改 IGV 复核包，也不代替人工判读；准确候选包及原生
绘图验收与下文较早的离线 CLI 验收分开记录。

这是已有 Phase 10 的 **RC5 候选增量**，不在 RC1–RC4 或早期 RC5 包内。
它已包含于[指定的 10 月 15 日开发包](../rc5_candidate.md#validated-development-payload--已验收开发包)，
并通过离线 CLI 的受限 Linux、回归和包检查。收到新包后仍应核对源码清单和 wheel
SHA256，不能仅看版本号。[English](../review_cli.md)。

## 准备：一条命令生成待复核包

先准备已有比对文件、reference FASTA/FAI、样本表和明确选定的 variant TSV。
原始 FASTQ/uBAM 不能当作 IGV alignment。参考 build 必须明确声明。

```yaml
reference:
  fasta: /your/reference/GRCh38.fa
  build: GRCh38
samples:
  sheet: /your/analysis/samples.tsv
review:
  enabled: false
  selection_file: /your/analysis/selected.tsv
  igv_executable: /your/tools/igv.sh
  flank_bp: 500
  overwrite: false
```

`review.enabled` 控制独立 Snakemake 分支，不是这个显式准备命令的必需开关。
YAML 内 selection 的相对路径以 YAML 为基准；selection 内 VCF 路径以 TSV 为基准；
alignment 路径以样本表为基准。先规划，再使用同一个尚未创建的输出目录准备：

```bash
hifivar --config review.yaml --dry-run review prepare --output-dir /your/project/review/bundle01
hifivar --config review.yaml review prepare --output-dir /your/project/review/bundle01
```

dry-run 读取轻量输入并输出计划，不写文件。实际准备会写：

- `review.igv.batch`：IGV 批处理文件；
- `review_manifest.json/.yaml/.tsv`：全部保持 `NOT_REVIEWED` 的原始复核清单；
- `manual_decisions.tsv`：可填写的人工判定模板；
- `effective_config.yaml`、`review_bundle.json`：配置与最后完成凭据；
- 空 `screenshots/`：预留位置，不表示截图已经生成。

不会执行 IGV、caller、自动选高可信位点或修改 VCF。输出目录必须全新；拒绝覆盖。
若写入失败，保留现场，无完成凭据的目录不能作为成功 bundle；不会自动清理重试。

selection 的必需列是 `review_id/sample/variant_id/variant_type/contig/start/end/source_vcf/source_caller/evidence_class`。
必须使用实际制表符。可选列为 `mate_contig/mate_position/flank_bp/trgt_visualization`。
坐标为 1-based、两端包含；BND 需要 mate，并为两个断点分别规划截图。
位点坐标检查 FAI 边界，但沿用原 Phase 10 的 flanking window，右侧窗口可能超过 contig 尾端，须在 viewer 中核对。

## 人工复核：必须由真实 reviewer 完成

在已安装并具备显示环境的 IGV 中，选择 **Tools > Run Batch Script**，打开生成的
batch 文件，或使用准备结果给出的启动命令。官方说明见
[IGV batch guide](https://igv.org/doc/desktop/UserGuide/tools/batch/)。

查看两条等位基因、覆盖、比对质量、soft/supplementary alignment、局部复杂度及其他
证据后，填写 reviewer、时间和说明。没有 reviewer 真实判读，不能因截图存在或多个
caller 支持就标记完成。本命令也不会代替此前暂缓的 24 个 IGV 位点作判定。

## 导入：判定与原始 caller 结果分开存储

可复制模板到 `reviewer.tsv` 后填写。列必须精确为：

```text
review_id sample variant_id status reviewer reviewed_at comment
```

以上只展示列名，实际文件使用 tabs。状态为：

- `NOT_REVIEWED`：未复核，reviewer/时间/comment 留空；
- `SUPPORT`：观察到支持；
- `NOT_SUPPORT`：未观察到支持；
- `UNCERTAIN`：证据不明确。

这不是 truth、pathogenicity、ACMG 或 PCR 验证结果；不会隐式把 PASS/REJECT 转为另一种语义。
已经作判定的行必须提供 reviewer、非空 comment 和带时区的 ISO 时间。
例如 `2026-10-14T10:00:00+08:00` 只是格式示例，不是真实复核日期记录。

```bash
hifivar --dry-run review import-decisions --bundle-dir /your/project/review/bundle01 --decisions /your/project/review/reviewer.tsv --output-dir /your/project/review/adjudicated01
hifivar review import-decisions --bundle-dir /your/project/review/bundle01 --decisions /your/project/review/reviewer.tsv --output-dir /your/project/review/adjudicated01
```

review_id、sample、variant_id 必须与 prepared bundle 精确匹配。重复/未知 ID、非法
状态、缺少 reviewer/comment、无时区日期或被改动的 sealed 文件会失败。
允许部分复核，缺失行仍保留 NOT_REVIEWED；空表不代表全部复核完成。
结果写入独立的新目录，保留原声明顺序、原始 manifest 和所有调用结果。

## 必须保留的限制

1. 准备只复用 Phase 10 轻量校验，不证明 BAM header/sample、完整 records 或选中
   VCF allele/ID 的身份正确，也不能替代准确率 benchmark。
2. 大型输入记录的是大小/mtime，**不是全文件 SHA256**；selection 和生成元数据有
   hashes。metadata 不证明原始大文件逐字节不可变，hash 也不是签名或身份认证。
3. 判定导入不重新验证大型原始文件或截图，也不认证 reviewer 本人是否进行了复核。
   它记录用户提交的独立声明，不把导入成功冒充独立人工复核证据。
4. 目录含样本名、绝对路径和 comment，默认作为私有证据保存。配置敏感键脱敏不等于
   任意自由文本安全；不要在路径/说明中写入凭据，不要默认上传 bundle。
5. 离线 CLI 的 Linux 验收与实际 IGV 截图、人工判读分开；验收通过不代表 reviewer
   已完成复核。没有新增科研指标、calling 算法或 Phase 16。

## BND 双端来源校验

开发版BND复核绑定真实选中VCF记录及ALT两端，单端与虚拟端点明确保留，不补造截图。派生BND来源需要显式bnd_bundle_manifest列和完整包验证；不得只提供孤立派生VCF。详见[双端消费者](../bnd_junction_consumers.md)。准备清单不等于人工已判读。
