# Indexed BAM batch launcher / 已比对 BAM 批量入口

[简体中文](zh_CN/indexed_bam_run.md)

## Availability / 可用范围

Use the [capability/evidence matrix](capabilities.md) for the current staged
launcher. This guide retains historical payload boundaries explicitly.

For the RC5 opt-in benchmark consumer bridge, see
[owned benchmark execution](benchmark_run.md). Historical payload limits below
remain historical; the same version string does not identify all later source.

RC5 可显式衔接 benchmark，见上述专门说明；以下历史安装包的范围不变。

This is an RC5 candidate feature, included in the
[identified 15 October RC5 payload](rc5_candidate.md#validated-development-payload--已验收开发包)
but **not in the initial security-only RC5 package or RC1–RC4 assets**. Install
an audited candidate with the `workflow` extra before using this command.
Linux is required for execution; Windows remains a configuration/unit-test host.

这是 RC5 候选功能，已包含于上面指定的 10 月 15 日开发包，**不在最初的安全修复型
RC5 包或 RC1–RC4 发布包内**。请核对交付清单与包 SHA256，安装 `workflow` 依赖。
正式运行平台为 Linux；同版本的不同开发包不能混同。

The first launcher accepts one coordinate-sorted, indexed BAM per sample. It
routes existing DeepVariant, read-SV caller, TRGT and optional HiPhase rules.
With alignment disabled it does not align FASTQ/uBAM; the new opt-in
[raw-input preparation](raw_input_run.md) provides that separate first stage.
It does not accept CRAM. Later development revisions add explicit assembly
FASTQ contexts and ready-input annotation/cohort routing; these are not features
of the historical payload identified above. Current benchmark execution uses
the optional consumer bridge linked above; final reporting remains separate.
RC5 additionally routes opt-in
[within-sample read-SV harmonization](sv_harmonization_run.md). Other modules retain their separate documented
workflow entry points. Small variants, SV and TR remain separate outputs.

首版要求每个样本已有一个坐标排序、带索引的 BAM，可选择 SNV/Indel、read-SV、TRGT
及 HiPhase 分支。FASTQ/uBAM 可另行启用[原始输入准备阶段](raw_input_run.md)；CRAM、组装、
队列分析和注释并非首版功能；后续开发修订支持独立组装 FASTQ 上下文及已有输入的
注释/队列入口；最新benchmark执行桥可按显式来源计划衔接，最终报告仍单独生成。
最新开发增量可显式启用
[同样本 read-SV 整合](sv_harmonization_run.md)，不支持跨样本 SV 聚类。
不同变异类别保持独立输出。

## Inputs and configuration / 输入和配置

Create `samples.tsv` using actual tab characters. Sample IDs must exactly match
**every** BAM read-group `SM` field; IDs cannot be invented to label anonymous BAMs.
Multiple read groups are allowed if their IDs are unique and all SM values agree.
The example paths and sample IDs below must be replaced with your own.

```text
sample_id	input	sex
SAMPLE_A	/your/data/SAMPLE_A.bam	female
SAMPLE_B	/your/data/SAMPLE_B.bam	male
```

样本表必须用真正的制表符分隔。`sample_id` 必须与对应 BAM 中每个 `@RG` 的 `SM`
精确相等。索引必须是唯一且非空的常规 BAI/CSI；BAM 参考序列名称及长度必须与 FASTA
索引完全相同。不会自动重命名样本、转换 `chr1/1`、创建索引或改写输入。

Minimal DeepVariant configuration (`analysis.yaml`):

```yaml
reference:
  fasta: /your/reference/GRCh38.fa
  build: GRCh38
samples:
  sheet: samples.tsv
runtime:
  threads: 8
workflow:
  memory_mb: 64000
  samtools_executable: /your/envs/samtools/bin/samtools
small:
  enabled: true
  execution_mode: apptainer
  deepvariant_image: /your/containers/deepvariant_1.10.0.sif
  threads: 4
  memory_mb: 64000
  max_concurrent_samples: 1
  overwrite: false
```

Relative reference, sample-sheet, catalog, container and executable paths resolve
beside the user YAML; BAM paths resolve beside the original TSV. Bare executable
names use PATH. Images must already exist locally; the launcher does not install
tools or download resources. Configure all enabled tool executables/dependencies
in their existing YAML sections and environment before launching.

按需启用 `sv.enabled` 和相应 caller 的 `enabled`；TR 需要明确 catalog、
`catalog_reference_build` 和 karyotype/样本 sex；HiPhase 需要同时启用 small 分支。
各工具仍通过既有 wrapper 调用，科学参数不由启动器自动优化。

`runtime.threads` is the total scheduler core budget. Each enabled rule's threads
and memory must fit within it and `workflow.memory_mb`, respectively. No silent
thread downscaling is allowed. These are **scheduler reservations, not OS memory
or kernel-thread limits**. DeepVariant framework thread pools and HPC UID/cgroup
limits still require host-specific checks. The launcher makes no WGS runtime or
resource-sufficiency guarantee.

资源预算用于调度，不能替代节点内存、用户线程上限和容器资源检查。示例的 64 GB 不构成
所有数据规模都足够的承诺。各分析参数来自 YAML，不经命令行任意透传。

## Launch and resume / 启动与续跑

### Optional bounded small-variant calling / 可选小区域调用

In RC5, `small.region: null` (default) leaves the
DeepVariant command unchanged. To request one bounded interval, set for example:

```yaml
small:
  region: chr21:43770001-43785000
```

This is **1-based, inclusive**, not BED coordinates. Only one explicit
`contig:start-end` string is supported, with exact reference contig naming and
FAI length checks; BED paths, multiple intervals and automatic `chr` renaming
are rejected. The same interval applies to every sample and is frozen in the
launch configuration/identity. Changing it requires a new run directory.
The Phase 3 request/report records `calling_region` and
`requested_calling_scope`; the actual command contains `--regions` only when set.

当前开发源码支持以上可选配置；默认 `null` 保留原来的全参考调用行为。
区间为 **1-based、两端包含**，必须与配置的参考基因组匹配，不会自动改名或
推断基因组版本。改变区域必须建立新的运行目录；旧发布包不包含该新增能力。

This option limits **DeepVariant only**. It does not filter the input BAM,
restrict SV callers, replace the TR catalog, or change phasing parameters.
For a bounded joint acceptance test, prepare a separate indexed regional BAM
and an explicitly pinned regional TR catalog without changing original inputs.
The default small-only HiPhase route phases the supplied regional small VCF.
The later explicit [joint route](joint_phasing.md) can also consume selected
SV/TR sources; regional execution still does not establish whole-genome scope
or phasing accuracy. See the [complete batch template](batch_multivariant_run.md).
VCF/gVCF names stay unchanged;
retain the launch manifest/effective YAML with them so regional outputs are
never misrepresented as whole-genome calls or whole-genome cohort inputs.
Requested scope is not a coverage or accuracy guarantee; allele representations
and interval-edge effects remain native DeepVariant behavior. No VCF is clipped
or rewritten by HiFiVar. Phase 15 rejects this option rather than silently
ignoring it. Benchmark scope must be independently matched to truth/BED.

该参数只限制 DeepVariant，不会限制 SV、替换 TR catalog 或修改原始 BAM。
区域验收不等于全基因组准确率验证；输出必须连同配置及 manifest 保存，不能
把区域 gVCF 当成全基因组结果用于队列分析。Phase 15 暂不支持此选项。

Choose a **new** run directory inside your authorized project root. Leave
`paths.workdir`, `paths.outdir`, `runtime.tmpdir`, `logging.file` null; the launcher
assigns `work`, `results`, `tmp`, `logs/hifivar.log` beneath that directory.

```bash
# Plan: performs real read-only BAM preflight and writes launch metadata/logs.
hifivar --config analysis.yaml --dry-run run --run-dir /your/project/runs/batch01

# Execute the same frozen configuration and input bytes after reviewing the plan.
hifivar --config analysis.yaml run --run-dir /your/project/runs/batch01 --resume

# Only after diagnosing a failed job and checking no scheduler/children remain:
hifivar --config analysis.yaml run --run-dir /your/project/runs/batch01 --resume --rerun-incomplete
```

Global options (`--config`, `--dry-run`) precede `run`. A dry-run is **not** zero-I/O:
it records metadata, computes streaming SHA256 of inputs and probes BAM headers,
EOF and one indexed region. It does not call variant callers. This preflight is
not a full record-level BAM correctness or scientific-accuracy audit.

先 dry-run 查看完整计划，再用相同配置和目录加 `--resume` 执行。预检会完整流式校验
输入 SHA256，较大的 BAM/reference/container 会增加读取时间；不加载整个大文件到内存。
实际执行后再次核对输入内容。保持输入只读，避免其他进程同时修改输入。

Resume requires matching configuration, input bytes, interpreter and Python/
workflow payload hashes. Changed inputs/config/source require a new run directory.
No overwrite or automatic retry is enabled. An interrupted process retains a lock:
**do not remove it until the scheduler and descendants are independently confirmed
stopped**. Snakemake may need manual failure-specific remediation; this command
does not kill jobs or remove partial work directories.

续跑不是自动修复。出现失败时保留所有日志和中间产物，不自动删除、重启或重试。
不要把“文件存在”当成完成，也不要在原配置中改动参数后强行沿用运行目录。

The command runs in the foreground. Use your site's persistent session or batch
scheduler for long jobs. This version uses the local Snakemake executor; it does
not submit distributed Slurm jobs or promise survival after a broken SSH session.

## Outputs and evidence / 输出和证据

- `launch_manifest.json`: input checksums, runtime payload identity and configuration.
- `effective.yaml` and `samples.tsv`: derived, absolute-path launch configuration.
- `attempts/attempt_<id>/`: unique preflight/process logs and `status.json` per invocation.
- `results/`: existing per-module VCF/gVCF/TR outputs; no forced mixed-class VCF.
- `work/`, `tmp/`, `logs/`, `cache/`, `config/`: task-local working facilities.

`DRY_RUN_PASS` means plan acceptance. `WORKFLOW_SUCCEEDED` means the requested DAG
exited successfully and input bytes remained unchanged; module wrappers retain
their own output checks. Neither status establishes benchmark accuracy, clinical
validity, complete genome coverage, or success of modules not enabled in the run.

运行日志和 manifest 含本地样本名及绝对路径，默认属于私有分析证据，不应直接上传公开
仓库。发布性能结论仍需匹配 reference/truth 的独立 benchmark 和明确的数据范围。

## Explicit downstream handoff / 显式下游衔接

The current **RC5 interface** adds
`workflow.downstream_source_plan`. Its source tests are separate from installed
Linux acceptance; do not assume an older wheel with the same `0.1.0rc5` version
contains it. Use the delivery manifest and measured package hash.

本节为 RC5 的显式接口，不能用旧候选包的验收替代本增量的 Linux 安装验收。
一个命令可以先运行所选 caller，再将已通过本次执行回执校验的产物交给注释或队列
分析。它不是把未来文件提前标成成功，也不会将不同变异类别混成一个 VCF。

Add the following to the analysis YAML after configuring the ordinary caller
and downstream tools. Keep legacy `annotation.input_manifest` and
`cohort.input_manifest` null when using a source plan.

```yaml
workflow:
  downstream_source_plan: downstream_sources.yaml
cohort:
  enabled: true
  cohort_id: EXAMPLE_COHORT
  input_manifest: null
  small_variants:
    enabled: true
    threads: 8
    memory_gb: 8
```

The 8 GiB setting is an explicit **bounded-region example**, not a WGS memory
recommendation. Choose and validate resources for the actual cohort. Existing
DeepVariant/tool configuration, sample identity, build and index requirements
still apply; this fragment does not install tools or constrain genomic scope.

`downstream_sources.yaml` explicitly chooses the gVCF producer for every sample:

```yaml
schema_version: 1
cohort:
  - sample_id: SAMPLE_A
    track: small_variants
    origin: GENERATED
    producer: small.deepvariant
    output_role: gvcf
    reference_build: GRCh38
  - sample_id: SAMPLE_B
    track: small_variants
    origin: GENERATED
    producer: small.deepvariant
    output_role: gvcf
    reference_build: GRCh38
```

将样本名替换为真实样本表中的名称，并显式列出每个启用 cohort track 的全部样本。
GLnexus 只接受真正的 gVCF；普通 VCF 不可替代。SV/TR 队列保留分类型表和缺失状态，
不执行跨样本 SV 聚类。TR 来源还须通过同一 catalog 内容 SHA256 的核对。

Annotation rows use `track: small`, `sv` or `tr`, with an explicit producer and
`output_role: vcf`. Enable and configure ANNOVAR/VEP/overlap separately. ANNOVAR
and VEP retain separate results and annotate the whole selected source VCF;
candidate `variant_ids` are not a hidden record filter for either backend.
Region overlap requires explicit selected IDs for SV/TR. External rows instead
use `origin: EXTERNAL`, existing VCF/index paths and explicit source metadata;
they are validated, not rewritten or silently normalized.

For TRGT records whose raw VCF ID is `.`, choose the INFO/TRID explicitly:

```yaml
schema_version: 1
annotation:
  - sample_id: SAMPLE_A
    track: tr
    origin: GENERATED
    producer: tr.trgt
    output_role: vcf
    reference_build: GRCh38
    variant_ids: ["TRID:EPM1_CSTB"]
```

This source-plan fragment requires an explicitly enabled TRGT branch with a
compatible catalog and configured annotation/BED resources. It neither installs
resources nor proves that an arbitrary sample has a called allele. The repeat
overlap excludes the VCF left anchor; raw VCF ID, TRID and source-record hashes
remain distinct in provenance. ANNOVAR/VEP still annotate the entire VCF.

TRGT 的位点名称常保存在 INFO/TRID，而 VCF ID 列为 `.`。TR overlap 请显式使用
`TRID:位点名称`，不修改原始 VCF 来补 ID。混用两种选择方式、TRID 缺失或重复均
报错；普通 SV/BND 的 ID 选择方式不变。重复区间为一基 `POS+1..END`，不包括
前方锚点。零重叠仍保留来源信息；overlap 不代表致病、可信度或真值一致。

执行时先运行 producer DAG；只有本次所需的全部产物、样本、参考及回执通过后，才生成
下游清单并运行第二个 DAG。下游配置关闭 caller，避免重复调用。阶段之间及运行结束
均核对输入；失败不自动重试或丢弃失败样本。

Disabled producers must also be absent from rule registration, not merely from
the final target list. This matters when a ready VCF/gVCF resides inside the
normal output directory: the consumer must use that artifact without trying
to resolve the disabled producer's inputs. Missing ready sources remain errors.
The current source correction covers small/read-SV/TR/phasing/assembly/
assembly-SV/harmonization rules. Later exact-payload regional native handoff
acceptance is recorded in the capability matrix; it does not validate every
route or backend. No calling algorithm or representation is changed.

关闭模块不仅意味着不列入最终目标，也意味着不注册该模块的生产规则。否则下游
读取运行目录内已有 VCF/gVCF 时，调度器仍可能尝试解析已关闭 caller 的输入。
本次源码修复保持缺失来源时报错，不自动重启 caller；外部来源的 harmonization
仍使用独立开关。后续指定候选已完成分项原生衔接验收；不能将某一路径通过推广到
全部后端组合，也不能用此前同版本包的测试代替当前候选验收。

`--dry-run` with a generated-source plan returns `GENERATED_SOURCE_PLAN_ONLY`:
it previews the producer DAG but cannot execute or certify consumers of future
outputs. Real execution uses the same `hifivar --config ... run` entry. Retain
`attempts/<id>/execution_context.json`, `receipts/`, and
`downstream/materialization.json` with their corresponding attempt.

New final reports can verify explicit step receipts and output SHA256. Missing
receipts remain PARTIAL even when Snakemake exits zero. Cross-attempt receipt
reuse is not implemented in this increment: an old output is not stamped as a
new execution, and a failed attempt cannot be hidden by an older success.
Neither a COMPLETE execution report nor a caller support count proves accuracy.

上述回执只证明规定范围内的执行与产物字节一致性，不等同于金标准准确率、人工 IGV
判读或临床结论。caller 的跨 attempt 新调用回执复用尚未实现，不能将历史结果伪装成本次新调用；比对准备的显式验证后复用见下文。

Consumer receipts additionally bind the actual request's VCF/index selection,
the normalized input manifest, effective configuration and materialization
SHA256. In planned mode a consumer can execute only in
`DOWNSTREAM_DAG_RUNNING`, never in the producer stage. Sealing and reading the
receipt recheck that binding. Full byte hashes cover the source roles actually
consumed by the unit; unrelated large producer outputs are not repeatedly read
by each consumer. This is integration/provenance validation, not a sandbox
against arbitrary code in the same Python process.

下游回执进一步绑定真实请求所用的 VCF/索引、标准化输入清单、生效配置和交接记录
SHA256。计划模式下，producer 阶段不能签发下游完成回执；写入与读回都复核来源。
源文件完整字节检查限定于该单元实际消费的角色，不对每个下游单元重复扫描所有大型
gVCF。已有 plain VCF 注释仍可显式使用无索引来源；不自动补建或伪造索引。

For a terminal FAILED attempt, missing or changed recorded artifacts can be
delivered in a FAILED diagnostic report. They never count as validated outputs.
An allegedly successful attempt with the same artifact problem still raises
an error; frozen-control or context/receipt identity violations also remain
fail-closed. Standalone workflows without a launcher context do not acquire a
launcher completion claim. Installed Linux acceptance of this RC5
increment remains required; earlier public assets are unchanged.

终态 FAILED 的任务即使产物缺失，也能交付失败诊断；缺失或改变的产物不计为验收
完成。声称成功的任务遇到同样问题仍严格报错，控制文件、上下文及回执身份校验也
不会放宽。独立 workflow 不会凭空取得 launcher 完成证明。本增量仍需独立 Linux
安装及原生交接验收，不能据此宣称旧公开资产已经包含这些能力。

For the explicit generated harmonization native/derived BND role and required
manifest/receipt bindings, see [BND generated source selection](bnd_junction_consumers.md#generated-bnd-role--本次生成的bnd派生角色).
原生 `vcf` 与 `bnd_vcf` 必须明确选择；派生包不适用于替换队列SV来源或truth benchmark。

## Alignment preparation identity / 比对准备与复用身份

When `alignment.enabled=true`, the launcher creates one sample-specific
preparation receipt in the selected attempt. `ALIGN_NEW` records the observed
mapper and BAM/index validation commands. `VALIDATE_EXTERNAL_REUSE` and
`VALIDATE_PRIOR_ALIGNMENT_REUSE` record fresh validation only: `VALIDATION_ONLY`
is not a new alignment or a new caller run. A previous verifiable producer ID
is preserved; legacy or external origin remains `NOT_PROVEN_LEGACY_OR_EXTERNAL`.
No previous completion record is rewritten during reuse.

The attempt seals `alignment_handoff.json`, including the exact ordered
`calling.samples.tsv`, `calling.yaml`, preparation receipt identities, and
selected BAM/index hashes. Calling configuration may differ from frozen input
configuration only in its explicit sample-sheet path. Consumer step contexts
recheck these small controls; launcher boundaries and alignment report readback
also hash the selected BAM/index bytes. This is not a fresh full FASTQ checksum
audit by every downstream task. The launcher still performs its full original
input checks at stage boundaries. A report distinguishes alignment actions;
files alone cannot establish preparation completion. Legacy attempts without
the handoff remain incomplete for alignment, not retroactively PASS.

启用比对后，每个样本在本次 attempt 中单独签发准备回执。新比对、外部 BAM 验证后
复用、既有比对结果验证后复用明确区分；后两者只证明本次校验，不算新比对或新
calling。旧记录没有可信 producer 身份时保留“来源未证明”，不补写历史成功。
实际交给 calling 的样本顺序、配置、BAM/索引和回执均按 SHA256 绑定；不得通过
修改样本表换入其他 BAM。报告会重新检查选中 BAM/索引，但不宣称重新审计原始
FASTQ，也不代表准确率、相位正确性、人工复核或临床有效性已被证明。
