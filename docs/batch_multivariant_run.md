# One-invocation batch small/SV/TR analysis / 一次调用的批量多类型分析

[简体中文](zh_CN/batch_multivariant_run.md)

## Availability and scope / 可用范围

This guide targets **RC5 release source**, including generated-source
handoff and optional joint HiPhase. Follow the [RC5 installation guide](installation.md)
and check the SHA256 of the downloaded asset. Earlier 0.1.0rc5 development builds
have different manifests and validation scopes. See [release identity](rc5_candidate.md).

本指南使用 GitHub RC5 发行包，涵盖调用结果交接及可选的联合定相。请按[安装指南](zh_CN/installation.md)
核对包的 SHA256。历史同版本开发包以各自的来源清单和验证记录为准。

The complete [example](../examples/batch_multivariant/README.md) starts from one
coordinate-sorted indexed BAM per sample. One analysis invocation executes
selected producers then a distinct consumer DAG. It is not one all-module DAG:

1. DeepVariant small VCF/gVCF; Sawfish, Sniffles2, pbsv and cuteSV VCFs; TRGT VCF.
2. Per-sample Jasmine harmonization/evidence and joint small/SV/TR HiPhase,
   with independent output categories.
3. GLnexus from this attempt's gVCFs; separate SV/TR source/state summaries.
4. After acceptance, generate the final report with a separate explicit command.

中文：从已有BAM批量分析；small、SV、TR始终分开。一条分析命令按配置串接producer和
consumer阶段，不把不同DAG或不同历史任务伪装成一个全模块任务。报告仍单独生成。
本模板不启用组装、注释、benchmark、AlphaGenome或IGV；启用其他分支需对应独立资源和
指南。原始FASTQ/uBAM先使用[可选准备阶段](raw_input_run.md)；主入口CRAM仍延期。

## 1. Install and provision / 安装与资源

Install an explicitly identified matching candidate with its workflow extra,
using [installation instructions](installation.md). The wheel does not bundle
examples, caller executables, SIF images, reference, TR catalog or databases.
Obtain examples from the matching trusted source/sdist, not from an unrelated
newer checkout. Confirm CLI/workflow resources and provision every enabled tool.

- Reference FASTA + FAI: declared GRCh38, exact contig names and lengths matching BAMs.
- BAM + unique BAI/CSI: coordinate sorted, all RG IDs valid, all SM equal to the sample ID.
- DeepVariant: pinned PACBIO container already present; Apptainer available.
- Four read-SV callers, Jasmine/Truvari and their Java/runtime dependencies.
- TRGT/bcftools/samtools and a pinned compatible TR catalog shared across samples.
- HiPhase, GLnexus, bgzip and tabix. This Jasmine VCFv4.4 route requires tabix >=1.23.1.

中文：各工具可用独立Conda环境，但YAML中的executable必须指向正确入口；不要把多个
环境的bin目录无序堆入PATH。Jasmine的Java和VEP的Perl等解释器依赖要独立验证。
该模板不启用VEP，不应为了运行此例加载任何API凭据。

## 2. Prepare three files / 准备三份文件

Copy config.yaml, samples.tsv and source_plan.yaml from examples/batch_multivariant
into a new analysis directory under your project root. Replace every /work path
with your own reference, catalog, container and BAM locations. The sample TSV
contains real tabs, not literal backslash-t sequences.

The example declares SAMPLE_A female and SAMPLE_B male only to demonstrate
metadata. Replace both with actual declared sample metadata. karyotype:auto uses
that metadata; it is not inferred from reads. Unknown/nonstandard karyotype needs
an explicit supported configuration, not guessed sex. Every BAM RG SM must match.

For each sample keep exactly three cohort source selectors:

| Track | Producer | Native output role | Meaning |
|---|---|---|---|
| small_variants | small.deepvariant | gvcf | GLnexus input, not ordinary small VCF |
| sv | sv.harmonization | vcf | That sample's harmonized SV source |
| tr | tr.trgt | vcf | TR source with the configured catalog content identity |

GENERATED is intent, not proof of execution. Do not add future output file paths,
copy old receipts or declare old VCFs as newly called. The launcher verifies native
producer receipts, bytes, sample/build, frozen controls and category before
materializing the downstream manifest. Generated TR selectors acquire the actual
configured catalog SHA256; an optional explicitly supplied hash must match.

中文：增删样本时同时修改TSV和source_plan的三行，不能只改样本表。
GLnexus必须接gVCF；定相后的VCF不会静默替代caller来源。

## 3. Budget and preflight / 资源与预检

The supplied scheduler example uses 8 total cores, a 256-GiB memory budget,
4 DeepVariant shards and one DeepVariant sample slot; GLnexus requests and caps
192 GiB. These are examples, not measured requirements or guarantees for your
cohort. Change all relevant per-rule allocations consistently, after checking
node availability, UID/cgroup thread limits, temporary disk, quotas and inputs.
The scheduler budget is not an operating-system memory/thread limit.

中文：不能因示例写了192GiB就认为所有队列都够，也不能把旧小区域8GiB配置照搬WGS。
默认small.region:null表示不主动限制DeepVariant到区域；只在有界测试时显式设置区域，
并使用独立区域BAM/catalog。区域参数不限制SV caller，不能以此冒充全基因组输出。
Sawfish CNV默认保留，不照搬稀疏区域案例中关闭CNV的权宜参数。

Inspect configuration without analysis, then opt in to native metadata probes:

~~~bash
hifivar --config config.yaml config validate
hifivar --config config.yaml doctor --preflight --run-dir /your/project/runs/batch01 --json
hifivar --config config.yaml doctor --preflight --run-dir /your/project/runs/batch01 --probe-tools --timeout-seconds 30 --json
~~~

The preflight run directory must not exist. Metadata-only preflight normally
returns INCOMPLETE/exit2 for unprobed scope; foundation doctor PASS is not caller
readiness. A preflight probe does not start containers/calling or certify every
cohort backend. Read each NOT_CHECKED field and provision consumer dependencies
separately. See [preflight details](preflight.md).

## 4. Plan, run, report / 计划、运行与报告

Use a new stable run directory within your authorized project root. Leave
paths.workdir/outdir, runtime.tmpdir and logging.file null; launcher-owned
working paths stay beneath run-dir. Launch from a credential-free environment.

~~~bash
# Reads real inputs and writes plan metadata; no calling.
hifivar --config config.yaml --dry-run run --run-dir /your/project/runs/batch01

# After reviewing the exact plan, execute the same frozen controls.
hifivar --config config.yaml run --run-dir /your/project/runs/batch01 --resume

# Use the ACTUAL immediate-child attempt ID recorded by that invocation.
hifivar report --run-dir /your/project/runs/batch01 --attempt-id attempt_ACTUAL_ID --output-dir /your/project/reports/batch01
~~~

Do not copy attempt_ACTUAL_ID literally. A generated-source dry-run may report
GENERATED_SOURCE_PLAN_ONLY: future consumer inputs are not yet certified. Config,
input or runtime changes require a new run. Resume after a plan is not permission
to retry a failed analysis automatically. On failure preserve logs/partials and
check scheduler/descendant identity before any separately authorized remediation.

The launcher is foreground/local-executor. Use a site-approved persistent session
or submit one appropriately resourced single-node Slurm job. It does not submit
a distributed Slurm DAG. Do not simultaneously launch duplicate schedulers or
promise survival after SSH loss without a tested persistent/batch mechanism.

中文：先计划，审核后同目录resume；report选本次真实attempt。SSH断开不等于进程终止，
必须核对原进程身份，不重复启动。失败不默认删文件、重启、改参数或自动重试。

## 5. Read the outputs correctly / 正确解读结果

| Output family | Location below results | Interpretation |
|---|---|---|
| Small VCF/gVCF | small/ | Per-sample calls/reference blocks are distinct |
| Four caller VCFs | sv/ | Caller-specific records overlap; do not sum into unique SVs |
| Harmonized SV/evidence | sv_harmonized/<sample>/ | Within sample; support is descriptive |
| TR | tr/ | Catalog-dependent repeat genotypes, independent category |
| Phased small/SV/TR | phasing/ | Separate VCF/index per selected source; not new calls |
| Joint small variants | cohort/COHORT01/small/ | GLnexus; declared and observed orders retained |
| SV/TR summaries | cohort/COHORT01/sv/ and tr/ | Source/GT/state tables, not new cross-sample calling |

Preserve launch_manifest.json, effective.yaml, samples.tsv and the actual
attempt's execution_context.json, receipts and downstream/materialization.json.
WORKFLOW_SUCCEEDED does not by itself establish every output or accuracy;
validate artifacts, native receipt source/command identities and input hashes.
Report keeps missing/failed/planned/executed states separate.

中文：不要把不同caller计数、phased/unphased文件或重复证据相加成“总变异数”。
SV汇总不是跨样本聚类，TR汇总不是联合分型；missing不填成hom-ref。GLnexus可能改变
物理样本列顺序，必须按名称映射QC，不能用列位置猜测身份。

## 6. What has actually been demonstrated / 已验证的范围

An independently checked installed-source case completed **one analysis invocation**
for two technical HG002 labels: fresh small variants, four read-SV callers, TRGT,
within-sample Jasmine, joint small/SV/TR HiPhase and all three cohort tracks.
It used bounded regional BAMs and one catalog locus; biological N=1, not two
independent people. Native joint phasing was also checked in a separate single-label
case. These executions used an identified earlier source snapshot; they do not
certify every later same-version package, all tools jointly at WGS scale, speed,
sequencing-cost savings, phasing accuracy, TR accuracy or independent clinical validity.

中文：这些是真实工程执行证据，但范围是特定包、区域输入和技术标签。WGS small与SV
金标准研究另有固定truth/BED/过滤/分母和数据冻结，不能从本例计数推导准确率。
较早的两项独立批量案例仍是历史独立运行；新同次案例不能倒写其身份。

See [capabilities](capabilities.md), [joint phasing](joint_phasing.md),
[within-sample integration](sv_harmonization_run.md), [benchmarking](benchmarking.md)
and [VEP output/deployment safety](vep_output_safety.md). Human review remains
NOT_REVIEWED without reviewer/date/decision/comment; no screenshot or software
execution replaces that evidence.
