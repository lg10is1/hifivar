# FASTQ and native uBAM batch execution / 原始数据批量运行

[简体中文](zh_CN/raw_input_run.md)

This development-source extension adds opt-in alignment preparation to
`hifivar run`. It is included in the
[identified 15 October RC5 payload](rc5_candidate.md#validated-development-payload--已验收开发包),
**not in RC1–RC4 or the initial security-only RC5 package**. Use the
delivery manifest and artifact hashes; a version string alone is insufficient.

这是 RC5 候选功能。上述指定开发包已完成相应的受限 Linux 验收；不能根据 `0.1.0rc5`
版本号认定所有旧候选包都有此入口。源码或文档改动后的新包仍需自己的交付身份与验收。

## Execution contract / 执行边界

One command performs two explicit stages: sequential alignment preparation,
then the existing multi-sample calling DAG. This is **not a single raw-input
Snakemake DAG** and does not submit distributed Slurm jobs. It accepts:

| Input | Action |
| --- | --- |
| HiFi FASTQ, including gzip and ordered files separated by `;` | pbmm2 or minimap2 alignment, then BAM validation |
| Native PacBio uBAM, `input_type: bam` in TSV | Header/record checks, native pbmm2 alignment; no FASTQ conversion |
| Coordinate-sorted BAM with one BAI/CSI | Strict reference/sample/index checks, then reuse unchanged |
| CRAM or aligned BAM without coordinate sorting/index | Explicit error; no automatic repair |

uBAM requires consistent `@RG` IDs, every `SM` matching `sample_id`, `PL:PACBIO`,
no `@SQ`, non-coordinate sorting, and nonempty primary, unpaired, unmapped
records. Missing/undeclared record RG tags fail. These checks do not certify
HiFi read accuracy. pbmm2 receives the uBAM directly without HiFiVar rewriting
read groups or discarding tags through FASTQ conversion. Tag retention by an
external tool is not a general methylation-analysis validation claim.

uBAM 不隐式转换为 FASTQ，以免丢失 PacBio 附加标签；不自动改 `SM`。因此样本表必须使用
原数据的真实样本名。已有 BAM 不重比对、不补索引；FASTQ 的样本身份由样本表声明。
原始文件、参考序列和旧分析结果保持只读。调用仍可选 small、read-SV、TRGT 和 HiPhase，
不自动启用组装、队列、注释或 benchmark。最新开发增量允许在同一配置中显式
开启[同样本 read-SV 整合](sv_harmonization_run.md)；默认仍关闭，不进行跨样本聚类。

## Configuration / 配置

Follow [indexed BAM execution](indexed_bam_run.md) for the calling branches.
Add the following to the same user YAML, replacing all example paths:

```yaml
runtime:
  threads: 8
workflow:
  memory_mb: 64000
  samtools_executable: /your/envs/samtools/bin/samtools
alignment:
  enabled: true
  tool: pbmm2
  pbmm2_executable: /your/envs/pbmm2/bin/pbmm2
  minimap2_executable: /your/envs/minimap2/bin/minimap2
  output_format: bam
  threads: 8
  pbmm2_sort_threads: 1
  index_threads: 2
  memory_mb: 32000
  runtime_minutes: 1440
  overwrite: false
  pbmm2_preset: CCS
  bam_index_format: auto
```

The pbmm2 mapping thread count is `threads - pbmm2_sort_threads`; sorting gets
the explicit remaining allocation. Indexing runs separately. Memory is a
planning budget, not an OS-enforced limit. External tools may start auxiliary
threads. Samples are aligned sequentially to avoid multiplying this footprint.
minimap2 uses its existing `map-hifi` wrapper and supports BAI only here; select
pbmm2 when CSI is required. Native uBAM plus `tool: minimap2` fails explicitly.

TSV example (real tab separators):

```text
sample_id	input	input_type	sex
FASTQ_SAMPLE	/your/data/part1.fastq.gz;/your/data/part2.fastq.gz	fastq	female
UBAM_SAMPLE	/your/data/unaligned.bam	bam	male
ALIGNED_SAMPLE	/your/data/aligned.bam	bam	female
```

## Plan, run and resume / 预检、运行与续跑

```bash
hifivar --config analysis.yaml --dry-run run --run-dir /your/project/runs/raw_batch01
hifivar --config analysis.yaml run --run-dir /your/project/runs/raw_batch01 --resume
```

With `alignment.enabled: true`, dry-run returns `RAW_INPUT_PLAN_ONLY`: it reads
and hashes inputs and inspects BAMs, writes routing evidence, but **does not run
alignment or claim that the downstream calling DAG has been dry-run**. The
calling DAG can only be constructed after real indexed BAMs exist. No empty
placeholder BAMs or invented indexes are created. Omitting dry-run executes
both stages in one invocation at a new run directory.

启用比对后的 dry-run 是输入路由预检，不是完整 calling DAG 验收。查看
`attempts/attempt_*/alignment_plan.json` 后执行同一配置。续跑会校验原始输入、配置、
解释器和源码身份；已完成的比对还会核验 BAM/索引完整 SHA256 并重新检查格式。失败的
比对仅在显式指定 `--resume --rerun-incomplete` 时，才进入新的独立尝试目录；历史尝试
不覆盖、不删除。手动处理中断锁前必须确认所有相关进程已停止。

## Outputs / 输出

- `results/alignment/<sample>/attempt_<id>/`: generated BAM/index and per-command evidence.
- `results/alignment/<sample>/completed.json`: successful alignment receipt and full checksums.
- `calling.samples.tsv`, `calling.yaml`: derived handoff to the existing calling DAG.
- `attempts/attempt_<id>/status.json`: invocation outcome, keeping failures separate.
- Existing per-class calling outputs remain under `results/`; they are never forced into one VCF.

`WORKFLOW_SUCCEEDED` is engineering execution success, not truth-set accuracy,
whole-genome coverage, or success of disabled branches. This extension does not
change the scientific algorithms in pbmm2/minimap2/callers. Keep sample IDs,
absolute paths and command evidence private. For native-format details see the
[pbmm2 documentation](https://github.com/PacificBiosciences/pbmm2).
