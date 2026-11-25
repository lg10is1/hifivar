# Existing haplotype FASTAs: explicit reuse contract

[简体中文](zh_CN/ready_assemblies.md)

## Status and scope

RC5 supports explicit ready-assembly input in
`hifivar run`, as well as the standalone wrapper API. It is not part of RC1–RC4
or earlier RC5 artifacts: identify a candidate by its complete source manifest
and package hash, not its version string alone. The input-validation record is
distinct from a new native assembly-SV receipt. Do not enable hifiasm with dummy
reads, substitute an unrelated FASTA at a generated output path, or manufacture
a native assembly receipt to bypass workflow guards.

The original FASTQ → hifiasm → PAV/SVIM-asm route is unchanged. This addition
does not alter calling, normalization, filtering, benchmark policy or phasing.

## One invocation, explicit input mode

Keep the real sample TSV/reference used for the project. Its existing FASTQ or
BAM input is an identity context, **not an instruction to analyze those reads**
when alignment, small, read-SV, TR and phasing are disabled. No dummy reads are
needed or accepted. The primary sample model remains FASTQ/BAM/CRAM; this does
not add a FASTA-only primary sample-sheet format. CRAM remains deferred.

Use an independent ready TSV with the exact same sample-ID set, then select:

```yaml
assembly:
  enabled: false
assembly_sv:
  enabled: true
  input_manifest: ready.tsv
  overwrite: false
  pav:
    enabled: false
  svim_asm:
    enabled: true
    threads: 2
    memory_mb: 32768
```

This excerpt overlays a complete project configuration: set the real reference,
sample table, scheduler capacity and installed tool paths separately. Resource
values must fit the configured global budget and be appropriate to the input;
the example is not a whole-genome minimum-memory recommendation. To run PAV,
explicitly enable it and configure its supported 2.4.6 version, real Snakefile
and adapter/finalizer executables. Do not select an unresolved placeholder.

```bash
hifivar --config analysis.yaml run --run-dir /your/project/runs/ready01
hifivar report --run-dir /your/project/runs/ready01 --output-dir /your/project/reports/ready01
```

Use fresh authorized directories. Report generation is a separate invocation,
not an automatically included DAG step. A repeat/resume gets a new attempt;
ready caller work directories are attempt-specific. Existing caller results do
not acquire new execution receipts merely because they still exist. Missing
assembly-SV receipts remain PARTIAL even if a scheduler exits zero.

The launcher freezes the original manifest and every selected FASTA identity.
It writes `ready_assembly_input_validation.json` in the current attempt and
seals that record's identity into `execution_context.json`. This is
`VALIDATION_ONLY`, not a native hifiasm step. PAV/SVIM-asm scripts bind actual
FASTA roles to this input validation before and after new native execution.
The final report keeps `assembly` DISABLED, reports `assembly_input_reuse`
separately, and certifies a caller only through its new validated native
receipt. Producer inputs remain frozen when a later generated annotation or
cohort stage disables calling; the consumer configuration clears the ready
manifest instead of reactivating assembly-SV.

For annotation or category-specific cohort summaries, a downstream source plan
can select the **new caller output** with `origin: GENERATED`,
`producer: assembly_sv.svim_asm` (or `assembly_sv.pav`), `track: sv`, and
`output_role: vcf`. Here GENERATED means a new PAV/SVIM-asm output from this
attempt, not a new hifiasm assembly. Keep the selected caller enabled and choose
exactly one assembly input mode: enabled hifiasm, or the explicit ready manifest
with hifiasm disabled. Intent parsing does not validate a missing ready TSV or
certify an output in advance; the launcher and current-attempt native receipt
checks are still required. SV cohort tables retain source-native records and
genotypes; they do not perform cross-sample SV clustering or joint genotyping.

## Manifest

A UTF-8 (optional BOM) tab-separated file requires exactly these four columns:

```text
sample_id	role	fasta	sha256
S1	haplotype2	assemblies/S1.h2.fa	<64-hex SHA256 of the actual file>
S1	haplotype1	assemblies/S1.h1.fa	<64-hex SHA256 of the actual file>
```

Optional columns: `source_tool`, `source_version`. Values are **user-declared**,
not command-probed versions. Omitted origins remain unknown. `fasta` paths are
resolved relative to the manifest, not the working directory. Only readable,
uncompressed `.fa`, `.fasta`, `.fna` inputs are accepted; no automatic conversion
or chromosome renaming occurs. Compute each expected checksum independently;
the placeholders above are not runnable values.

One or two explicit roles per sample are supported, matching the existing
haploid/diploid wrapper contract. Two roles require distinct physical files;
same-file paths, symlinks and hard-link aliases cannot represent different
roles or samples. The declared sample set must equal the caller's expected
sample set. Empty/unsafe IDs, unknown columns, duplicate roles and missing or
extra samples fail. An 8 MiB manifest / 10,000-row safety limit applies.

Every FASTA byte is hashed and its complete structure is checked using bounded
chunks, including unwrapped chromosome-length sequence lines. Headers require
unique nonempty ASCII identifiers and at most 64 KiB each. Each record must
contain IUPAC DNA; blank lines, sequence whitespace, empty records, duplicate
IDs and invalid symbols are rejected. LF, CRLF and a missing terminal newline
are supported. Contig identifiers and sequence bytes are never rewritten.

## Library use

```python
from pathlib import Path
from hifivar.ready_assemblies import load_ready_assembly_manifest

ready = load_ready_assembly_manifest(Path("ready.tsv"), expected_samples=("S1",))
assemblies = ready.for_sample("S1")
# Pass assemblies to the existing AssemblySvRequest for a standalone wrapper.
# Select an independent caller work directory and explicit resources as usual.
```

Declared row and sample order are retained. `for_sample` provides h1 then h2
argument order so diploid caller inputs do not depend on TSV row order. This
ordering is not evidence that individual variants belong to either haplotype.

The assembly input dictionaries record `EXTERNAL_READY`, `VALIDATION_ONLY`,
the expected SHA256, actual byte/sequence/base counts, and
`USER_DECLARED_NOT_PHASE_TRUTH`. `source_gfa`, `hifiasm_version` and
`reference_independent` are null rather than fabricated. No new hifiasm
command, version observation or native producer receipt is generated.

PAV/SVIM-asm standalone requests can use these artifacts. A request cannot mix
generated hifiasm and externally ready artifacts. Full FASTA/manifest content
is rechecked before the caller is planned and before its output is accepted;
manifest sample/role/source fields are bound as well as file bytes. This is
not a filesystem sandbox against adversarial concurrent writers: callers
should retain read-only inputs and isolated task directories.

In ready mode, caller work directories must be new or empty, not reused native
workspaces, and must not contain the source FASTAs or manifest.
Output/index/log aliases of those inputs are rejected before tool execution,
including hard links, existing symlink aliases and the stdout log automatically
derived by `CommandRunner` from a stderr-only log path. Native work, logs and final
outputs must use separately owned paths. The caller still executes through
its existing `CommandRunner` and output validation. A successful **new caller
run** is distinct from **validation-only input reuse** and must be reported as
such. No raw input is copied, moved or altered by the validation layer.

## 中文说明

当前开发源码已接入一键启动器的现成组装模式，同时保留独立 wrapper API。
这不是 RC1–RC4 或所有早期 rc5 包都有的功能，必须核对具体候选的完整源码清单
与包哈希。保持真实的样本表／参考上下文，设 `assembly.enabled: false`，并用
`assembly_sv.input_manifest` 显式选择 ready TSV；二者不能与新 hifiasm 组装
同时启用。未增加“FASTA 作为主样本表 input_type”的格式，不得创建占位 reads。

启动器把原清单、FASTA 哈希与本次 attempt 绑定，输入校验记录只标记
`VALIDATION_ONLY`。真正完成新调用后才生成 PAV／SVIM-asm 各自的 native 回执。
报告中 `assembly` 仍为 DISABLED，`assembly_input_reuse` 单独呈现验证后复用；
缺少 caller 回执时，即使调度器退出为0也不能写成调用完成。报告需单独执行，
不是一个命令内自动完成的 DAG 步骤。旧结果不能自动获得新 attempt 的执行身份。

注释或分类型汇总的 downstream source plan 可用 `origin: GENERATED`，选择
本次新执行的 `assembly_sv.pav`／`assembly_sv.svim_asm` 的 `vcf` 输出。
这里 GENERATED 指新 caller 产物，不表示重新执行 hifiasm。组装输入仍须二选一：
新 hifiasm，或显式 ready 清单；所选 caller 必须启用。意图解析不等于输入验收，
后续仍核对完整清单、FASTA 和当前 attempt 的真实 native 回执。SV 汇总表保留
来源位点和基因型，不是跨样本 SV 聚类或联合分型。

清单逐行声明样本、`haplotype1`／`haplotype2`、真实 FASTA 路径及完整 SHA256。
相对路径以清单目录为基准；样本集合必须完全一致。角色／样本不可共用同一物理
文件。清单顺序保留，调用参数固定按 h1、h2 排列；这些是用户声明的输入角色，
不是生物分相真值，更不是逐变异的单倍型支持。

校验器流式检查整个 FASTA，不只检查首条记录。原始 FASTA、contig 名称、序列
都不修改。不存在的 GFA、工具版本和 reference-independent 属性保留为空。
输入复用标记为 `VALIDATION_ONLY`，不写成新执行了 hifiasm。新的 PAV／SVIM-asm
调用若真实完成，应单独记录其新执行；输入验证不能替代工具版本和输出验收。

工作目录不得包住原输入，输出、索引、日志不得通过路径、软链或硬链覆盖输入。
所有 scientific accuracy、最低深度、成本和分相准确率结论均不由此功能建立。
