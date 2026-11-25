# Within-sample SV integration / 同样本 SV 整合

[简体中文](zh_CN/sv_harmonization_run.md)

This RC5 candidate extension connects generated read-SV outputs to
the existing Phase 9 workflow inside `hifivar run`. It does not change any
caller or Jasmine matching parameter defaults. It is included in the
[identified 15 October RC5 payload](rc5_candidate.md#validated-development-payload--已验收开发包)
and has passed a fresh two-sample regional calling-to-harmonization acceptance.
RC1–RC4 and earlier RC5 candidates remain unchanged; select a package using its
exact source manifest and artifact SHA256.

本开发增量把本次 read-SV calling 的输出接到既有 Phase 9：同一个命令可为多个
样本分别完成调用和整合，不把不同样本放进同一次 Jasmine。版本号仍为开发中的
`0.1.0rc5`；历史同版本候选包不因此自动具备本功能。

## Configuration / 配置

Start with the [indexed BAM guide](indexed_bam_run.md) or
[raw-input guide](raw_input_run.md). Add this section to the same YAML, replacing
all executable paths with your provisioned deployments:

```yaml
runtime:
  threads: 8
workflow:
  memory_mb: 64000
sv:
  enabled: true
  overwrite: false
  sawfish: {enabled: true, threads: 4, memory_mb: 16000}
  sniffles2: {enabled: true, threads: 4, memory_mb: 16000}
  pbsv: {enabled: true, threads: 4, memory_mb: 16000}
  cutesv: {enabled: true, threads: 4, memory_mb: 16000}
  harmonization:
    enabled: true
    backend: jasmine
    jasmine_executable: /your/envs/jasmine/bin/jasmine
    truvari_executable: /your/envs/truvari/bin/truvari
    bgzip_executable: /your/envs/htslib/bin/bgzip
    tabix_executable: /your/envs/htslib/bin/tabix
    threads: 4
    memory_mb: 32000
    runtime_minutes: 1440
    max_dist: 1000
    distance_type: linear
    overwrite: false
    input_vcfs: {}
```

Generated read and assembly callers and explicitly selected ready VCFs can be
inputs. Every sample requires at least one source; one source is allowed but
supplies no multi-caller agreement. Generated inputs precede external inputs.
Within each group the order is sawfish, sniffles2, pbsv, cuteSV, PAV, SVIM-asm,
with unselected callers omitted. This order does not depend on YAML key order. `SUPP_VEC` membership follows this actual input-list order; `IDLIST`
retains source-record provenance, not caller-membership inference.

仅启用的 caller 参与整合，顺序固定如上；输入 caller 失败时，下游规则不继续，
不会偷偷忽略它、添加外部来源或制造空 VCF。`threads`/`memory_mb` 必须在总调度预算
内；预算不等于 OS 资源上限。Java 必须与所部署 Jasmine 兼容，并在显式工具环境中
可用。Jasmine 输出 VCFv4.4 时需 tabix/HTSlib ≥1.23.1；HiFiVar 不降级 header。

The existing Phase 9 rule also runs Truvari comparisons between the harmonized
VCF and each source. These are concordance comparisons, not GIAB benchmarks.
They require Truvari even when `benchmark.enabled` is false. Tool versions are
probed through the existing wrappers. With probes enabled, preflight rejects
tabix below the shared VCFv4.4 compatibility floor of 1.23.1 before analysis;
it does not execute harmonization or certify the complete Java/HTSlib deployment.

## Check and execute / 检查与执行

```bash
hifivar --config analysis.yaml doctor --preflight --run-dir /your/project/runs/sv_batch01 --probe-tools --json
hifivar --config analysis.yaml --dry-run run --run-dir /your/project/runs/sv_batch01
hifivar --config analysis.yaml run --run-dir /your/project/runs/sv_batch01 --resume
```

With alignment disabled, dry-run shows the calling/harmonization DAG after
read-only BAM checks. With alignment enabled it still returns
`RAW_INPUT_PLAN_ONLY`; the downstream DAG needs real aligned BAMs first.
Resume requires the same input/config/runtime identity. No retry, overwrite, cross-sample pooling or source discovery is automatic.
External inputs must be declared explicitly as below.

同一份样本表中，每个样本独立产生一套结果；并非 SV cohort 联合分型或跨样本聚类。
每个来源必须显式选择。实际启用 assembly/assembly-SV 时仍需要独立 FASTQ
和已配置的工具部署；只使用现成 assembly-SV VCF 时不要启用 assembly caller。
benchmark 与 review 不由本命令自动执行。

## Outputs and interpretation / 输出及解释

For each sample, retain:

- `results/sv/<sample>.<caller>.sv.vcf.gz` and its index, unchanged source calls.
- `results/sv_harmonized/<sample>/<sample>.harmonized.sv.vcf.gz` and `.tbi`.
- `results/sv_harmonized/<sample>/<sample>.sv.evidence.tsv`.
- `results/sv_harmonized/<sample>/<sample>.phase9.provenance.json`.
- `results/sv_harmonized/<sample>/truvari/<caller>/` source comparisons.

The evidence table records caller membership and read-source class, not call
confidence, pathogenicity or truth. SNV/Indel, SV and TR VCFs remain separate.
Raw caller VCFs are not rewritten; Jasmine's derived raw/finalized outputs and
the input list remain separate provenance. A successful DAG proves execution
within its specified scope, not biological accuracy or wet-lab validation.

结果应连同 launch manifest、有效配置、来源顺序、工具版本和日志保留。caller 支持数
不是 read 数，也不是准确率；Truvari 来源比较不是匹配金标准的独立 benchmark。
这些运行文件通常含样本名和绝对路径，默认保存在私有项目目录，不自动公开。

## Optional BND bundle in development source / 开发源码中的可选 BND 数据包

The later development source can explicitly enable the separate Phase 9 BND
representation bundle. It is default-disabled, per-sample and no-overwrite.
Native VCF/evidence outputs and existing concordance inputs retain their roles;
derived VCFs are not silently selected by other consumers. See the exact
[policy and YAML contract](bnd_representation.md#phase-9-opt-in--显式启用).
This is source support, not a claim that RC1–RC4 assets contain the option or
that BND accuracy has been validated.

后续开发源码可显式启用独立 BND 派生包。默认关闭，按样本隔离、禁止覆盖，不改变
原生输出和既有 concordance 输入。历史 Release 不因此获得该功能，不能将表示兼容性
写成 BND 准确率。完整数据包而非单独派生 VCF 才具有原生记录可逆性。


## Explicit ready VCFs in current development source / 当前开发源码的现成来源

This later, uncommitted development increment extends the launcher to the existing
Phase 9 external-input contract; a version string alone does not identify it.
Historical releases and earlier RC5 packages are unchanged. Its native acceptance
must be checked against the exact candidate manifest, not the earlier regional
read-calling acceptance quoted above.

For harmonization only, disable calling and list every sample's sources:

```yaml
sv:
  enabled: false
  harmonization:
    enabled: true
    overwrite: false
    input_vcfs:
      HG002:
        sawfish: /your/project/ready/HG002.sawfish.sv.vcf.gz
        pav: /your/project/ready/HG002.pav.assembly.sv.vcf.gz
        svim_asm: /your/project/ready/HG002.svim_asm.assembly.sv.vcf.gz
```

Keep the reference, sample sheet, tool paths and scheduling budget configured as
above. Sources must be single-sample BGZF VCFs with a readable conventional
`.tbi`, exact sample identity and matching declared contig lengths/build.
Unknown samples/callers and a caller both generated and external fail closed.
Sources are never renamed, filtered, normalized or silently discovered.
Missing/unrecognized build metadata remains unverified; matching names and
lengths do not prove reference-allele sequence identity.

现成来源按样本显式配置，严格核对样本名、染色体名称/长度和已声明的参考版本。
缺失版本信息不能视为已证实匹配。文件名遵循样本与 caller 的既有规则，原始文件和
索引均纳入启动身份；不自动添加、删减或混合不同样本。

A harmonization step binds its actual canonical VCF/index pairs to frozen launch
inputs (external) or exact same-attempt producer receipts (generated). The
completion receipt retains source order/origin, observed commands and output
identities; source bytes are checked before execution, at seal and on receipt
readback. A synthetic test context is not a native launcher validation.
Runtime fingerprints include JSON policy resources as well as Python, YAML and
workflow bytes. Changing a policy or source requires a new run directory.

执行回执绑定实际使用的来源，而非只登记输出文件。可选来源身份视图另见
[SV source identity](sv_source_identity.md)：仅改变重复/缺失 ID 的显式派生视图，
不重写 canonical caller VCF，也不把支持数解释成准确率或置信度。

## Source assembly roles are not variant phase / 来源组装角色不等于变异相位

An external PAV/SVIM-asm filename does not establish that haplotype1 and
haplotype2 were both inputs. Without a bound provenance record, the bridge
retains empty source roles and UNVERIFIED_ASSEMBLY_INPUT_ROLES in source JSON;
raw FORMAT/GT and other VCF fields remain untouched. Empty role metadata is
unknown, not absence of a variant or evidence for an unphased biological sample.

For current-attempt generated assembly-SV sources, bridges bind the actual
request's role/FASTA identities to the same-attempt assembly receipt, recheck
them after execution and retain that lineage in the assembly-SV receipt.
Harmonization resolves only these bound roles and rechecks the lineage at seal
and receipt readback. Missing, swapped or changed inputs/receipts fail closed.
The evidence TSV's legacy haplotypes column describes source input roles, not
per-variant haplotype support, genotype correctness or biological phasing accuracy.
Standalone workflow/API artifacts can still carry explicit declared roles, but
a missing attempt binding cannot be invented from a filename or phased GT.

现成组装来源若缺少可核验角色绑定，不再自动填入“两套单倍型”。来源 JSON 明确
保留未知；VCF 原始 GT/phase/allele 不改。本次生成来源则关联真实 assembly
回执及角色对应的 FASTA SHA，文件或回执在执行中、完成后被替换均拒绝。
这里的角色是来源输入描述，不是某个变异同时被两套单倍型支持或已通过相位真值验证。

## Optional caller-support selection / 可选支持数筛选

Current development source accepts `sv.harmonization.minimum_caller_support`
(default 1). A value above 1 produces a separate per-sample support-ge-k bundle
after full native harmonization; complete native outputs and existing consumers
keep their roles. See the [bilingual configuration and output contract](sv_support_selection.md).

当前开发源码可设置最低 caller 支持数；完整结果始终保留，筛选输出另存，不能把
支持数解释成置信度。历史候选包与 Release 不因此自动具备此选项。
