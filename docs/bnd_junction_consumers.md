# BND junction-aware annotation and review / 双端断点注释与复核

[简体中文](zh_CN/bnd_junction_consumers.md)

This is an RC5 release source increment. It does not alter caller
algorithms, redirect an existing benchmark or establish BND accuracy.

这是 RC5 中的独立后续处理入口，不改变caller、不追认历史结果、不建立BND准确率。

## Coordinates and identity / 坐标与身份

A paired BND is consumed as two VCF anchors parsed from CHROM/POS/REF/ALT, with
its exact bracket orientations retained. END is not a local event-spanning
interval. A single breakend has one known endpoint and an unknown mate.
Coordinates 0 and N+1 remain virtual anchors, not ordinary reference bases.

双端分别注释与截图；不会将两个断点间的基因全都称为受影响基因。普通坐标按一基闭区间
与零基半开BED比较。虚拟端点不参与普通BED重叠、不被平移为1或N，也不生成虚构截图；
它们仍出现在结果JSON的端点状态中。方向为现有解析器的local_side/mate_side及shape编码，
不推断基因融合、致病性或相位准确率。

Only explicitly selected record IDs are materialized. The complete VCF is
streamed; duplicate selected IDs fail as ambiguous. Record ordinal, source and
raw-record SHA256, REF/ALT, coordinate semantics and declared MATEID QC remain
traceable. Missing or duplicate mate IDs are retained as UNRESOLVED or
AMBIGUOUS_DUPLICATE_ID; ALT mate coordinates are not invented record links.
Coordinate inequality alone does not prove a biologically invalid mate.

## Region-overlap input / 区域注释输入

The ordinary five-column annotation manifest remains supported:

`sample variant_category source_vcf source_tool source_variant_ids`

The actual format is TSV. For an explicit BND-derived VCF, append exactly one
optional column, `bnd_bundle_manifest`, after those five columns. Use `sv` as
the category, the exact derived VCF as source_vcf, and the matching bundle's
manifest.json. Relative paths are interpreted beside the input manifest.

Do not pass an isolated derived VCF. The consumer verifies the pinned policy,
mandatory native/derived/index/ledger/junction files, exact reference identities
and a complete native-record round trip. It rejects changed/missing files or a
bundle attached to a different VCF. The launcher freezes the bundle files and
binds the explicit role to its attempt receipt; readback revalidates them.
This proves internal consistency, not cryptographic authenticity against a
malicious party changing the entire bundle.

启用现有annotation.overlap_enabled及版本化BED配置即可；不新增隐式选源策略。
输出TSV每个重叠带LOCAL/ALT_MATE端点、方向、来源键和junction_json；没有命中的端点、
单端未知和虚拟端点保存在Phase11 provenance的junction_annotations中。普通SV/TR跨度
行为保留，small/SV/TR仍分开。ANNOVAR与VEP不会因这一增量自动改用派生文件。

## Review selection / 复核清单

Use the existing review prepare entry and TSV selection. For BND, contig/start/end
must match the exact source CHROM/POS/POS. Optional mate_contig/mate_position may
be omitted: the known ALT endpoint is then read directly from that selected
record. If supplied, both must exactly match ALT. Single breakends must not
attach an invented mate. Append optional bnd_bundle_manifest only for an explicit
verified BND-derived source.

复核清单会保留两端方向、来源记录与未解析配对QC。普通端点生成独立IGV窗口，显示窗口
限制在参考内；原始断点坐标不变。单端只生成一端截图计划。虚拟端点不生成普通goto，
并明确记录为NO_REFERENCE_BASE_VIRTUAL_ANCHOR；未知端记录UNKNOWN_SINGLE_BREAKEND。
准备批处理文件不代表已运行IGV，更不代表已完成24个人工复核。

All initial decisions remain NOT_REVIEWED. Reviewer, timestamp and comment are
still required for imported decisions; declarations do not rewrite a VCF or
become authenticated truth. Existing output and source bytes remain protected.

## Limitations / 限制

No fusion, pathogenicity, BND truth accuracy, or assembly phase-accuracy inference.
No automatic lifting/renaming of contigs, reconstruction of missing mate IDs,
filtering, or canonical VCF rewrite. Native and derived roles require explicit
selection; the source plan never chooses a representation automatically. Large
review selections currently repeat source validation per target; use a bounded
explicit review list rather than all WGS calls. Generated-source selection
is described below and does not imply scientific accuracy.

Grammar follows the pinned [VCFv4.4 specification](https://github.com/samtools/hts-specs/blob/c101c797106d8a1006b9b06f75e8bcc6062cf81e/VCFv4.4.tex),
not a new rearrangement caller. See [BND bundle policy](bnd_representation.md).

## Generated BND role / 本次生成的BND派生角色

The RC5 source plan can explicitly select the native role `vcf` or
`bnd_vcf` from `sv.harmonization`. The native role remains unchanged; the user
must actually declare one role. There is no automatic native/derived swap.
The derived selector requires the existing exact enabled BND policy, including
its policy id/version/hash/mode. It is an SV annotation source, not a cohort
SV replacement or benchmark redirect. Configure the independent harmonization
switch and actual sources even when read-SV calling is disabled.

```yaml
schema_version: 1
annotation:
  - sample_id: SAMPLE_A
    track: sv
    origin: GENERATED
    producer: sv.harmonization
    output_role: bnd_vcf
    reference_build: GRCh38
    variant_ids: [EXACT_GENERATED_RECORD_ID]
```

This is a fragment, not a complete run configuration. IDs are exact future
harmonized record IDs, not guessed raw-caller IDs; absent/ambiguous selection
fails. For exploratory selection after seeing results, use a separate explicit
ready-input invocation rather than fabricating a completed producer receipt.

Only a same-attempt completed harmonization receipt may materialize this role.
The consumer independently checks the derived VCF/index, matching manifest and
every mandatory bundle member against sealed producer outputs, fixed policy,
exact reference and native-record round trip. Copied native VCF/index/evidence
must equal the unchanged canonical producer bytes. The downstream annotation
manifest carries its explicit sixth bundle column; receipt sealing/readback
repeat the binding. Missing bundles, policy/role changes, altered outputs and
cross-attempt substitutions fail, without rebuilding or silently using native.

The consumer-only configuration disables the independent harmonization switch
and its producer-side BND/source-identity policies, not merely sv.enabled. The
original policies remain frozen in the launch identity and producer receipts;
the annotation consumer still verifies the explicit bundle. Existing-bundle
no-overwrite guards remain unchanged. Otherwise a second DAG could register a
producer again and incorrectly block on the already produced bundle.

Final reporting verifies the same explicit producer/output/index/bundle selector
and preserves the derived manifest in source rows. A report must not reject a
valid bnd_vcf as an unknown old-style VCF role or silently substitute its native
index. This provenance check does not turn a failed attempt into a success.

External source-plan rows may instead use `origin: EXTERNAL` with existing
VCF/index paths, source metadata and `bnd_bundle_manifest`. Those full bundle
bytes are frozen as external inputs. Generated rows must not supply a future
bundle path: its path is learned exclusively from the selected producer receipt.

原生与派生角色显式区分；下游只消费本次已验收产物，不重新调用上游工具。不能把缺失
manifest、错误reference或旧attempt产物悄悄降级为可用。派生角色仍不证明BND真假、
基因融合、临床意义或人工复核完成。该入口的工程/原生验收与科学准确率分别记录；
旧RC1–RC4及相同版本号的其他未提交候选包不能替代当前字节身份。
