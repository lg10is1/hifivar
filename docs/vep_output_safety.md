# VEP output safety / VEP输出保护

[简体中文](zh_CN/vep_output_safety.md)

This describes a later **RC5 release source**, not historical RC1-RC4
or every package named `0.1.0rc5`. Use exact source/package identities and their
acceptance scope. No new calling/annotation algorithm or accuracy claim is added.

## Output transaction

VEP annotation runs in an exclusively owned `.vep-*` directory below the final
TSV's parent, after version detection succeeds. Its native TSV and ordinary
statistics sidecars are written there, not at an existing final-output prefix.
Only a successful command with a structurally validated native TSV can publish
the final artifact. Nonempty content alone is insufficient: an error page or
partial table must not replace a previous annotation result.

The later source format gate streams UTF-8 text before publication. It requires
a single tab-separated column header containing VEP's 13 default fields,
distinct nonempty column names, matching row widths and nonempty default-field
cells (native missing values use `-`). `##` metadata is permitted before the
header; extra named columns and their ordering are preserved. Header-only
output is legal, as are multiple feature rows per variant. The gate neither
rewrites bytes nor requires one output row per input variant. It does not prove
annotation accuracy, cache completeness, input/output identity or a clinical
interpretation. Custom narrowed `--fields` output is outside this wrapper's
current default-tab contract.

This contract follows the official [tabular output implementation](https://github.com/Ensembl/ensembl-vep/blob/release/116/modules/Bio/EnsEMBL/VEP/OutputFactory/Tab.pm)
and [default fields](https://github.com/Ensembl/ensembl-vep/blob/release/116/modules/Bio/EnsEMBL/VEP/Constants.pm).

| Event | Final TSV | Owned staging / native sidecars |
| --- | --- | --- |
| Existing output and `overwrite: false` | Rejected before probing; unchanged | None created |
| Version failure | Previous TSV unchanged | None created |
| Command failure / missing, empty or malformed output | Previous TSV unchanged | Preserved for diagnosis |
| Successful explicit overwrite | Atomic replacement after validation | Sidecars preserved |
| Successful non-overwrite | Exclusive hardlink publication | Validated native TSV/sidecars retained |
| Another result appears before publication | Rejected; other result unchanged | Candidate retained |

Non-overwrite publication requires same-filesystem hardlink support. If that
operation is unavailable, fail clearly and keep staging; never silently fall
back to a replacing rename or enable overwrite. The staging directory is on
the output filesystem, not a shared system temporary directory. No automatic
deletion of historical results or failed evidence is performed.

The wrapper checks the configured output, explicit/implicit stdout and stderr
destinations, source VCF/index, reference, cache resources and configured program
for aliases. It rechecks after version probing and native execution, including
the actual staged TSV before promotion. These are fail-closed path guards, not
an OS sandbox or proof against malicious same-UID concurrent code.

Dry-run previews the requested final destination and creates no staging. Actual
execution receipts record the real unique staged `--output_file` command and
the separately validated final artifact path. Do not rewrite command records to
pretend the native tool wrote directly to the final path.

中文：先探测版本，再在最终输出父目录下创建专有暂存目录，成功且产物校验通过后才提交
最终TSV。版本、工具、产物校验或提交失败均不得先删除旧结果。VEP统计附属文件留在
自有暂存目录中，不沿最终文件前缀的既存链接写入。未允许覆盖时，晚出现的其他结果也
不能被替换；文件系统不支持安全提交时明确失败，不擅自开启覆盖。

中文：后续源码增量在提交前逐行检查UTF-8、唯一且非空的表头、13个默认字段和行列数；
默认字段缺值须按原生格式写为`-`。允许新增列、表头前的`##`元数据、一个位点的多条
feature记录和只有合法表头的零注释结果。错误页、截断行、重复表头和控制字符不能
当作成功结果覆盖旧TSV。这里只验证文件结构，不证明注释准确率、完整性或临床意义；
不因此自动执行任何被暂缓的VEP原生注释。

## Interpreter deployment is a separate concern

The configured `annotation.vep_executable` can be an explicitly prepared tool
entry. The same entry is used for `--help` version detection and annotation.
If PATH selects a Perl ABI inconsistent with the installed VEP modules, prepare
and verify a project-local entry that uses VEP's own interpreter/modules. Do not
blindly install DBI into an unrelated Perl, modify protected environments, skip
version checks or call an existing failed annotation attempt successful.

Check that entry with a credential-free environment and bounded `--help` probe
before any separately authorized annotation run. A passing version probe proves
only that deployment probe, not cache content, real annotations or a previously
failed workflow. Tool entries and plugin/library paths must be documented for
the actual deployment; no private host path is a universal installation default.

中文：Perl/PATH问题和上述TSV保护是两件事。既有`vep_executable`可以指向项目内显式工具
入口，使版本探测与实际运行使用同一入口；正确解释器的`--help`通过不等于注释整链通过。
不因本次代码修复修改旧失败状态、自动重试VEP/PAV或冒充新的原生验收。

See [known limitations](known_limitations.md), [installation](installation.md)
and the historical [Phase 11](phases/phase11.md). Annotation remains separate
from caller quality, functional prediction and human review.
