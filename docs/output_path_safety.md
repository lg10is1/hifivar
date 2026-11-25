# Output path safety / 输出路径保护

[简体中文](zh_CN/output_path_safety.md)

These checks are part of the unreleased source candidate. They do not modify
historical release assets or imply a new caller/benchmark result.

## TRGT

Before version probes, command execution or output cleanup, TRGT validates its
six output paths and the actual per-command stderr/stdout companion paths.
They must not alias each other, an input BAM/index, reference FASTA/FAI,
repeat catalog, or a resolved tool executable. Resolved paths, existing
hardlinks and file/directory ancestry are checked. An explicit `overwrite`
does not authorize replacement of inputs. Catalog BED coordinates must lie
within the current reference contig length; an end equal to that length is
valid. Invalid intervals are rejected, never silently clipped.

## Standalone TR plotting

The separate unreleased [TR plotting consumer](trgt_plot.md) uses only an
exclusively created fresh directory; it has no overwrite/resume/cleanup option.
Inputs, plan and tool files are checked before mkdir/version probing and after
execution. Static SVG, logs and the last-published completion receipt are bound
by hashes; publication is followed by readback, not deletion or rollback.
An interrupted call can retain a receipt or `failure.json`, so acceptance needs
normal exit 0 plus matching evidence and no `failure.json` present.
This finite single-writer contract is not an OS sandbox or a multi-file
physically atomic transaction. No foreign or failed files are removed.

## Release bundles

Bundle-item destinations are portable relative POSIX paths. Backslashes,
drive-qualified paths, traversal and nonportable Windows file components are
rejected. The `manifests`, `reports`, `configs` and `provenance` namespaces are
reserved for generated metadata; selected result items belong outside them.

Before the first report is written, the complete write set is checked for
duplicate/ancestor destinations, pointer-name collisions, links escaping the
bundle, existing file/directory conflicts and aliases to protected sources.
Large primary files still become pointers unless explicitly included. Missing
selected inputs or a missing supplied sample sheet fail before publication;
unselected absent artifacts in failed/not-run tracks remain representable.
Normal explicitly requested overwrite of an owned bundle is supported.

These are pre-execution safety checks for a single writer, not an OS sandbox
or a promise of atomic publication of the entire directory. Do not run another
writer against the same bundle or change paths while it is being assembled.
Disk or process failure during copying remains a failed operation requiring
inspection; a partial directory is not an accepted final bundle.

## Cohort SV/TR tables

The two tables are fully validated in private sibling staging before either
formal path is published. Each destination is created exclusively; an existing
or late-arriving file is never overwritten. On a captured publication error,
only unchanged files owned by this call are rolled back. Changed external bytes
and failed staging are retained with diagnostic paths. Published bytes are
checked again against their validated snapshots before success is returned.

The TR scratch-database argument supplies a parent-directory hint. Its actual
SQLite file lives in uniquely owned staging and is logged; an older diagnostic
database at the requested name is not replaced or removed. Successful cleanup
removes only unchanged owned files and empty owned staging directories.

This is not a physically atomic two-path transaction. SIGKILL, power loss and
arbitrary concurrent path replacement can still leave incomplete artifacts;
acceptance receipts remain necessary. Table schemas, source-native identifiers,
genotypes, missing states and clustering policy are unchanged. Repeated native
SV IDs can still yield repeated source-row keys; this separate interface issue
is not silently resolved by the publication-safety change.

## 中文说明

TRGT在版本探测、日志创建和输出清理前检查产物及实际日志路径。
即使启用覆盖，也不允许覆盖输入BAM/索引、参考、catalog或工具文件。
catalog的BED终点不得超过当前参考染色体长度，等于长度合法；不会自动截断。

交付包在写入首份报告前检查完整输出集合，包括大文件pointer的实际名称。
生成的报告、配置和provenance目录为保留命名空间，不能用自选结果覆盖。
路径越界、重复／祖先冲突及输入别名会明确失败；未选择的缺失产物仍可在
NOT_RUN／失败报告中记录。该保护不等于并发文件系统隔离或整个目录的事务提交，
失败后的部分目录不能当作已验收交付物。

SV/TR cohort双表先在独占staging中完成全部验证，再逐文件排他发布。
可捕获异常只回滚本次创建且身份、内容未变的文件；外部变更、失败表和SQL保留，
便于诊断。正式产物在返回成功前再次核对完整内容。scratch参数只指定临时资源的
父目录，实际唯一SQL路径写入日志，不覆盖既有诊断数据库。
这不是物理原子的双文件事务，SIGKILL／断电后的partial仍须由回执验收识别。
本次不改变GT、missing、表字段或来源key，也没有修复重复SV原ID的来源key碰撞。
