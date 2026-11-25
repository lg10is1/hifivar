# SV depth-response study / SV测序深度响应研究

[简体中文](zh_CN/sv_depth_study.md)

**RC5 release source.** This command is not contained in all older packages with the same version number. Verify the new candidate manifest/artifact identity.

## Scientific purpose / 科学目标

The existing Phase15 DeepVariant/hap.py depth study evaluates small variants. It is **not** SV depth evidence. The new entry runs read-SV callers on explicitly listed existing depth/replicate BAMs, derives auditable INS/DEL query views and benchmarks them against one fixed matched-reference truth protocol.

研究目标是分别描绘每个caller在实测深度下的precision、recall、F1及变异类型/长度表现，再研究来源组合和支持门槛的取舍。SV支持数不是概率置信度，不用AlphaGenome评分替代calling真值或深度实验。

## Inputs and dependency / 输入与依赖

- One BAM/BAI plus an existing measured-coverage JSON for each declared unit.
- Exact sample/RG-SM, coordinate sorting and reference/index compatibility.
- One matching truth VCF/TBI and confident BED, with source/build/version declared.
- Explicit fixed caller settings and Truvari protocol for the entire study.
- Linux native callers, samtools/bcftools; same Python runtime containing Truvari5.4.0 and compatible pysam.

```bash
python -m pip install '.[sv-benchmark]'
hifivar benchmark sv-depth --help
```

The Python optional extra does not install all native callers or reference/truth data. An executable from another Conda environment does not make its Python extensions importable in the current interpreter. Never copy compiled pysam across Python ABIs.

## Execute only a declared study / 声明式运行

Use an explicitly prepared YAML manifest; unknown/duplicate keys are rejected. See the [complete schema template](../examples/sv_depth_study/manifest.template.yaml) and the study-specific frozen manifest. The schema records reference/truth/tools/thresholds/sv/design and units. Each unit retains its nominal sampling label, sampling fraction, replicate, actual BAM path and **measured observed depth**. Fractions must agree with nominal/full-input labels; nominal labels are not substituted for observed sequencing coverage.

```bash
hifivar --dry-run benchmark sv-depth --manifest sv_depth.yaml --output-dir /your/project/sv_depth/study01
hifivar benchmark sv-depth --manifest sv_depth.yaml --output-dir /your/project/sv_depth/study01 --unit-id 5x-r1 --caller sniffles2
```

The second command executes only the explicitly selected unit/caller. It is not proof that all requested study cells completed. After review, the same declared matrix can be continued:
```bash
hifivar benchmark sv-depth --manifest sv_depth.yaml --output-dir /your/project/sv_depth/study01 --resume
```

Resume reuses only complete, unchanged, reconciled output bundles. Failed/incomplete attempts are retained and are not silently retried or overwritten. A new authorized attempt/output root is needed for changed settings or a failure recovery. One study writer owns the study lock. Busy locks do not authorize deleting locks or terminating processes.

## Fixed scientific protocol / 固定比较口径

Raw caller VCFs retain all emitted SV classes. A separate query copy is restricted by native pinned Truvari classification to INS/DEL; it does not normalize raw calls or assign confidence.

**Preserve the declared original truth.** The accepted HG002 v5.0q research protocol uses the original mixed truth, including small/star records. Do not prefilter truth merely because query views are INS/DEL. Native Truvari size/filter/matching/refine rules determine the actual scored denominator. Renaming contigs does not change genome build.

Always retain TP-base and TP-comp separately: recall uses TP-base/(TP-base+FN), precision uses TP-comp/(TP-comp+FP). Reproduce F1 from those precision/recall values. No INV/BND/CNV accuracy follows from an INS/DEL benchmark.

For this project's existing inputs, nominal5/10/15/20 correspond approximately to **7.95/15.90/23.86/31.81× observed**, not5/10/15/20× measured. The three full-input units share the same read set and are technical execution repeats. Nested subsamples share molecules and are not independent biological samples.

Fixed assembly-SV callsets must not be reused at all depths and relabelled as low-depth assembly evidence. An assembly-depth study requires assemblies generated from the corresponding sampled reads and remains a separate experiment.

## Output and interpretation / 输出与解释

The study persists per-cell caller/query/benchmark identity, failures, provenance, metrics and immutable output receipts, plus a master study report and summary TSV. A bounded selection can succeed while the whole declared study remains INCOMPLETE. A present VCF or a zero exit alone is not full end-to-end acceptance.

No minimum SV depth or sequencing-cost saving is established by implementation tests or a regional smoke gate. Predeclare accuracy criteria and evaluate actual coverage, caller settings, truth scope and replicate outcomes. If the lowest tested depth passes, it is the lowest **tested** acceptable point, not a universal minimum. Generalization requires independent samples/evaluation data.

The command is a dedicated analysis entry; it does not silently alter the existing small-variant Phase15 DAG, rerun alignment, repair BAM hard clips/RG tags or change depth-dependent caller parameters. Expensive whole-genome execution follows bounded engineering acceptance and a frozen study protocol. Native process/resource limits must be configured explicitly; caller configuration values are not an OS memory cap or distributed-scheduler guarantee.

## Fixed defaults are a baseline, not a low-depth optimum / 默认参数不是低深度最优

The current default cuteSV minimum_support is10 while the lowest existing observed depth is about7.95×; fixed support can therefore limit recall. Sniffles support is automatic by default. Preserve and label these baseline settings. A tuned support profile must be a separate frozen comparison, not a hidden depth-dependent change or an assertion that a caller intrinsically cannot work at low coverage.

## Valid zero-refinement result / 合法的无需细化结果

With Truvari 5.4.0, a successful bench --refine may report zero eligible
refinement regions and exit without writing refined summaries. The new
development contract recognizes only the narrowly verified perfect-match case:
the current command exits zero; native parameters match the request; native
summary counts and rates reconcile with all four readable, indexed assignment
VCFs; both TP sides are positive and FP/FN are zero; the native logs establish the zero-region
exit; and no partial refined or phab outputs exist.

Such a result records refinement_status: NOT_NEEDED and uses the actual
summary.json as **bench** metrics. refined_metrics stays empty. No refined
summary is manufactured or copied from the bench summary. Receipts retain the
proof and its file identities; reuse rechecks them rather than trusting a status
flag. Ordinary completed refinement continues to use its native refined summary.

该分支只处理有完整证据的 Truvari 5.4.0 完全匹配情形，不把一般的缺文件视为成功。
非零退出、部分 refined 输出、日志或参数不符、分区计数不符等仍失败；其他原生
“无需细化”情形暂不放宽。这里的完全匹配仅指该比较口径的计数，不代表基因型完全一致。
比较不同单元时保留每行实际指标来源，不能把 bench
行统一标成 refined，也不能据此增加“算法提高准确率”的结论。

This is an RC5 output-contract fix. Previously failed attempts remain
failed in their historical reports; validate a newly identified candidate in a
separate attempt. A matching package version string alone is insufficient.

## Depth × caller support comparison / 深度与caller支持数联合比较

A separate development command, `benchmark sv-depth-support`, consumes accepted
outputs from `benchmark sv-depth`. It does not call variants, resample reads,
realign BAMs or change the original depth command. The
[support-plan template](../examples/sv_depth_study/support.template.yaml) names
the original study manifest and the exact source study root for each
unit/caller. Complete a four-caller unit before requesting its support analysis.

这是单caller深度研究之后的独立后处理入口。相同深度、重复中的Sawfish、
Sniffles2、pbsv、cuteSV必须全部具备通过校验的完整运行记录；只有同名样本或
现成VCF并不足够。不能把其他深度、组装来源、失败尝试或仅补跑benchmark的
结果混入。不同来源目录必须明确声明，程序不会选择最新目录或F1最高的尝试。
跨深度比较还要求同一caller的实际工具版本一致；仅使用同一份配置不足以证明
不同运行目录使用了相同二进制版本。

```bash
hifivar --dry-run benchmark sv-depth-support --manifest support.yaml --output-dir /your/project/sv_depth/support_attempt01
hifivar benchmark sv-depth-support --manifest support.yaml --output-dir /your/project/sv_depth/support_attempt01 --unit-id 5x-r1
```

The first command checks readiness without starting native tools or writing
analysis outputs. The second is a new computational analysis of existing calls.
Use a fresh output directory; this entry has no automatic retry or in-place
resume. Its source plan can map a subset of study units, but successful selected
units never imply that the full parent depth matrix completed.

每个通过准入的unit只运行一次Jasmine。完整合并结果保留不变，再分别生成
caller支持数≥1、≥2、≥3、≥4的派生VCF和证据表。所有阈值来自同一次合并；
SUPP_VEC按固定的Sawfish→Sniffles2→pbsv→cuteSV输入顺序解释，重复caller
不能重复计票。来源身份策略必须显式开启并绑定随包发布的V1策略。

The support count is the number of distinct supporting callers, not supporting
reads, biological replicates, a probability or a calibrated confidence score.
Jasmine `max_dist` and `distance_type` control clustering; the support threshold
controls retention after clustering. Truvari `refdist` controls truth matching
and is a separate parameter. Record all three layers. This entry compares four
support thresholds under one explicitly fixed Jasmine configuration; changing
the merge configuration requires a separately identified plan/output.

每个阈值随后独立派生INS/DEL评估视图，使用母研究中同一原始mixed truth、BED、
参考和Truvari参数重新匹配。完整合并VCF和原caller VCF均不改写。派生支持集合不额外保证BND配对完整性，
本比较也不评估BND准确率。不能直接过滤旧TP/FP文件来估计新的precision或recall；阈值改变后，真值与query的匹配关系也
可能改变。原truth不能因某阈值query较少而重新过滤。

### Results and provenance / 结果与溯源

- `support_study_report.json`: declared and selected units, exact source
  receipts, source identities, merge and selection artifacts, actual native
  benchmark branch, failures and per-threshold outcomes.
- `support_study_summary.tsv`: one row per unit/support threshold, with nominal
  label, observed coverage, replicate type, TP-base, TP-comp, FP, FN, precision,
  recall, F1 and native summary identity. Failed or unavailable metrics remain
  missing, not zero.
- Protected caller VCFs and shared reference/truth inputs are checked. BAM/BAI content
  hashes are inherited from the accepted calling receipts and distinguished from
  current file-stat checks; this post-processing entry does not claim a fresh
  full-BAM/BAI hash or remeasurement of coverage.

每份阈值收据必须结合最终主报告及其完整性状态核验；单个收据或VCF存在不代表
整个矩阵通过。全局来源或输出完整性失败时，不使用其中的指标绘制成功曲线。

较高阈值可能产生空query。空筛选集合本身合法，但不保证固定的Truvari refinement
协议能够处理该集合。不能把缺少refined文件自动改写成PASS，也不能静默改为另一个
比较协议；失败必须保留原因及空指标。已证实的完全匹配NOT_NEEDED分支仍使用真实
bench指标，不能统一标记为refined。

### How to use the response curves / 如何解释响应曲线

Plot observed depth on the horizontal axis, with precision, recall and F1 in
separate panels and support thresholds encoded consistently. Show every
replicate. Nominal depth labels remain sampling labels; the three full-input
runs are technical execution repeats, not independent biological replicates.
Comparisons within each nested-sampling replicate are paired.

支持门槛越高不保证precision或F1必然提高。若研究更重视precision，应先规定可接受
的precision下限，再在满足要求的设置中比较recall；若强调综合表现，可在预先指定
的开发集比较F1，同时公开两项原始指标。不能在同一HG002评估集挑出最佳k后，再把
该点宣传为经过独立验证的推荐值。此前六来源研究的后验最佳阈值不能直接套用到本
四read-caller研究。当前入口不自动推荐最佳k，不宣称最低测序深度或降低测序成本。

The new entry requires its own software/package acceptance and later native
Jasmine/Truvari execution on four accepted caller outputs. Implementation tests,
a successful dry-run, and the existing small-variant depth curve are not results
of this SV depth-by-support study.
