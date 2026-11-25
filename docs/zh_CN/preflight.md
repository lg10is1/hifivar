# 批量运行前检查：先知道缺什么，再开始计算

[English](../preflight.md)。本功能是 RC5 候选增量，不在 RC1–RC4、早期 RC5
或上一轮 review-only 候选包中。它已包含于[指定的 10 月 15 日开发包](../rc5_candidate.md#validated-development-payload--已验收开发包)，
并完成相应受限 Linux 验收。请核对源码清单和包 SHA256，不能只看版本号；文档修订后的
新包也须保留独立身份。

## 默认命令保持不变；新检查明确开启

`hifivar doctor` 仍只报告基础环境，不代表所有 caller 已安装或可以运行。
新增模式针对当前“一键运行”的新任务范围：

```bash
hifivar --config analysis.yaml doctor --preflight --run-dir /your/project/run01 --json
```

目录必须尚不存在；命令不会创建它。只读检查启用分支、配置资源预算、参考
FASTA/FAI 与明确的 build、轻量输入、所需 BAM 唯一索引、工具路径、workflow
资源，以及 TR catalog 的字段、contig、坐标边界、声明 build 和 karyotype。
没有启用的 caller 不强制安装。启用了 alignment，就检查所选 alignment 工具，
即使部分样本后续可能直接复用 BAM。
显式启用同样本 read-SV 整合时，另检查 Jasmine、bgzip/tabix 和既有 Phase 9
来源比较所需的 Truvari；整合关闭时不强制安装或探测这些依赖。
开启探测后，此 launcher 路由要求 tabix 至少 1.23.1，以支持 Jasmine 的 VCFv4.4；
旧版本在分析前明确 BLOCKED，不自动升级，也不降级 VCF header。这是该路由的
保守部署下限，不改变独立工作流处理较早 VCF 格式的原有行为。

不会安装、下载、修复输入、生成索引、比对、calling、启动 Snakemake 或重试。
这个模式也不打开配置中的日志文件；需要保存时，请显式把标准输出重定向到
你允许的项目位置。路径按已有 launcher 的规则解析，不做 chr1/1 或 SM 改写。

## 工具版本和 BAM 身份需要显式探测

```bash
hifivar --config analysis.yaml doctor --preflight --run-dir /your/project/run01 --probe-tools --timeout-seconds 30 --json
```

只允许有超时限制的 native 版本命令和 BAM 元数据查询。调用命令由既有 wrapper
负责，统一通过 CommandRunner 执行，不使用 shell。Snakemake 版本由当前 Python
的隔离模式检查；BAM 复用既有检查器核对 EOF、coordinate 声明、所有 RG/SM、
精确参考名字和长度，以及一个区域的索引访问，不遍历整个 BAM 或认证 read 质量。

HiFiVar 不写 probe 日志或生成分析目录；外部 native 程序可能管理自己的缓存，
这不是操作系统沙箱。请在可信环境中使用项目内已准备好的独立临时/缓存位置。
探测不会继承凭据或 Python 注入变量，不会代替用户确认凭据已撤销/轮换。

全局参数如 --config 在 doctor 前；doctor 参数在后。--dry-run 与 --probe-tools
不能组合，避免“dry-run”实际执行工具。额外参数必须与 --preflight 配套。

## 如何看结果

| 结果 | 退出码 | 含义 |
| --- | --- | --- |
| CHECKED_SCOPE_PASS | 0 | 本次明确范围的必需检查通过 |
| BLOCKED | 1 | 存在失败项，需明确修复配置、输入或工具环境 |
| INCOMPLETE | 2 | 没有失败，但有必需项尚未检查 |

默认仅元数据检查通常返回 2，因为版本和 BAM 身份尚未探测；这不是分析任务失败。
首个不合法 launcher 契约会阻止更深探测，不把空 phase0 任务当作批量分析完成。
JSON 保留每项状态、工具路径、版本、命令、退出码和耗时，不生成准确率指标。

资源部分只报告可用磁盘字节、可读的 inode/CPU/进程与文件描述符限制。配额、
同 UID 的真实剩余线程空间、分配内存和 WGS 所需磁盘空间仍明确为 NOT_CHECKED
或 NOT_ESTIMATED。不会自动减少线程/shard、调整科学参数或保证大型任务一定成功。

## 不能冒充完成的部分

- Apptainer image 只做可读路径元数据检查，不确认 digest、模型或容器实际部署。
  不启动/拉取任何 Apptainer/Docker 容器；即使开启探测，容器版本仍为 NOT_CHECKED。
- 原始 uBAM 只检查 header 是否适合明确的 pbmm2 路由，不做昂贵 flag/RG 全记录
  扫描，因此仍为 INCOMPLETE；实际 raw launcher 原有完整契约保持不变。
- FASTQ 只检查前缀，不认证所有记录、molecule、覆盖或最低深度。
- 本模式不支持 CRAM、resume 或 Slurm 验证；benchmark、review、low_depth
  仍受 launcher 范围限制。后续源码支持现成 SV 及部分 assembly/cohort/annotation
  路由，但受限 probe 注册表不等于这些后端的全部部署验收（见下节）。整合仍按样本
  独立执行，不引入跨样本 SV 聚类。
- scoped pass 不是 calling 成功、全流程就绪、原始数据逐字节不可变、准确率、临床
  可靠性或人工复核已完成的证明。已有受限 Linux 验收不覆盖未检查或未接入的范围。

报告可能包含样本名和绝对路径，默认作为私有运行证据保存，不自动上传。不要把
凭据写入路径/自由文本。旧包、旧报告和旧版本验证不因新入口而被重新解释。

## 后续开发源码修正

整合依赖按 `sv.harmonization.enabled` 独立检查，不再受 `sv.enabled`
控制。即使 read caller 全部关闭，现成来源仍须具备 Jasmine、Truvari 和配置的
BGZF/TBI 工具。只作样本表上下文、不进行比对的 FASTQ 不强制要求无关 BAM
检查工具；实际比对准备或 BAM 元数据探测仍检查相应依赖。显式来源身份视图
启用时，bgzip 与其 wrapper 共用 >=1.23.1 的兼容性下限。不开启版本探测仍是
INCOMPLETE，不以文件存在或可执行路径代替实际版本证据。

后续开发源码新增以下原生版本预检；命令仍由各自既有 wrapper 构建，不执行分析。

| 已启用分支 | 仅版本探测范围 |
| --- | --- |
| assembly / hifiasm | hifiasm |
| assembly_sv / SVIM-asm | SVIM-asm、minimap2、samtools、bgzip、tabix |
| annotation / VEP | VEP 程序，不是全部 cache 内容 |
| cohort / small variants | GLnexus、bcftools，不做联合分型 |

父分支和子后端开关必须同时开启。缺少依赖时阻断该组，某组版本失败不会掩盖
其他组的检查结果。不开启 --probe-tools 时，文件存在不替代版本证据。
实际成功探测保留 version_source=command 与仅版本验证范围，不认证数据库和结果。

PAV 和 ANNOVAR 仍保留必需的 NOT_CHECKED：PAV 配置版本不同于 Snakemake
adapter 版本，ANNOVAR 声明的源码日期也不是原生版本命令实测。不得借声明或
adapter 探测升级为原生 PASS。不启动容器、不下载 cache、不组装、不注释、不做
联合分型。版本探测前仍执行真实 launcher 的配置／输入／资源契约校验。
scoped PASS 仍不等于全流程就绪或科学准确率通过。

这些增量属于明确标识的未提交开发源码，不在 RC1-RC4 中，也不代表所有标为 RC5
的包都包含它。请核对新源码 manifest 和包哈希；旧包／旧原生验收不追认新载荷。
