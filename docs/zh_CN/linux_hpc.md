# Linux/HPC 运行

[English](../linux_hpc.md)

在计算节点通过站点调度器分配线程、内存和临时存储。激活 HiFiVar 环境并按分支配置真实可执行程序或镜像。预演工作流、检查预计任务和资源后再执行，记录软件、运行参数及输出日志。登录节点适合准备配置和检查状态。

## 运行命令

下面的示例使用项目路径占位符。完整配置字段和参数以英文参考为准。

~~~bash
conda activate hifivar
hifivar --version
hifivar doctor
snakemake --version
apptainer --version
ulimit -n
~~~

~~~bash
snakemake \
  --snakefile "$WORKFLOW_ROOT/Snakefile" \
  --configfile effective_config.yaml \
  --cores 1 \
  --dry-run \
  --printshellcmds
~~~

~~~bash
snakemake \
  --snakefile "$WORKFLOW_ROOT/Snakefile" \
  --configfile effective_config.yaml \
  --profile /path/to/site-profile \
  --rerun-incomplete \
  --printshellcmds
~~~

[完整参数、配置和边界](../linux_hpc.md)
