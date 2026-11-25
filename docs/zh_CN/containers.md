# Conda、Docker 和 Apptainer

[English](../containers.md)

按固定 tag 建立 Conda/Mamba 源码环境，或在本地构建核心镜像。镜像主要包含包和工作流，调用器、参考、数据库及运行路径由部署者提供。运行前试验非 root 权限、挂载路径、版本与配置；新容器与调用器组合需独立验证。

## 运行命令

下面的示例使用项目路径占位符。完整配置字段和参数以英文参考为准。

~~~bash
mamba env create -f environment.yml
conda activate hifivar
hifivar --version
~~~

~~~bash
mamba create -n hifivar-build -c conda-forge conda-build conda-verify
conda activate hifivar-build
conda build conda-recipe
~~~

~~~bash
docker build --pull -t hifivar:0.1.0rc5 .
docker run --rm hifivar:0.1.0rc5 --version
docker run --rm hifivar:0.1.0rc5 config validate
~~~

~~~bash
docker run --rm \
  --user "$(id -u):$(id -g)" \
  -v "$PWD:/work" \
  hifivar:0.1.0rc5 \
  --config /work/config.yaml config validate
~~~

~~~bash
docker run --rm \
  --entrypoint /bin/sh \
  -v "$PWD:/work" \
  hifivar:0.1.0rc5 \
  -c 'WORKFLOW_ROOT=$(python -c "from hifivar.package_resources import installed_workflow_root; print(installed_workflow_root())"); snakemake --snakefile "$WORKFLOW_ROOT/Snakefile" --configfile /work/effective_config.yaml --cores 1 --dry-run'
~~~

~~~bash
apptainer build hifivar_0.1.0rc5.sif containers/hifivar.def
apptainer test hifivar_0.1.0rc5.sif
apptainer run hifivar_0.1.0rc5.sif --version
~~~

[完整参数、配置和边界](../containers.md)
