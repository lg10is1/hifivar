# HiFiVar 0.1.0rc5 安装

[English](../installation.md)

实际分析使用 Linux 或 HPC 计算节点。建议 Python 3.12 和 Snakemake 8/9。另行部署启用分支所需程序、镜像、参考及目录。

## GitHub wheel

在新目录下载 [RC5 发行文件](https://github.com/lg10is1/hifivar/releases/tag/v0.1.0-rc5) 并核对校验值。

~~~bash
base=https://github.com/lg10is1/hifivar/releases/download/v0.1.0-rc5
wheel=hifivar-0.1.0rc5-py3-none-any.whl
for asset in "$wheel" hifivar-0.1.0rc5.tar.gz SHA256SUMS; do
  curl -fL "$base/$asset" -o "$asset" || exit 1
done
sha256sum -c SHA256SUMS || exit 1
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install "./$wheel[workflow]"
hifivar --version
hifivar config validate
hifivar doctor
~~~

预期版本为 `hifivar 0.1.0rc5`。早期同版本开发包有不同校验值。历史 0.1.0rc4 仍可从原 tag 使用。

## 固定 tag 源码

~~~bash
git clone --branch v0.1.0-rc5 --single-branch https://github.com/lg10is1/hifivar.git
cd hifivar
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install '.[workflow]'
~~~

## 检查已安装的工作流

~~~bash
WORKFLOW_ROOT="$(python -c 'from hifivar.package_resources import installed_workflow_root; print(installed_workflow_root())')"
hifivar --preset standard config dump-effective --output effective_config.yaml
snakemake --snakefile "$WORKFLOW_ROOT/Snakefile" --configfile effective_config.yaml --cores 1 --dry-run
~~~

此计划检查工作流基础结构。按[快速入门](quickstart.md)配置实际样本和工具。需要自行构建 Conda 环境或容器时查阅[容器与 Conda](containers.md)。`review`、`sv-benchmark`、`alphagenome` 是对应功能所需的可选 Python 依赖。
