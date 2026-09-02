# HiFiVar 0.1.0rc4 中文安装指南

[English](../installation.md) | [简体中文](installation.md)

HiFiVar 的正式生产平台是 Linux/HPC。Windows 适合配置、Python/mock 测试、
打包和 Snakemake dry-run，不作为真实 external-tool 工作流的正式平台。

## 1. 核心要求

- Python 3.10、3.11 或 3.12；
- 工作流执行需要 Snakemake 8 或 9；
- PyYAML 6 或更高版本；
- 启用分支所需的外部工具、容器、参考和数据库。

外部工具按照[部署矩阵](../deployment.md)独立安装。不要使用 `latest` 代替已验证
版本，也不要认为 core wheel 包含 DeepVariant、PAV 等大型工具。

## 2. 推荐方式：GitHub Release wheel

从 [`v0.1.0-rc4`](https://github.com/lg10is1/hifivar/releases/tag/v0.1.0-rc4)
下载 wheel、sdist 和校验文件：

```bash
mkdir -p hifivar-0.1.0rc4-install
cd hifivar-0.1.0rc4-install

curl -fLO https://github.com/lg10is1/hifivar/releases/download/v0.1.0-rc4/hifivar-0.1.0rc4-py3-none-any.whl
curl -fLO https://github.com/lg10is1/hifivar/releases/download/v0.1.0-rc4/hifivar-0.1.0rc4.tar.gz
curl -fLO https://github.com/lg10is1/hifivar/releases/download/v0.1.0-rc4/SHA256SUMS
sha256sum -c SHA256SUMS
```

只有两个 artifact 都显示 `OK` 后再安装：

```bash
python3 -m venv hifivar-env
source hifivar-env/bin/activate
python -m pip install --upgrade pip
python -m pip install './hifivar-0.1.0rc4-py3-none-any.whl[workflow]'
```

如果 HPC 无法访问 GitHub，可以在可联网的 Windows 机器下载三个文件，校验后
通过 SCP/SFTP 上传；不要关闭 TLS 校验，也不要使用来源不明的镜像文件。

## 3. 固定 tag 的源码安装

```bash
git clone https://github.com/lg10is1/hifivar.git
cd hifivar
git checkout v0.1.0-rc4
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install '.[workflow]'
```

## 4. Conda/Mamba 源码环境

```bash
git clone https://github.com/lg10is1/hifivar.git
cd hifivar
mamba env create -f environment.yml
conda activate hifivar
```

本地 `conda-recipe/` 可用下列命令检查：

```bash
conda build -c conda-forge -c bioconda conda-recipe
```

这不表示项目已经发布到 Bioconda/conda-forge。目前也尚未发布到 PyPI，因此不要
直接执行 `pip install hifivar`。

## 5. Docker/Apptainer 边界

仓库中的 core container 只包含 HiFiVar 核心和工作流，不是所有 caller 的
单体 WGS 镜像。DeepVariant、PAV 等工具使用独立的、固定版本的部署。详细说明见
[containers](../containers.md) 和[deployment](../deployment.md)。

## 6. 安装后验证

```bash
hifivar --version
hifivar --help
hifivar config validate
hifivar doctor
python -c 'from hifivar.package_resources import installed_workflow_root; print(installed_workflow_root())'
```

预期版本为 `hifivar 0.1.0rc4`。最后一条命令输出的目录应包含
`Snakefile`、`rules/`、`scripts/` 和 `envs/`。

## 7. 下一步

继续阅读[完整中文快速入门](quickstart.md)。先完成 core smoke 和 Snakemake
dry-run，再部署并逐个启用真实工具分支。
