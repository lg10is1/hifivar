# HiFiVar

[简体中文](README.zh-CN.md) · [Documentation](docs/index.md)

HiFiVar coordinates PacBio HiFi small-variant, structural-variant (SV) and tandem-repeat analysis with separate source and run records. **0.1.0rc5** is a research release candidate for Linux and HPC. Provision external callers, indexes and reference resources for each enabled branch. [GitHub release](https://github.com/lg10is1/hifivar/releases/tag/v0.1.0-rc5)

1. Install the versioned wheel. [Installation](docs/installation.md)
2. Begin with one small sample and one configured caller. [Quick start](docs/quickstart.md)
3. For multiple samples see [batch execution](docs/batch_multivariant_run.md) and [Slurm](docs/slurm_multi_sample.md).

The [`hifivar run`](docs/indexed_bam_run.md) launcher takes indexed BAM. [Raw-input preparation](docs/raw_input_run.md) aligns supported FASTQ or native uBAM before calling when configured. CRAM remains deferred at the main launcher. Supported independent Snakemake branches can consume an existing indexed BAM/CRAM. Historical 0.1.0rc4 remains available at its release tag. [Capabilities and evidence](docs/capabilities.md)

Native calls, derived artifacts, and manual decisions have separate records. Check each new external-tool and reference combination in a small real run. [Outputs](docs/outputs.md) · [Benchmarking](docs/benchmarking.md) · [Known limitations](docs/known_limitations.md)

[Apache-2.0](LICENSE). External tools and resources retain their own licenses.
