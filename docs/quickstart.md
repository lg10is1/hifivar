# HiFiVar 0.1.0rc5 quick start

[简体中文](zh_CN/quickstart.md)

1. Follow [installation](installation.md).
2. Copy the [minimal config](../examples/minimal/config.yaml) and [sample sheet](../examples/minimal/samples.tsv). Replace all `/data/...` and `/work/...` paths. Prepare FASTA with FAI and coordinate-sorted BAM with an index. Keep sample IDs consistent with BAM RG/SM and match contig names and lengths.
3. Enable one caller and set its executable or image and resources. The template starts with callers disabled. [Deployment](deployment.md)
4. Run the first two lines to check the plan. On a compute node run the last line using the same configuration, input and run directory.

~~~bash
hifivar --config config.yaml config validate
hifivar --config config.yaml --dry-run run --run-dir ./runs/sample01
hifivar --config config.yaml run --run-dir ./runs/sample01 --resume
~~~

## FASTQ boundary

The launcher accepts indexed BAM. [Optional raw preparation](raw_input_run.md) supports selected FASTQ and native uBAM. Raw-input dry-run returns `RAW_INPUT_PLAN_ONLY`. Inspect the complete calling DAG after producing the alignments. Main-launcher CRAM remains deferred.

Check the native VCFs, indexes and run receipts. [Outputs](outputs.md) · [Single-sample Bash](single_sample_bash.md) · [Batch](batch_multivariant_run.md) · [Slurm](slurm_multi_sample.md)
