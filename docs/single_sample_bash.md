# One-sample Linux Bash workflow

[简体中文](zh_CN/single_sample_bash.md)

[简体中文详细指南](zh_CN/single_sample_bash.md)

On a compute node, prepare a coordinate-sorted/indexed BAM and matching reference. For FASTQ or native uBAM, use [raw input preparation](raw_input_run.md) or align separately with pbmm2. Verify the BAM reference lengths, read-group sample name and index. Enable a provisioned caller and allocate resources for a single sample. Review the plan before executing.

~~~bash
hifivar --config config.yaml config validate
hifivar --config config.yaml --dry-run run --run-dir ./runs/sample01
hifivar --config config.yaml run --run-dir ./runs/sample01 --resume
~~~

Check the output VCF and its index and record the run receipt. The linked Chinese guide includes a separate alignment and DeepVariant example.
