# Troubleshooting

## Command not found

Activate the environment and verify `python -m pip show hifivar`. For a wheel
install, confirm the environment's `bin` directory is on `PATH`.

## Packaged Snakefile not found

Do not assume `workflow/Snakefile` exists outside a source checkout:

```bash
python -c 'from hifivar.package_resources import installed_workflow_root; print(installed_workflow_root())'
```

## Config validation fails

Global options must precede the subcommand:

```bash
hifivar --config config.yaml config validate
```

Check YAML types, required reference/sample fields, file existence and indexes.

## Reference contig mismatch

HiFiVar intentionally does not translate `chr1` and `1`. Use mutually
compatible reference, BAM/CRAM, VCF, BED and catalog resources. Normalize
outside HiFiVar only through an explicit, recorded process.

## External command failure

Preserve command, stdout, stderr, effective config and tool version. Confirm
container binds, writable output/log parents, executable availability, resource
limits and expected output names. Do not delete successful upstream artifacts.

## DeepVariant file descriptor error

Check `ulimit -n`. A minimum of 4096 is required by the wrapper preflight and
65536 is recommended where site policy permits.

## Docker or Apptainer build unavailable

Container build may be restricted by the site. Use a prebuilt, checksum-verified
image or build on an authorized host. Do not use `sudo` unless the site owner has
explicitly approved it.

## FASTQ does not enter the small-variant DAG

This is the `0.1.0rc2` boundary: small/SV/TR Snakemake rules consume indexed
BAM/CRAM. Align FASTQ through the Phase 2 Python API or an explicitly managed
pbmm2 step, then provide the resulting indexed alignment. A unified alignment
rule/CLI is not present in this release candidate.
