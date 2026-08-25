# Minimal existing-BAM example

This example is a template, not bundled scientific data. Replace every `/data`
and `/work` path with real Linux paths before validation.

Requirements:

- `GRCh38.fa` and `GRCh38.fa.fai`;
- coordinate-sorted `HG002.bam` and `HG002.bam.bai`;
- exact contig compatibility between reference and alignment;
- external tools for every enabled branch.

Copy the files into a writable project directory:

```bash
cp examples/minimal/config.yaml ./config.yaml
cp examples/minimal/samples.tsv ./samples.tsv
```

Update `samples.sheet` to the copied sheet and edit all paths. Then:

```bash
hifivar --config config.yaml config validate
hifivar --config config.yaml config dump-effective --output effective_config.yaml
WORKFLOW_ROOT="$(python -c 'from hifivar.package_resources import installed_workflow_root; print(installed_workflow_root())')"
snakemake --snakefile "$WORKFLOW_ROOT/Snakefile" --configfile effective_config.yaml --cores 1 --dry-run --printshellcmds
```

With every biological branch disabled, execution writes only the Phase 0 smoke
marker. Enable and validate one branch at a time.
