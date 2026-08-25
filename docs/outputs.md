# Output layout

Configured `paths.outdir` contains retained results; `paths.workdir` contains
rebuildable intermediates. Logs remain under `logs/` unless a module documents
another path.

Representative retained artifacts include:

```text
results/
├── phase0/
├── small/<sample>.small.vcf.gz
├── small/<sample>.g.vcf.gz
├── sv/<sample>.<caller>.sv.vcf.gz
├── tr/<sample>.tr.vcf.gz
├── phasing/<sample>.phased.vcf.gz
├── assembly/
├── assembly_sv/
├── harmonization/
├── review/
├── annotation/
├── cohort/
├── benchmark/
└── reports/
```

VCF families remain separate. Indexes, manifests, tool versions, commands and
provenance accompany the corresponding track. Raw caller outputs are immutable;
annotation, benchmark and review records never overwrite them.

Large FASTQ/BAM/CRAM/reference artifacts are pointers in a release bundle by
default. A `COMPLETE` report status means the configured track completed its
contract; it does not mean a variant is pathogenic, clinically valid or true.
