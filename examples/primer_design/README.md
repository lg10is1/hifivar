# Synthetic primer-design example

[简体中文](README.zh-CN.md)

The reference and three variants are deterministic synthetic data for software
acceptance only, not a human dataset, accuracy estimate or wet-lab result.

Run the plan through `hifivar primers design --plan examples/primer_design/plan.yaml --output-dir /work/new_primer_demo`. Use a fresh directory.

Primer3 is an external dependency. Candidates remain genome-wide specificity
NOT_ASSESSED and experimental status NOT_RUN. SNV primers are for PCR/Sanger;
ordinary gels do not distinguish the equal-length products.
