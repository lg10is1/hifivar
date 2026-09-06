# Benchmark Data Dictionary

## File: hg002_small_variant_benchmark.tsv

| Column | Type | Description |
|--------|------|-------------|
| variant_class | string | Variant class (SNP or INDEL) |
| precision | float | Precision score (0.0-1.0) |
| recall | float | Recall score (0.0-1.0) |
| f1 | float | F1 score (0.0-1.0) |

## File: hg002_small_variant_stratified.tsv

| Column | Type | Description |
|--------|------|-------------|
| stratification | string | GIAB genome difficulty category |
| variant_class | string | Variant class (SNP or INDEL) |
| precision | float | Precision score (0.0-1.0) |
| recall | float | Recall score (0.0-1.0) |
| f1 | float | F1 score (0.0-1.0) |

### Stratification categories

| Category | Description |
|----------|-------------|
| all_difficult_regions | Union of all difficulty categories |
| not_in_all_difficult_regions | Complement of all_difficult_regions |
| low_complexity | Low-complexity regions |
| homopolymer | Homopolymer runs |
| tandem_repeats | Tandem repeat regions |
| segmental_duplications | Segmental duplication regions |
| mappability | Low-mappability regions |
| GC_content | Extreme GC content regions |
| functional_regions | Functional genomic regions |

## File: hg002_cmrg_sv_benchmark.tsv

| Column | Type | Description |
|--------|------|-------------|
| caller | string | SV caller name |
| version | string | Caller version |
| tp | integer | True positives |
| fp | integer | False positives |
| fn | integer | False negatives |
| precision | float | Precision score (0.0-1.0) |
| recall | float | Recall score (0.0-1.0) |
| f1 | float | F1 score (0.0-1.0) |
| gt_concordance | float | Genotype concordance (0.0-1.0) |
