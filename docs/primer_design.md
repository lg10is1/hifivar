# Candidate PCR / Sanger primer design

[简体中文](zh_CN/primer_design.md)

This development module generates computational candidates for an explicitly
selected sequence-resolved biallelic SNV, or simple left-anchored insertion or
deletion. It is a preparation step, not experimental validation. It does not
support symbolic insertions, INV, BND, DUP, uncertain breakpoints, complex
replacement alleles, or unphased reconstruction of nearby variants.

The native adapter requires `primer3_core`; Primer3 is not bundled. All external
calls use HiFiVar's `CommandRunner`, with a finite timeout, version probe, input
record, stdout, stderr and command metadata. Parameters are explicit and are
saved with each result. A successful candidate must bind unchanged flanks on
both REF and the selected ALT. Product sequences and lengths are reconstructed
from those alleles rather than inferred from the SVTYPE label.

## Python interface

The public source CLI accepts a strict YAML plan, checks the exact biallelic VCF
record, selected sample ALT genotype and reference contig/REF, then fetches only
the bounded template through the FASTA index:

```bash
hifivar --dry-run primers design --plan examples/primer_design/plan.yaml \
  --output-dir /work/primer_candidates_01
hifivar primers design --plan examples/primer_design/plan.yaml \
  --output-dir /work/primer_candidates_01
```

The example reference/variants are **synthetic**, not human benchmark or
experimental results. Replace them with matched local inputs for real work.
`primer_sequences.tsv` is the sequence handoff. It lists each pair's forward
and reverse oligos in 5′→3′ orientation with predicted REF/ALT product lengths
and assessment status. The completed CLI command returns its path in
`primer_sequences_file`.

`primers.tsv`, `predicted_amplicons.fasta` and `predicted_bands.tsv` summarize
the candidates and REF/ALT length classes permitted by the source GT. No band
intensity, gel resolution or successful amplification is predicted. Local
reference access is bounded; no sequence is sent to an online service.

Settings include requested count, positional diversity, product-size intervals,
primer lengths, Tm/GC bounds and finite execution timeout. A custom thermodynamic
parameter directory can be declared; it is currently recorded as a path, not
sealed as a complete database. Default parameter contents depend on deployment.
See the [Chinese guide](zh_CN/validation_preparation.md) for interpretation.

The lower-level library interface is:

```python
from hifivar.primer_design import VariantTemplate, PrimerDesignSettings
from hifivar.primer3 import Primer3Wrapper

# reference_sequence is a bounded local sequence fetched from the exact
# reference at sequence_start, not an entire genome loaded into memory.
target = VariantTemplate(
    sample="HG002", contig="chr21", position=position_1based,
    ref=selected_ref, alt=selected_alt,
    sequence_start=local_start_1based, reference_sequence=local_sequence,
    source_variant_id=source_record_id,
    source_record_sha256=selected_record_sha256,
)
settings = PrimerDesignSettings(
    requested_pairs=3, product_size_ranges=((200, 800),),
    timeout_seconds=120,
)
result = Primer3Wrapper("primer3_core").design(
    target, "results/primers/unique_variant_unit", settings,
)
```

The output directory must be fresh. A dry run creates no output and makes no
native tool probe. Failures preserve the input and native logs created before
the failure; an existing unit is never overwritten. The caller is responsible
for validating the source VCF, sample, reference build and paths before passing
an exact selected allele and fetched template to this interface.

## Interpretation

- The target is three positionally distinct candidates. Primer3 first searches
  more pairs; validation may retain zero, one or two. Shortfall is reported.
- A candidate is **not genome-wide specificity-certified** merely because
  Primer3 selected it. The core records `specificity_status=NOT_ASSESSED` and
  `sample_nearby_variants_status=NOT_ASSESSED`. A separately recorded reference
  screen cannot exclude unknown sample-specific off-targets or primer-site
  variants.
- SNVs use `assay=PCR_SANGER`: REF and ALT products have the same length.
  Ordinary gel electrophoresis cannot identify the single-base substitution.
  Allele-specific PCR/ARMS and restriction assays require their own design
  logic and are not implemented by this module.
- INS/DEL products have explicit REF and ALT bp lengths. Their separability
  on a real gel, amplification, allele dropout and genotype remain untested.
- Primer3's thermodynamic metrics are retained when returned. They do not
  establish experimental dimer absence or wet-lab success. Conditions and the
  optional native thermodynamic parameter path are recorded.
- Forward/reverse sequences are 5′→3′. Genomic binding intervals are 1-based
  closed. The reverse oligo is the reverse complement of the interval ending at
  Primer3's RIGHT coordinate. Only the selected allele is reconstructed;
  unknown or unphased neighbouring variants are not silently combined.

Native outputs per unit are `primer3.input.txt`, `primer3.output.txt`,
`primer3.stderr.log` and, only after successful parsing, `primer_design.json`.
The JSON contains identity, exact products, primer binding coordinates, Tm,
GC, reported structure scores, candidate rejections, count shortfall and
`experimental_status=NOT_RUN`. Native logs can contain local sequence and are
private unless separately approved for sharing.

The local tests use synthetic sequences and fake native output to check
coordinate, orientation, allele and failure contracts. They are not a real
Primer3 installation acceptance, a PCR experiment, or an estimate of human
variant validation success.

Primer3 command and coordinate definitions:
[official manual](https://primer3.org/manual.html).
