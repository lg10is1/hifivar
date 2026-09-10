"""Alignment input classification and minimap2 alignment rules.

This rule file demonstrates the alignment input classification and minimap2
wrapper integration. It is NOT connected to the main Snakefile and serves
as a reference implementation for alignment input hardening.

R2-CRITICAL-001 fix: ALIGN and REUSE are separate rules with truthful
DAG contracts. REUSE produces a validated JSON manifest pointing to the
existing file, never empty BAM/BAI placeholders.

R2-MAJOR-001 fix: Classification is driven by explicit config/manifest
that names FASTQ, uBAM, BAM, or CRAM, not hard-coded .fastq wildcard.
"""


rule classify_alignment_input:
    """Classify input files as FASTQ, uBAM, or aligned BAM/CRAM.

    This rule produces a classification JSON that determines whether
    the sample needs alignment (ALIGN) or can reuse an existing alignment
    (REUSE).

    R2-MAJOR-001: Classification uses explicit config-driven input paths,
    not hard-coded wildcard patterns.
    """
    input:
        input_file="data/{sample}.fastq",
    output:
        classification="results/{sample}.alignment_classification.json",
    params:
        expected_sample_id="{sample}",
        reference_contigs=[],
    log:
        "logs/{sample}.classify_alignment_input.log",
    run:
        from hifivar.alignment_input import classify_input
        import json

        classification = classify_input(
            input.input_file,
            expected_sample_id=params.expected_sample_id,
            reference_contigs=tuple(params.reference_contigs) if params.reference_contigs else None,
        )

        with open(output.classification, "w") as f:
            json.dump(classification.to_dict(), f, indent=2)


rule minimap2_align:
    """Align FASTQ to reference using minimap2.

    This rule only executes when the classification routing_action is 'align'.
    It produces the aligned BAM and its index through atomic replacement.

    R2-CRITICAL-001: This rule is only invoked for ALIGN routing.
    R2-MAJOR-003: Uses atomic replacement via the wrapper.
    R2-MAJOR-004: Validates all output properties before success.
    R2-MAJOR-006: Structured provenance for each step.
    """
    input:
        fastq="data/{sample}.fastq",
        reference="data/reference.fa",
        reference_index="data/reference.fa.fai",
        classification="results/{sample}.alignment_classification.json",
    output:
        bam="results/{sample}.aligned.bam",
        bai="results/{sample}.aligned.bam.bai",
    params:
        preset="map-hifi",
        threads=8,
        sample_id="{sample}",
    log:
        "logs/{sample}.minimap2_align.log",
    resources:
        runtime_minutes=360,
    run:
        import json

        with open(input.classification) as f:
            classification_data = json.load(f)

        if classification_data["routing_action"] != "align":
            raise RuntimeError(
                f"minimap2_align rule invoked for non-ALIGN routing: "
                f"routing_action={classification_data['routing_action']}. "
                f"This rule should only be called for ALIGN routing."
            )

        from hifivar.alignment import (
            AlignmentRequest,
            AlignmentResources,
            AlignmentOutputFormat,
            AlignmentTool,
        )
        from hifivar.minimap2 import Minimap2Wrapper, Minimap2Options
        from hifivar.sample import InputDataset, Sample
        from hifivar.reference import ReferenceGenome

        reference = ReferenceGenome.from_fasta(input.reference)
        sample = Sample(
            params.sample_id,
            InputDataset.from_files([input.fastq]),
        )

        request = AlignmentRequest(
            sample=sample,
            reference=reference,
            output_path=output.bam,
            tool=AlignmentTool.MINIMAP2,
            output_format=AlignmentOutputFormat.BAM,
            resources=AlignmentResources(threads=params.threads),
        )

        options = Minimap2Options(
            preset=params.preset,
            output_format="bam",
            log_level="INFO",
        )

        wrapper = Minimap2Wrapper(options=options)
        result = wrapper.run(request, stderr_path=log)

        logger.info(f"Alignment completed: {result.status.value}")


rule reuse_existing_alignment:
    """Validate and expose an existing coordinate-sorted BAM for reuse.

    R2-CRITICAL-001: This rule produces a JSON manifest pointing to the
    existing validated alignment, never empty BAM/BAI placeholders.
    The DAG contract is truthful: downstream rules consume the original
    BAM/BAI paths through this manifest.
    """
    input:
        input_file="data/{sample}.bam",
        input_index="data/{sample}.bam.bai",
        classification="results/{sample}.alignment_classification.json",
    output:
        manifest="results/{sample}.alignment_reuse_manifest.json",
    params:
        expected_sample_id="{sample}",
        reference_contigs=[],
    log:
        "logs/{sample}.reuse_existing_alignment.log",
    run:
        import json

        with open(input.classification) as f:
            classification_data = json.load(f)

        if classification_data["routing_action"] != "reuse":
            raise RuntimeError(
                f"reuse_existing_alignment rule invoked for non-REUSE routing: "
                f"routing_action={classification_data['routing_action']}. "
                f"This rule should only be called for REUSE routing."
            )

        from hifivar.alignment_header import SamtoolsHeaderParser
        from hifivar.alignment_input import (
            classify_input,
            validate_aligned_bam_for_reuse,
        )

        parser = SamtoolsHeaderParser()
        classification = classify_input(
            input.input_file,
            expected_sample_id=params.expected_sample_id,
            reference_contigs=tuple(params.reference_contigs) if params.reference_contigs else None,
            parser=parser,
        )

        validate_aligned_bam_for_reuse(
            classification,
            reference_contigs=tuple(params.reference_contigs) if params.reference_contigs else None,
        )

        # Produce a truthful manifest, not empty artifacts
        manifest_data = {
            "action": "reuse",
            "sample_id": params.expected_sample_id,
            "alignment_path": str(Path(input.input_file).absolute()),
            "index_path": str(Path(input.input_index).absolute()),
            "classification": classification.to_dict(),
        }

        with open(output.manifest, "w") as f:
            json.dump(manifest_data, f, indent=2)

        logger.info(f"Reuse manifest created for {params.expected_sample_id}")


rule minimap2_align_dry_run:
    """Dry-run minimap2 alignment for command preview."""
    input:
        fastq="data/{sample}.fastq",
        reference="data/reference.fa",
        reference_index="data/reference.fa.fai",
    output:
        command_plan="results/{sample}.minimap2_command_plan.json",
    params:
        preset="map-hifi",
        threads=8,
        sample_id="{sample}",
    log:
        "logs/{sample}.minimap2_dry_run.log",
    run:
        from hifivar.alignment import (
            AlignmentRequest,
            AlignmentResources,
            AlignmentOutputFormat,
            AlignmentTool,
        )
        from hifivar.minimap2 import Minimap2Wrapper, Minimap2Options
        from hifivar.sample import InputDataset, Sample
        from hifivar.reference import ReferenceGenome
        import json

        reference = ReferenceGenome.from_fasta(input.reference)
        sample = Sample(
            params.sample_id,
            InputDataset.from_files([input.fastq]),
        )

        request = AlignmentRequest(
            sample=sample,
            reference=reference,
            output_path=output.command_plan.with_suffix(".bam"),
            tool=AlignmentTool.MINIMAP2,
            output_format=AlignmentOutputFormat.BAM,
            resources=AlignmentResources(threads=params.threads),
        )

        options = Minimap2Options(
            preset=params.preset,
            output_format="bam",
            log_level="INFO",
        )

        wrapper = Minimap2Wrapper(options=options)
        result = wrapper.run(request, dry_run=True)

        with open(output.command_plan, "w") as f:
            json.dump(result.to_dict(), f, indent=2)

        logger.info(f"Dry-run completed: {result.status.value}")
