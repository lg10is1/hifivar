from pathlib import Path

from hifivar.command import CommandResult
from hifivar.igv import IgvWrapper
from hifivar.phase10 import Phase10Status, run_phase10
from hifivar.review import read_review_selection


class FakeIgvRunner:
    def require_executable(self, executable):
        return Path(executable)

    def run(self, command, **kwargs):
        args = tuple(str(value) for value in command)
        if "--version" in args:
            return CommandResult(args, 0, "IGV 2.19.8", "", 0.1, None, True)
        batch = Path(args[-1])
        snapshot_root = None
        for line in batch.read_text(encoding="utf-8").splitlines():
            if line.startswith("snapshotDirectory "):
                snapshot_root = Path(line.split(" ", 1)[1].strip('"'))
                snapshot_root.mkdir(parents=True, exist_ok=True)
            elif line.startswith("snapshot "):
                (snapshot_root / line.split(" ", 1)[1]).write_bytes(b"PNG")
        return CommandResult(args, 0, "", "", 0.2, None, True)


def test_phase10_explicit_phase3_to_phase9_sources_end_to_end(tmp_path):
    reference = tmp_path / "reference.fa"
    reference.write_text(">chr1\n" + "A" * 1000 + "\n", encoding="utf-8")
    Path(f"{reference}.fai").write_text("chr1\t1000\t6\t1000\t1001\n", encoding="utf-8")
    bam = tmp_path / "S1.bam"
    bam.write_bytes(b"BAM")
    Path(f"{bam}.bai").write_bytes(b"BAI")
    classes = (
        ("small", "small_variant", "SNV"),
        ("sawfish", "read_sv", "DEL"),
        ("trgt", "tandem_repeat", "TR"),
        ("hiphase", "phased_variant", "INDEL"),
        ("pav", "assembly_sv", "INS"),
        ("jasmine", "harmonized_sv", "DUP"),
    )
    header = (
        "review_id\tsample\tvariant_id\tvariant_type\tcontig\tstart\tend\t"
        "source_vcf\tsource_caller\tevidence_class\tmate_contig\tmate_position\t"
        "flank_bp\ttrgt_visualization\n"
    )
    rows = []
    originals = {}
    for index, (caller, evidence, variant_type) in enumerate(classes, 1):
        vcf = tmp_path / f"S1.{caller}.vcf"
        vcf.write_text(
            "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
            f"chr1\t{100 + index}\tv{index}\tA\tT\t.\tPASS\t.\tGT\t0/1\n",
            encoding="utf-8",
        )
        originals[vcf] = vcf.read_bytes()
        rows.append(
            f"R{index}\tS1\tv{index}\t{variant_type}\tchr1\t{100 + index}\t{105 + index}\t"
            f"{vcf.name}\t{caller}\t{evidence}\t\t\t25\t\n"
        )
    selection = tmp_path / "selected_variants.tsv"
    selection.write_text(header + "".join(rows), encoding="utf-8")
    output = tmp_path / "review"
    targets = read_review_selection(
        selection,
        alignments={"S1": bam},
        reference_fasta=reference,
        output_directory=output,
    )
    report = run_phase10(
        targets,
        output_directory=output,
        igv_wrapper=IgvWrapper(runner=FakeIgvRunner()),
    )
    assert report.status is Phase10Status.COMPLETED
    assert len(report.manifest.results) == 6
    assert {result.target.evidence_class.value for result in report.manifest.results} == {
        "small_variant", "read_sv", "tandem_repeat", "phased_variant",
        "assembly_sv", "harmonized_sv",
    }
    assert all(result.status.value == "NOT_REVIEWED" for result in report.manifest.results)
    assert all(path.read_bytes() == content for path, content in originals.items())
    assert report.to_dict()["scientific_policy"]["raw_variant_artifacts_modified"] is False
