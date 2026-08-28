from __future__ import annotations

import gzip
import json
from pathlib import Path

from hifivar.cohort import (
    CohortDefinition, CohortManifest, CohortSampleInput, CohortTrack,
    SampleCallState,
)
from hifivar.cohort_tracks import build_sv_cohort_tables, build_tr_cohort_tables
from hifivar.command import CommandResult
from hifivar.glnexus import GLnexusRequest, GLnexusResources, GLnexusWrapper
from hifivar.reference import Contig, ReferenceGenome


def write_vcf(path: Path, sample: str, record: str, info_headers: str = "") -> None:
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(f"##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n{info_headers}#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n{record}")
    Path(f"{path}.tbi").write_bytes(b"index")


class Runner:
    def __init__(self, samples): self.samples = samples
    def require_executable(self, executable, **kwargs): return Path(str(executable))
    def run(self, command, **kwargs):
        args = tuple(str(item) for item in command)
        stdout = "glnexus_cli release v1.4.1" if "--help" in args else "bcftools 1.21" if "--version" in args else ""
        if args[0] == "glnexus_cli" and "--help" not in args:
            Path(kwargs["stdout_path"]).parent.mkdir(parents=True, exist_ok=True); Path(kwargs["stdout_path"]).write_bytes(b"BCF")
        elif len(args) > 1 and args[1] == "view":
            with gzip.open(Path(args[args.index("-o") + 1]), "wt", encoding="utf-8") as handle:
                handle.write("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(self.samples) + "\nchr1\t10\tv1\tA\tC\t.\tPASS\t.\tGT\t0/1\t0/0\n")
        elif len(args) > 1 and args[1] == "index":
            target = Path(args[-1]); Path(f"{target}.tbi" if "--tbi" in args else f"{target}.csi").write_bytes(b"index")
        return CommandResult(args, 0, stdout, "", 0.1, None, True, kwargs.get("stdout_path"), kwargs.get("stderr_path"))


def test_phase12_three_tracks_and_manifest_remain_independent(tmp_path: Path) -> None:
    samples = ("S1", "S2")
    reference = ReferenceGenome(tmp_path / "ref.fa", tmp_path / "ref.fa.fai", "GRCh38", (Contig("chr1", 1000),), "b" * 64)
    cohort = CohortDefinition("family", samples, reference)
    small_inputs = []
    for sample in samples:
        path = tmp_path / f"{sample}.g.vcf.gz"
        write_vcf(path, sample, "")
        small_inputs.append(CohortSampleInput(sample, SampleCallState.NO_CALLS, path, Path(f"{path}.tbi"), "deepvariant", "1.10.0", "GRCh38"))
    request = GLnexusRequest(cohort, tuple(small_inputs), tmp_path / "glnexus.DB", tmp_path / "small" / "family.small.bcf", tmp_path / "small" / "family.small.vcf.gz", resources=GLnexusResources(2, 4))
    small = GLnexusWrapper(runner=Runner(tuple(reversed(samples)))).run(request).as_track_result()

    sv_path = tmp_path / "S1.sv.vcf.gz"
    write_vcf(sv_path, "S1", "chr1\t20\tsv1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=30\tGT\t0/1\n")
    sv = build_sv_cohort_tables(cohort, (
        CohortSampleInput("S1", SampleCallState.CALLED, sv_path, Path(f"{sv_path}.tbi"), "jasmine", "1.1.5", "GRCh38"),
        CohortSampleInput("S2", SampleCallState.FAILED),
    ), site_table=tmp_path / "sv" / "sites.tsv", sample_matrix=tmp_path / "sv" / "matrix.tsv")

    tr_path = tmp_path / "S1.tr.vcf.gz"
    write_vcf(tr_path, "S1", "chr1\t40\t.\tA\t<STR>\t.\tPASS\tTRID=L1;MOTIFS=CAG;END=45\tGT\t0/1\n", "##INFO=<ID=TRID,Number=1,Type=String,Description=x>\n##INFO=<ID=MOTIFS,Number=1,Type=String,Description=x>\n")
    tr = build_tr_cohort_tables(cohort, (
        CohortSampleInput("S1", SampleCallState.CALLED, tr_path, Path(f"{tr_path}.tbi"), "trgt", "5.1.0", "GRCh38", "catalog-1"),
        CohortSampleInput("S2", SampleCallState.NOT_RUN),
    ), locus_table=tmp_path / "tr" / "loci.tsv", sample_matrix=tmp_path / "tr" / "matrix.tsv", scratch_database=tmp_path / "tr.sqlite")

    manifest = CohortManifest(cohort, (small, sv, tr), "0.0.1.dev0")
    manifest.write(tmp_path / "cohort.json", tmp_path / "cohort.yaml")
    payload = json.loads((tmp_path / "cohort.json").read_text(encoding="utf-8"))
    assert [track["track"] for track in payload["tracks"]] == ["small_variants", "sv", "tr"]
    assert payload["tracks"][1]["sample_states"][1]["state"] == "FAILED"
    assert payload["tracks"][2]["sample_states"][1]["state"] == "NOT_RUN"
    assert payload["tracks"][0]["metrics"]["declared_sample_order"] == ["S1", "S2"]
    assert payload["tracks"][0]["metrics"]["output_sample_order"] == ["S2", "S1"]
    assert payload["tracks"][0]["metrics"]["sample_set_match"] is True
    assert payload["tracks"][0]["metrics"]["sample_order_match"] is False
    assert payload["tracks"][0]["metrics"]["per_sample_non_ref_count"] == {"S1": 0, "S2": 1}
    assert request.output_vcf.exists() and sv_path.exists() and tr_path.exists()
