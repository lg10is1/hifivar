from pathlib import Path
import struct
import zlib

from hifivar.assembly_sv import SVEvidenceSource
from hifivar.command import CommandResult
from hifivar.context import AnalysisContext
from hifivar.harmonization import EvidenceRunStatus, SVEvidenceSourceArtifact
from hifivar.jasmine import JasmineWrapper
from hifivar.phase9 import run_phase9
from hifivar.reference import ReferenceGenome
from hifivar.sample import InputDataset, Sample
from hifivar.sample_sheet import SampleRecord
from hifivar.truvari import TruvariWrapper


def bgzf(path, payload):
    compressor = zlib.compressobj(level=6, wbits=-15)
    data = compressor.compress(payload) + compressor.flush()
    size = 18 + len(data) + 8
    header = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00" + struct.pack("<H", size - 1)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(header + data + struct.pack("<II", zlib.crc32(payload), len(payload)))


def vcf(sample, records=""):
    return (
        "##fileformat=VCFv4.3\n"
        "##contig=<ID=chr1,length=1000>\n"
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n"
        + records
    ).encode()


class FakeRunner:
    def require_executable(self, executable):
        return Path(executable)

    def run(self, command, **kwargs):
        args = tuple(str(item) for item in command)
        dry_run = kwargs.get("dry_run", False)
        if "--version" in args or args[-1:] == ("version",):
            return CommandResult(args, 0, f"{args[0]} 1.1.5", "", 0.1, None, True)
        if not dry_run and args[0] == "jasmine":
            raw = Path(next(item.split("=", 1)[1] for item in args if item.startswith("out_file=")))
            raw.write_bytes(
                vcf(
                    "S1",
                    "chr1\t10\tJ1\tA\t<DEL>\t.\tPASS\t"
                    "SVTYPE=DEL;END=20;SUPP_VEC=111111;IDLIST=r1,r2,r3,r4,r5,r6"
                    "\tGT\t0/1\n",
                )
            )
        elif not dry_run and args[0] == "bgzip":
            bgzf(Path(kwargs["stdout_path"]), Path(args[-1]).read_bytes())
        elif not dry_run and args[0] == "tabix":
            bgzf(Path(f"{args[-1]}.tbi"), b"TBI\x01")
        elif not dry_run and args[:2] == ("truvari", "bench"):
            output = Path(args[args.index("-o") + 1])
            output.mkdir(parents=True)
            (output / "summary.json").write_text('{"comparison_only": true}\n')
        return CommandResult(args, 0, "", "", 0.1, None, not dry_run)


def test_phase9_all_six_sources_end_to_end_with_provenance(tmp_path):
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">chr1\n" + "A" * 1000 + "\n")
    Path(f"{fasta}.fai").write_text("chr1\t1000\t6\t1000\t1001\n")
    bam = tmp_path / "S1.bam"
    bam.write_bytes(b"BAM")
    Path(f"{bam}.bai").write_bytes(b"BAI")
    context = AnalysisContext(
        ReferenceGenome.from_fasta(fasta, build="GRCh38"),
        (SampleRecord(Sample("S1", InputDataset.from_files([bam]))),),
        {},
    )
    callers = (
        ("sawfish", SVEvidenceSource.READ),
        ("sniffles2", SVEvidenceSource.READ),
        ("pbsv", SVEvidenceSource.READ),
        ("cutesv", SVEvidenceSource.READ),
        ("pav", SVEvidenceSource.ASSEMBLY),
        ("svim_asm", SVEvidenceSource.ASSEMBLY),
    )
    sources = []
    for number, (caller, source_type) in enumerate(callers, 1):
        path = tmp_path / f"{caller}.vcf.gz"
        bgzf(path, vcf("S1", f"chr1\t{number}\tr{number}\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=20\tGT\t0/1\n"))
        bgzf(Path(f"{path}.tbi"), b"TBI\x01")
        sources.append(
            SVEvidenceSourceArtifact(
                "S1", source_type, caller, path, Path(f"{path}.tbi"),
                EvidenceRunStatus.COMPLETED,
                ("haplotype1", "haplotype2") if source_type is SVEvidenceSource.ASSEMBLY else (),
            )
        )
    config = {"sv": {"harmonization": {
        "enabled": True, "backend": "jasmine",
        "jasmine_executable": "jasmine", "truvari_executable": "truvari",
        "bgzip_executable": "bgzip", "tabix_executable": "tabix",
        "threads": 2, "memory_mb": 1000, "runtime_minutes": 30,
        "max_dist": 1000, "distance_type": "linear", "overwrite": False,
    }}}
    runner = FakeRunner()
    report = run_phase9(
        context, {"S1": tuple(sources)}, output_directory=tmp_path / "out",
        work_directory=tmp_path / "work", config=config,
        jasmine_wrapper=JasmineWrapper(runner=runner),
        truvari_wrapper=TruvariWrapper(runner=runner),
    )
    sample = report.sample_results[0]
    assert len(sample.sources) == 6
    assert len(sample.truvari) == 6
    assert sample.jasmine.artifact is not None
    evidence = sample.jasmine.artifact.evidence_table.read_text()
    assert "READ_AND_ASSEMBLY" in evidence
    assert "sawfish,sniffles2,pbsv,cutesv,pav,svim_asm" in evidence
    assert report.to_dict()["scientific_policy"]["harmonization_is_truth"] is False
    assert report.write_json(tmp_path / "phase9.json").exists()
