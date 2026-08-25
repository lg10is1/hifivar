"""Combine independent benchmark result JSON files without altering them."""
import csv
import json
from pathlib import Path
from hifivar.serialization import write_json_atomic

results = [json.loads(Path(str(path)).read_text(encoding="utf-8")) for path in snakemake.input]  # type: ignore[name-defined]
payload = {"benchmark_id": snakemake.config["benchmark"]["benchmark_id"], "reference_build": snakemake.config["reference"]["build"], "results": results, "scientific_semantics": {"performance_is_not_pathogenicity": True, "raw_queries_modified": False}}  # type: ignore[name-defined]
write_json_atomic(payload, Path(str(snakemake.output.json)), overwrite=bool(snakemake.config["benchmark"]["overwrite"]), artifact_name="Benchmark manifest")  # type: ignore[name-defined]
path = Path(str(snakemake.output.tsv)); path.parent.mkdir(parents=True, exist_ok=True)  # type: ignore[name-defined]
with path.open("w" if snakemake.config["benchmark"]["overwrite"] else "x", encoding="utf-8", newline="") as handle:  # type: ignore[name-defined]
    writer=csv.writer(handle,delimiter="\t",lineterminator="\n"); writer.writerow(("benchmark_id","sample","variant_class","status","tool","stratum","metric","value"))
    for result in results:
        for metric in result.get("metrics",[]): writer.writerow((payload["benchmark_id"],result.get("sample_id"),result.get("variant_class"),result.get("status"),result.get("tool"),metric.get("stratum"),metric.get("name"),metric.get("value")))
