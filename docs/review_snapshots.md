# Execute candidate IGV snapshots

[简体中文](zh_CN/review_snapshots.md)

This is an **unreleased development increment**. It is not a claim about the
contents of RC1–RC4 assets. Existing `review prepare` remains offline.

Install this source candidate with the review extra in a suitable environment:

```bash
python -m pip install '.[review]'
```

IGV Desktop, Java, samtools and a working desktop/display are external deployment
requirements. On a properly configured headless Linux host, `xvfb-run -a` may
wrap the command below. HiFiVar does not install or modify the display service.

Configure the existing `reference`, `samples`, `workflow.samtools_executable`
and `review.igv_executable` fields. Set `review.source_plan` to an explicit plan
as described in [review sources](review_sources.md); leave
`review.selection_file` null. This direct entry accepts existing **EXTERNAL**
sources. GENERATED sources continue to use the owned launcher DAG and its
producer receipts. No calling or automatic candidate selection is started.

```bash
hifivar --config config.yaml --dry-run review snapshots \
  --settings capture_settings.yaml --output-dir /work/review_candidate_01
hifivar --config config.yaml review snapshots \
  --settings capture_settings.yaml --output-dir /work/review_candidate_01
```

Dry-run does not reserve the directory or execute native tools. Execution requires
a fresh output directory. `capture_settings.yaml` is optional; supported defaults:

```yaml
view_mode: scientific
small_flank_bp: 40
breakpoint_flank_bp: 250
long_variant_bp: 1000
max_panel_height: 1000
sleep_interval_ms: 1000
timeout_seconds: 14400
```

Scientific views use a close interval for SNVs/INDELs/INS, an overview for other
ordinary variants, and additional left/right breakpoint images for long spans.
BND views retain the selected ALT's real endpoints. Unknown/virtual endpoints
are recorded but not rendered. Contig boundaries clip scientific views. Existing
legacy target construction defaults remain unchanged. Window/panel/delay settings
are saved; other IGV user preferences, quality thresholds and downsampling are
not fully frozen by this first increment. Their effect still needs native review.

IGV 2.19.8's batch `genome` parser retains literal quotes. This increment uses
an unquoted local reference path. The native CLI alternative did not select the
whitespace reference reliably in the tested 2.19.8 deployment. Reference paths
containing whitespace or quotes therefore fail explicitly before execution;
support for them is deferred. Whitespace output directories remain supported:
per-image sessions use safe relative names within the recorded working directory.
No input is moved, copied or substituted to bypass this limitation.

Before IGV, samtools checks the BAM header's sample/reference dictionary,
coordinate-sort declaration, EOF and one indexed-access region. This is not a
complete BAM record audit. Selected VCF records/index identities are bound to the
source plan. Large BAM/FASTA identities use stat metadata rather than full-file
hashes; smaller controls and sources receive checksums. No input is rewritten.

Outputs include a local `index.html`, screenshots, batch, source binding,
`snapshot_result.json`, capture settings, review manifests, and
`manual_decisions.tsv`. Full PNG structure/CRC/decode and dimensions are checked
on all wrapper execution paths. The direct entry additionally saves and validates
one IGV session XML per image: exact reference and selected BAM/VCF resource
paths, plus the requested locus (at most two bases of display rounding).
PNG checks alone cannot prove correct track loading,
visible read evidence, variant accuracy or human judgment. Empty-coverage views
may be legitimate observations. Images are **NOT_REVIEWED** initially.

To record actual human decisions, edit the separate TSV template with
`SUPPORT`, `NOT_SUPPORT` or `UNCERTAIN`, a reviewer, timezone-aware date and reason,
then import to another new directory:

```bash
hifivar review import-decisions --bundle-dir /work/review_candidate_01 \
  --decisions reviewed.tsv --output-dir /work/review_decisions_01
```

The snapshot bundle seals its images and execution record. Editing those files
or the per-image session XML
invalidates the bundle; the decisions template is editable. Moving a bundle does
not silently retarget its original absolute source/screenshot bindings. Human
review uses the calling BAM and is not independent experimental truth. Originals
and failed-attempt evidence are retained. For the generic offline path see
[review CLI](review_cli.md).

Official IGV commands: <https://igv.org/doc/desktop/UserGuide/tools/batch/>.
