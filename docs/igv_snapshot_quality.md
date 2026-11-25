# IGV screenshot integrity and review status

[简体中文](zh_CN/igv_snapshot_quality.md)

The current development increment validates all screenshots returned by
`IgvWrapper`, both in the independent review DAG and the owned launcher bridge.
Install `hifivar[review]` to provide Pillow. Planning and dry-run do not need
Pillow; native execution checks this dependency before launching IGV.

Technical acceptance checks PNG signature, chunk lengths and every chunk CRC,
required IHDR/IDAT/IEND structure, absence of trailing data, full image decoding,
dimensions, file-size/pixel limits, and expected screenshot destinations.
The default minimum image is 64 × 64 pixels; the maximum is 50 million pixels
and 100 MiB per file. These are integrity limits, not journal layout criteria.
Each accepted image records dimensions, mode, bytes and SHA256. The run report
records expected and validated screenshot counts. Invalid generated images
remain at their original paths as failure evidence and do not become completed
review evidence. Existing input and output-path protections remain active.

The batch planner optionally accepts `max_panel_height` (64–50000 pixels) and
`sleep_interval_ms` (0–60000 ms), records both, and emits the corresponding IGV
`maxPanelHeight` / `setSleepInterval` commands. These settings are omitted when
not explicitly supplied, leaving IGV's display defaults in place. See the [official batch
commands](https://igv.org/doc/desktop/UserGuide/tools/batch/). The native version
probe is bounded to 30 seconds; `IgvWrapper.run` accepts `timeout_seconds`
(positive finite seconds, default 14400) for the screenshot command. Waiting
between commands does not prove track loading has succeeded.

The exact scope is `PNG_CRC_AND_FULL_DECODE_NOT_TRACK_OR_HUMAN_QA`.
The check does not establish that a track loaded, a locus is visible, the
window has coverage, a call is correct, or a human reviewed it. An empty-looking
image may still decode; no pixel heuristic assigns variant truth. Image metadata
therefore retains `track_loading=NOT_ASSESSED`, `human_review=NOT_REVIEWED` and
`biological_accuracy=NOT_ASSESSED`. Read loading/log review, locus checks and
human judgments remain separate steps. A genuine no-coverage observation must
not be silently converted to either a false call or a broken screenshot.

## Reference, resource and locus binding

The explicit `hifivar review snapshots` entry additionally requires a native
saved session immediately before every screenshot. It checks the selected
reference, the registered BAM/CRAM and VCF resources, and the requested locus.
The session's reported coordinates may differ by at most two bases from the
requested window; the tolerance and both coordinates are recorded. This layer
is `REFERENCE_RESOURCES_LOCUS_MATCH`. It is stronger than image decoding, while
remaining separate from data coverage, genotype interpretation and human review.
Image and session SHA256 values are sealed in the bundle before subsequent
manual-decision imports. Editing or relocating a sealed bundle is not silently
accepted as the same source identity.

The batch starts with an explicit, unquoted local reference command. Native
IGV 2.19.8 retains literal quote characters in its `genome` argument; a quoted
reference can fail to load while the process still exits zero. Placing the
reference command first also avoids loading a default remote genome before
the requested reference. Reference paths containing whitespace, control
characters or quotes are therefore rejected before execution with
`IGV_REFERENCE_PATH_UNSUPPORTED`. Inputs are not moved, copied or aliased to
hide this limitation. Output directories may contain spaces: session commands
use safe relative filenames and IGV executes from the owned screenshot directory.

Native acceptance used a small, explicitly synthetic reference, real
samtools-built BAM/BAI and VCF/TBI files, and IGV 2.19.8. The SNV and deletion
views underwent image decoding, saved-session identity checks and visual
software-content inspection. These are software tests, not HG002 accuracy
estimates, human review of research loci or PCR validation.

Tests use synthetic decodable PNG images to validate the software contract.
They are not native IGV screenshots or evidence of real-data interpretation.
Neither this validation nor the test images certify historical RC1–RC4 assets.
