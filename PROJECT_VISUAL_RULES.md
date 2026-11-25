# HiFiVar Scientific Figure Rules

[简体中文](PROJECT_VISUAL_RULES.zh-CN.md)

These rules apply to every manuscript, benchmark, architecture, workflow, and
README figure produced for HiFiVar. They complement the frozen-data and
provenance requirements in `docs/reproducibility.md`.

## Nature-family production target

- Prepare figures at their intended print width: 89 mm for a single column or
  183 mm for a double column; keep the composed figure at or below 170 mm high.
- Use editable vector output (`PDF` and `SVG`) for plots, diagrams, text, lines,
  and arrows. PNG is a review/README derivative, not the manuscript master.
- Use RGB colour. Raster derivatives must be exported at 300 dpi or higher at
  final size; never upscale a low-resolution source to claim compliance.
- Use Arial or Helvetica where available, with a metrically safe sans-serif
  fallback. At final print size, ordinary text must remain between 5 and 7 pt.
  Multi-panel labels use 8 pt bold upright lower-case letters.
- Keep all text editable and embed TrueType fonts in PDF. Do not convert text
  to paths unless a publisher explicitly requires it.
- Use visible axes and tick marks, label quantities and units, and keep the
  thinnest line at least 0.5 pt at final size.
- Do not use background gridlines, drop shadows, 3D effects, gradients,
  decorative DNA helices, excessive boxes, or coloured body text.
- Do not truncate bar-chart axes in a way that exaggerates differences. A
  restricted axis is permitted for point/interval plots only when its limits
  are explicit and scientifically justified.

## Colour and accessibility

- Use the shared `project_visuals/visual_style.yaml`; do not introduce local
  ad-hoc palettes in individual scripts.
- The primary palette is colour-blind accessible and remains distinguishable
  when printed in grayscale. Never encode a critical comparison by red versus
  green alone.
- Use colour together with marker shape, line style, direct labels, or spatial
  grouping whenever identity must survive monochrome reproduction.
- Reserve vermillion for errors or false positives, amber for warnings or
  uncertainty, neutral gray for unavailable/disabled data, and black for text.

## Scientific integrity

- Quantitative figures must be generated from frozen, checksummed source data.
  Generative illustration tools must never invent or modify measurements.
- Every plotted metric needs a provenance record, denominator, benchmark
  scope, software/database version, and a caption that defines statistics and
  error bars. If no interval was calculated, state that the values are point
  estimates and do not draw error bars.
- Do not describe CMRG results as whole-genome SV performance. Do not use the
  non-canonical combined small-variant score as headline performance.
- Manual review, caller support, functional impact, and call confidence remain
  distinct concepts in both labels and visual encodings.
- Any unsupported claim or planned panel must be labelled `EVIDENCE REQUIRED`
  or `PLANNED`; it must not appear as a completed result.

## Pre-submission verification

1. Regenerate the figure from frozen data in a clean environment.
2. Verify data checksums and compare plotted values with the source table.
3. Inspect PDF/SVG at 100% final print size and confirm 5–7 pt readability.
4. Run a colour-vision-deficiency and grayscale check.
5. Confirm axes, units, symbols, abbreviations, sample sizes, and statistics.
6. Confirm PDF font embedding and editable SVG text.
7. Confirm captions match the visible data and do not overstate scope.
8. Retain plotting code, input tables, environment versions, and command line.

The target journal's current author instructions always take precedence at the
time of submission.
