# Contributing

HiFiVar welcomes reproducible bug reports and narrowly scoped changes.

Before opening a change:

1. describe the affected module and scientific/runtime boundary;
2. preserve raw caller outputs and existing public contracts;
3. add unit or fake-tool regression tests;
4. run `python -m pytest`, `python -m compileall src` and `python -m build`;
5. do not include WGS data, credentials, private paths, logs or proprietary
   databases;
6. mark real-tool claims as pending until independently validated on Linux.

External commands must remain in dedicated wrappers and execute through
`CommandRunner`. Do not introduce silent contig renaming, overwrite behavior or
clinical interpretation.

By intentionally submitting a contribution for inclusion, you agree that it is
provided under the Apache License 2.0 as described in `LICENSE`.
