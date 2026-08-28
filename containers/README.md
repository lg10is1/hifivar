# Core container definitions

`hifivar.def` builds the HiFiVar core from `python:3.12-slim-bookworm` and the
current source tree. Build it from the repository root:

```bash
apptainer build hifivar_0.1.0rc2.sif containers/hifivar.def
```

The core image contains no caller binaries, databases, references or PAV/
DeepVariant images. See `docs/containers.md` and `docs/deployment.md`.
