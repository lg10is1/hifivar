# Phase 0: Foundation

Phase 0 establishes the project skeleton and shared engineering infrastructure before any bioinformatics tool integration.

**Phase 0 status: COMPLETE**

## Phase 0.1–0.3

- Repository skeleton
- Python packaging metadata
- Base exception hierarchy

## Phase 0.4

- Shared console and UTF-8 file logging for the `hifivar` namespace
- Strict logging-level validation and repeatable configuration

## Phase 0.5

- Safe UTF-8 YAML loading with default, preset, and user layers
- Deep merge, lightweight schema validation, provenance, and effective YAML

## Phase 0.6

- Safe list-based external command execution with explicit results and errors
- Dry-run, timeout, environment overrides, redaction, and file redirection

## Phase 0.7

- Argparse console and `python -m hifivar` entry points
- Packaged configuration resources, config inspection, and foundation doctor

## Phase 0.8

- Non-mutating path, lightweight biological text, and index-presence validation
- Contig subset compatibility, output checks, and streaming SHA256 checksums

## Phase 0.9

- Modular Snakemake entry point with common and cross-platform smoke rules
- Effective-config handoff, work/output/log paths, and resource conventions

## Phase 0.10

- End-to-end package, CLI, config, logging, validation, CommandRunner, and
  Snakemake integration smoke
- Full regression, coverage baseline, sdist/wheel build, package-resource, and
  isolated wheel-install verification

## Definition of Done

- Repository skeleton and Python packaging
- Shared exception hierarchy and logging
- Layered configuration and effective YAML
- Safe CommandRunner
- Installed CLI and packaged configuration resources
- Lightweight biological input validation
- Modular Snakemake infrastructure smoke DAG
- Unicode and non-repository working-directory support
- Full Phase 0 end-to-end and regression tests

## Not yet implemented

Deep BAM/CRAM and BGZF/index integrity validation, workflow execution,
external-tool wrappers, biological rules, and variant analysis remain outside
Phase 0. The smoke DAG is infrastructure validation, not analysis.
