"""Two-sample generated-read-SV DAG contract, without running callers."""
from pathlib import Path
import os
import shutil
import subprocess
import sys

import pytest
import yaml

from hifivar.config import load_config
from hifivar.breakends import BndPolicy, POLICY, POLICY_SHA256
from hifivar.launcher import prepare_launch
from hifivar.serialization import write_yaml_atomic

ROOT = Path(__file__).resolve().parents[2]


@pytest.mark.parametrize('bnd_enabled', [False, True])
def test_snakemake_launcher_two_sample_harmonization(tmp_path, bnd_enabled):
    adjacent = Path(sys.executable).with_name('snakemake.exe' if sys.platform == 'win32' else 'snakemake')
    executable = str(adjacent) if adjacent.is_file() else shutil.which('snakemake')
    if not executable:
        pytest.skip('Snakemake is not installed.')
    fasta = tmp_path / 'ref.fa'
    fasta.write_text('>chr1\nACGT\n', encoding='utf8')
    Path(str(fasta) + '.fai').write_text('chr1\t4\t6\t4\t5\n', encoding='utf8')
    rows = ['sample_id\tinput\tinput_type']
    for sample in ('S2', 'S1'):
        bam = tmp_path / f'{sample}.bam'
        bam.write_bytes(b'SYNTHETIC_NOT_BAM')
        Path(str(bam) + '.bai').write_bytes(b'SYNTHETIC_NOT_INDEX')
        rows.append(f'{sample}\t{bam}\tbam')
    sheet = tmp_path / 'original_samples.tsv'
    sheet.write_text('\n'.join(rows) + '\n', encoding='utf8')
    user = tmp_path / 'analysis.yaml'
    user.write_text(yaml.safe_dump({'reference': {'fasta': str(fasta), 'build': 'CUSTOM_TEST'},
        'samples': {'sheet': str(sheet)}, 'runtime': {'threads': 4}, 'workflow': {'memory_mb': 4000},
        'sv': {'enabled': True,
            'sawfish': {'enabled': True, 'threads': 2, 'memory_mb': 2000},
            'sniffles2': {'enabled': True, 'threads': 2, 'memory_mb': 2000},
            'pbsv': {'enabled': False}, 'cutesv': {'enabled': False},
            'harmonization': {'enabled': True, 'threads': 2, 'memory_mb': 2000,
                'bnd_representation': BndPolicy(POLICY['policy_id'], 1, POLICY_SHA256, POLICY['mode']).to_dict()
                if bnd_enabled else {'enabled': False}}}}), encoding='utf8')
    cfg = load_config(ROOT / 'src/hifivar/resources/configs/default.yaml', user_config=user)
    run = tmp_path / 'owned_run'
    effective, _, derived_sheet, _ = prepare_launch(cfg, run)
    run.mkdir()
    (run / 'samples.tsv').write_text(derived_sheet, encoding='utf8')
    config_path = run / 'effective.yaml'
    write_yaml_atomic(effective, config_path)
    env = os.environ.copy()
    env.update(XDG_CACHE_HOME=str(tmp_path / 'cache'), XDG_CONFIG_HOME=str(tmp_path / 'config'))
    if sys.platform == 'win32':
        env.update(LOCALAPPDATA=str(tmp_path / 'local'), APPDATA=str(tmp_path / 'app'))
    result = subprocess.run([executable, '--snakefile', str(ROOT / 'workflow/Snakefile'),
        '--configfile', str(config_path), '--cores', '4', '--dry-run', '--printshellcmds',
        '--shared-fs-usage', 'input-output', 'persistence', 'software-deployment',
        'software-deployment-cache', 'sources', 'storage-local-copies', '--', 'all'],
        cwd=run, env=env, capture_output=True, text=True, encoding='utf8', errors='replace',
        timeout=60, check=False)
    output = result.stdout + result.stderr
    assert result.returncode == 0, output
    for sample in ('S2', 'S1'):
        assert f'{sample}.harmonized.sv.vcf.gz' in output
        assert f'{sample}.sv.evidence.tsv' in output
        assert f'{sample}.phase9.provenance.json' in output
        if bnd_enabled:
            assert (run / 'results/sv_harmonized' / sample / 'bndrepr-v1').as_posix() in output.replace('\\', '/')
        for caller in ('sawfish', 'sniffles2'):
            assert f'{sample}.{caller}.sv.vcf.gz' in output
    assert output.count('rule harmonize_sv:') == 2
    assert output.count('rule read_based_sv:') == 4
    assert 'pav_assembly_sv' not in output and 'cohort_small' not in output
    assert not (run / 'results').exists()
    assert ('bndrepr-v1' in output) is bnd_enabled
