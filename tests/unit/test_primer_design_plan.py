"""Synthetic plan integration; fake candidates do not certify native PCR design."""
from pathlib import Path
import csv
import json
import random

import pytest
import yaml

from hifivar.exceptions import InputValidationError, OutputValidationError
from hifivar.primer_design import evaluate_candidates
from hifivar.primer_design_plan import fetch_sequence, run_primer_design_plan
from test_primer_design import response


def fixture_plan(tmp_path):
    rng = random.Random(83240)
    sequence = ''.join(rng.choice('ACGT') for _ in range(900))
    ref = sequence[300]
    alt = next(x for x in 'ACGT' if x != ref)
    fasta = tmp_path / 'reference.fa'
    fasta.write_bytes(b'>chr21\n' + b''.join(sequence[i:i+60].encode() + b'\n' for i in range(0,900,60)))
    fai = tmp_path / 'reference.fa.fai'
    fai.write_text('chr21\t900\t7\t60\t61\n',encoding='utf-8')
    vcf = tmp_path / 'selected.vcf'
    vcf.write_text('##fileformat=VCFv4.2\n##contig=<ID=chr21,length=900>\n'
                   '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                   '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n'
                   f'chr21\t301\tSYNTHETIC_SNV\t{ref}\t{alt}\t50\tPASS\t.\tGT\t0/1\n',encoding='utf-8')
    data = {'schema_version':1,'reference':{'fasta':'reference.fa','build':'SYNTHETIC_TEST'},
            'targets':[{'target_id':'site1','vcf':'selected.vcf','sample':'HG002','contig':'chr21',
                        'position':301,'ref':ref,'alt':alt}]}
    plan = tmp_path / 'plan.yaml'
    save_plan(plan,data)
    return plan,data,fasta,fai,vcf,sequence


def save_plan(path,data):
    path.write_text(yaml.safe_dump(data,sort_keys=False),encoding='utf-8')


def fake_design(self,item,output_dir,settings):
    output_dir.mkdir()
    (output_dir/'fake-native-evidence.txt').write_text('Synthetic contract only',encoding='utf-8')
    return evaluate_candidates(item,settings,response(item))


def test_dry_run_validates_input_without_native_or_output(tmp_path,monkeypatch):
    plan,data,fasta,fai,vcf,sequence=fixture_plan(tmp_path)
    def forbidden(*args,**kwargs):raise AssertionError('dry run called native design')
    monkeypatch.setattr('hifivar.primer_design_plan.Primer3Wrapper.design',forbidden)
    before={p:p.read_bytes() for p in (plan,fasta,fai,vcf)}
    output=tmp_path/'dry'
    result=run_primer_design_plan(plan,output,dry_run=True)
    assert result['status']=='PLANNED' and result['candidate_pair_count'] is None
    assert 'primer_sequences_file' not in result
    assert not output.exists() and all(p.read_bytes()==b for p,b in before.items())


def test_plan_aggregation_retains_identity_sequences_and_unassessed_status(tmp_path,monkeypatch):
    plan,data,fasta,fai,vcf,sequence=fixture_plan(tmp_path)
    monkeypatch.setattr('hifivar.primer_design_plan.Primer3Wrapper.design',fake_design)
    output=tmp_path/'candidates'
    result=run_primer_design_plan(plan,output)
    assert result['status']=='CANDIDATE_DESIGN_COMPLETE'
    assert result['candidate_pair_count']==3 and result['targets_with_requested_count']==1
    assert result['experimental_status']=='NOT_RUN' and result['specificity_status']=='NOT_ASSESSED'
    assert result['manual_review_status']=='NOT_REVIEWED'
    assert result['targets'][0]['source_record']['record_id']=='SYNTHETIC_SNV'
    assert result['targets'][0]['source_record']['genotype']=='0/1'
    assert len(result['targets'][0]['source_record']['record_sha256'])==64
    saved=json.loads((output/'primer_design_summary.json').read_text(encoding='utf-8'))
    assert saved==result
    with (output/'primers.tsv').open(encoding='utf-8',newline='') as handle:
        rows=list(csv.DictReader(handle,delimiter='\t'))
    assert len(rows)==3 and {x['assay'] for x in rows}=={'PCR_SANGER'}
    assert {x['specificity_status'] for x in rows}=={'NOT_ASSESSED'}
    sequence_path = Path(result['primer_sequences_file'])
    assert sequence_path == output / 'primer_sequences.tsv'
    with sequence_path.open(encoding='utf-8', newline='') as handle:
        sequences = list(csv.DictReader(handle, delimiter='\t'))
    assert len(sequences) == 3
    assert 'forward_start_1based' not in sequences[0]
    for seq, full in zip(sequences, rows, strict=True):
        assert seq['forward_sequence_5to3'] == full['forward_sequence_5to3']
        assert seq['reverse_sequence_5to3'] == full['reverse_sequence_5to3']
        assert seq['candidate'] == full['candidate']
        assert seq['specificity_status'] == 'NOT_ASSESSED'
    products=(output/'predicted_amplicons.fasta').read_text(encoding='utf-8').splitlines()
    assert len(products)==12 and products[0]=='>site1.pair1.reference'
    assert len(products[1])==len(products[3])==501 and products[1]!=products[3]


def test_primer_cli_reports_copyable_sequence_table(tmp_path, monkeypatch, capsys):
    from hifivar.cli import main

    plan, *_ = fixture_plan(tmp_path)
    monkeypatch.setattr('hifivar.primer_design_plan.Primer3Wrapper.design', fake_design)
    output = tmp_path / 'cli_output'
    assert main(['primers', 'design', '--plan', str(plan), '--output-dir', str(output)]) == 0
    displayed = json.loads(capsys.readouterr().out.splitlines()[-1])
    sequence_path = output / 'primer_sequences.tsv'
    assert Path(displayed['primer_sequences_file']) == sequence_path
    assert sequence_path.is_file()


@pytest.mark.parametrize('which',['sample','position','ref','alt','contig'])
def test_exact_selected_vcf_identity_is_required(tmp_path,which):
    plan,data,*_=fixture_plan(tmp_path)
    target=data['targets'][0]
    target[which]={'sample':'OTHER','position':302,'ref':'N','alt':'<INS>','contig':'21'}[which]
    save_plan(plan,data)
    output=tmp_path/'rejected'
    with pytest.raises(InputValidationError):run_primer_design_plan(plan,output,dry_run=True)
    assert not output.exists()


@pytest.mark.parametrize('gt',['0/0','./.','0/.','0/2','1/2','x'])
def test_sample_requires_called_selected_alt_gt(tmp_path,gt):
    plan,data,fasta,fai,vcf,sequence=fixture_plan(tmp_path)
    vcf.write_text(vcf.read_text(encoding='utf-8').replace('\t0/1\n',f'\t{gt}\n'),encoding='utf-8')
    with pytest.raises(InputValidationError,match='genotype'):
        run_primer_design_plan(plan,tmp_path/'bad_gt',dry_run=True)


@pytest.mark.parametrize('gt',['1','1/1','0|1','1|0'])
def test_alt_genotypes_are_preserved_not_forced_diploid_heterozygous(tmp_path,gt,monkeypatch):
    plan,data,fasta,fai,vcf,sequence=fixture_plan(tmp_path)
    vcf.write_text(vcf.read_text(encoding='utf-8').replace('\t0/1\n',f'\t{gt}\n'),encoding='utf-8')
    monkeypatch.setattr('hifivar.primer_design_plan.Primer3Wrapper.design',fake_design)
    result=run_primer_design_plan(plan,tmp_path/'gt')
    assert result['targets'][0]['source_record']['genotype']==gt


def test_ref_mismatch_against_real_template_fails_before_output(tmp_path):
    plan,data,fasta,fai,vcf,sequence=fixture_plan(tmp_path)
    original=data['targets'][0]['ref']
    replacement=next(x for x in 'ACGT' if x not in (original,data['targets'][0]['alt']))
    vcf.write_text(vcf.read_text(encoding='utf-8').replace(f'\t{original}\t',f'\t{replacement}\t'),encoding='utf-8')
    data['targets'][0]['ref']=replacement;save_plan(plan,data)
    output=tmp_path/'mismatch'
    with pytest.raises(InputValidationError,match='REF allele'):
        run_primer_design_plan(plan,output,dry_run=True)
    assert not output.exists()


@pytest.mark.parametrize('change',['length','duplicate_record','duplicate_sample','missing_contig_length','multiallelic'])
def test_vcf_contract_rejects_ambiguous_or_incompatible_inputs(tmp_path,change):
    plan,data,fasta,fai,vcf,sequence=fixture_plan(tmp_path)
    raw=vcf.read_text(encoding='utf-8')
    if change=='length':raw=raw.replace('length=900','length=901')
    elif change=='missing_contig_length':raw=raw.replace(',length=900','')
    elif change=='duplicate_record':raw+=raw.splitlines()[-1]+'\n'
    elif change=='duplicate_sample':raw=raw.replace('\tHG002\n','\tHG002\tHG002\n').replace('\t0/1\n','\t0/1\t0/1\n')
    else:raw=raw.replace(f"\t{data['targets'][0]['alt']}\t",f"\t{data['targets'][0]['alt']},A\t")
    vcf.write_text(raw,encoding='utf-8')
    with pytest.raises(InputValidationError):run_primer_design_plan(plan,tmp_path/'bad_vcf',dry_run=True)


@pytest.mark.parametrize('change',['duplicate_key','unknown_key','duplicate_target','unknown_setting','unknown_target_key','bool_version'])
def test_strict_yaml_contract(tmp_path,change):
    plan,data,*_=fixture_plan(tmp_path)
    if change=='duplicate_key':
        plan.write_text(plan.read_text(encoding='utf-8')+'\nschema_version: 1\n',encoding='utf-8')
    else:
        if change=='unknown_key':data['unrecognised']=True
        elif change=='duplicate_target':data['targets'].append(dict(data['targets'][0]))
        elif change=='unknown_setting':data['settings']={'primer_really_good':True}
        elif change=='unknown_target_key':data['targets'][0]['confidence']='truth'
        else:data['schema_version']=True
        save_plan(plan,data)
    with pytest.raises(InputValidationError):run_primer_design_plan(plan,tmp_path/'bad_yaml',dry_run=True)


@pytest.mark.parametrize('fai_text',['chr21\t900\t999999\t60\t61\n',
                                    'chr21\t900\t7\t0\t1\n',
                                    'chr21\t900\t7\t60\t6000000\n',
                                    'chr21\t900\t7\t60\t61\nchr21\t900\t7\t60\t61\n'])
def test_invalid_fai_never_creates_output(tmp_path,fai_text):
    plan,data,fasta,fai,vcf,sequence=fixture_plan(tmp_path)
    fai.write_text(fai_text,encoding='utf-8')
    output=tmp_path/'bad_fai'
    with pytest.raises(InputValidationError):run_primer_design_plan(plan,output,dry_run=True)
    assert not output.exists()


def test_huge_fai_stride_fails_before_unbounded_read(tmp_path,monkeypatch):
    fasta=tmp_path/'tiny.fa';fasta.write_bytes(b'>c\nACGT\n')
    original=Path.open
    def guarded_open(self,*args,**kwargs):
        if self==fasta and args and args[0]=='rb':raise AssertionError('invalid stride must fail before reading')
        return original(self,*args,**kwargs)
    monkeypatch.setattr(Path,'open',guarded_open)
    with pytest.raises(InputValidationError,match='bounded|outside'):
        fetch_sequence(fasta,(1000000,3,1,1000000000),1,100000)


def test_existing_output_is_preserved_and_no_native_run(tmp_path,monkeypatch):
    plan,data,*_=fixture_plan(tmp_path)
    output=tmp_path/'existing';output.mkdir();keep=output/'keep.txt';keep.write_bytes(b'old result')
    def forbidden(*args,**kwargs):raise AssertionError('existing output reached native')
    monkeypatch.setattr('hifivar.primer_design_plan.Primer3Wrapper.design',forbidden)
    with pytest.raises(OutputValidationError):run_primer_design_plan(plan,output)
    assert keep.read_bytes()==b'old result'


def test_output_cannot_enclose_inputs(tmp_path):
    plan,data,*_=fixture_plan(tmp_path)
    with pytest.raises(OutputValidationError):run_primer_design_plan(plan,tmp_path,dry_run=True)


def test_changed_source_preserves_partial_evidence_without_completion_receipt(tmp_path,monkeypatch):
    plan,data,fasta,fai,vcf,sequence=fixture_plan(tmp_path)
    def mutate(self,item,output_dir,settings):
        result=fake_design(self,item,output_dir,settings)
        vcf.write_text(vcf.read_text(encoding='utf-8')+'##source_changed\n',encoding='utf-8')
        return result
    monkeypatch.setattr('hifivar.primer_design_plan.Primer3Wrapper.design',mutate)
    output=tmp_path/'partial'
    with pytest.raises(InputValidationError,match='changed'):
        run_primer_design_plan(plan,output)
    assert (output/'site1/fake-native-evidence.txt').is_file()
    assert not (output/'primer_design_summary.json').exists()
    assert not (output/'primers.tsv').exists()


@pytest.mark.parametrize('gt,expected_ploidy,products',[('0/1','2','reference:501;alternate:501'),
                                                     ('1','1','alternate:501'),
                                                     ('1/1','2','alternate:501')])
def test_snv_band_table_preserves_ploidy_and_no_gel_discrimination(tmp_path,monkeypatch,gt,expected_ploidy,products):
    plan,data,fasta,fai,vcf,sequence=fixture_plan(tmp_path)
    vcf.write_text(vcf.read_text(encoding='utf-8').replace('\t0/1\n',f'\t{gt}\n'),encoding='utf-8')
    monkeypatch.setattr('hifivar.primer_design_plan.Primer3Wrapper.design',fake_design)
    output=tmp_path/'bands'
    run_primer_design_plan(plan,output)
    with (output/'predicted_bands.tsv').open(encoding='utf-8',newline='') as handle:
        rows=list(csv.DictReader(handle,delimiter='\t'))
    assert len(rows)==3
    assert rows[0]['source_genotype']==gt and rows[0]['ploidy_from_gt']==expected_ploidy
    assert rows[0]['predicted_allele_products_bp']==products
    assert rows[0]['predicted_length_pattern']=='ONE_LENGTH_CLASS'
    assert rows[0]['gel_discrimination']=='NO_FOR_SNV'
    assert rows[0]['actual_bands']=='NOT_MEASURED'


def test_deletion_band_table_contains_two_computed_lengths_not_observed_bands(tmp_path,monkeypatch):
    plan,data,fasta,fai,vcf,sequence=fixture_plan(tmp_path)
    target=data['targets'][0];old_ref=target['ref'];old_alt=target['alt']
    target['ref']=sequence[300:304];target['alt']=sequence[300]
    vcf.write_text(vcf.read_text(encoding='utf-8').replace(f'\t{old_ref}\t{old_alt}\t',
                   f"\t{target['ref']}\t{target['alt']}\t"),encoding='utf-8')
    save_plan(plan,data)
    monkeypatch.setattr('hifivar.primer_design_plan.Primer3Wrapper.design',fake_design)
    output=tmp_path/'deletion'
    run_primer_design_plan(plan,output)
    with (output/'predicted_bands.tsv').open(encoding='utf-8',newline='') as handle:
        row=next(csv.DictReader(handle,delimiter='\t'))
    assert row['predicted_allele_products_bp']=='reference:501;alternate:498'
    assert row['distinct_product_lengths_bp']=='498;501'
    assert row['predicted_length_pattern']=='TWO_LENGTH_CLASSES'
    assert row['gel_discrimination']=='NOT_ASSESSED' and row['actual_bands']=='NOT_MEASURED'


def test_zero_candidates_produces_honest_empty_aggregate_not_successful_primers(tmp_path,monkeypatch):
    plan,data,*_=fixture_plan(tmp_path)
    def none(self,item,output_dir,settings):
        output_dir.mkdir()
        return evaluate_candidates(item,settings,response(item,()))
    monkeypatch.setattr('hifivar.primer_design_plan.Primer3Wrapper.design',none)
    output=tmp_path/'none'
    result=run_primer_design_plan(plan,output)
    assert result['candidate_pair_count']==0 and result['targets_with_requested_count']==0
    assert result['targets'][0]['status']=='NO_CANDIDATES'
    assert (output/'predicted_amplicons.fasta').read_bytes()==b''
    assert len((output/'predicted_bands.tsv').read_text(encoding='utf-8').splitlines())==1
    with (output / 'primer_sequences.tsv').open(encoding='utf-8', newline='') as handle:
        sequences = csv.DictReader(handle, delimiter='\t')
        assert list(sequences) == []
        assert 'forward_sequence_5to3' in sequences.fieldnames
        assert 'reverse_sequence_5to3' in sequences.fieldnames
