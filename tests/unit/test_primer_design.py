"""Synthetic sequence tests establish software contracts, not PCR success."""
from dataclasses import replace
import random

import pytest

from hifivar.exceptions import InputValidationError, OutputValidationError
from hifivar.primer_design import (VariantTemplate, PrimerDesignSettings, build_boulder_input,
                                  evaluate_candidates, parse_boulder_record, reverse_complement)


def template(ref=None, alt=None):
    rng = random.Random(83240)
    sequence = ''.join(rng.choice('ACGT') for _ in range(900))
    ref = ref or sequence[300]
    sequence = sequence[:300] + ref + sequence[300 + len(ref):]
    alt = alt or next(base for base in 'ACGT' if base != ref)
    return VariantTemplate('HG002', 'chr21', 1301, ref, alt, 1001, sequence,
                           'synthetic-record', 'a' * 64)


def response(item, pairs=((100, 600), (110, 610), (120, 620))):
    data = {line.split('=', 1)[0]: line.split('=', 1)[1]
            for line in build_boulder_input(item, PrimerDesignSettings()).splitlines() if line != '='}
    data['PRIMER_PAIR_NUM_RETURNED'] = str(len(pairs))
    for i, (left, right) in enumerate(pairs):
        data.update({f'PRIMER_LEFT_{i}': f'{left},20', f'PRIMER_RIGHT_{i}': f'{right},20',
                     f'PRIMER_LEFT_{i}_SEQUENCE': item.reference_sequence[left:left + 20],
                     f'PRIMER_RIGHT_{i}_SEQUENCE': reverse_complement(item.reference_sequence[right - 19:right + 1]),
                     f'PRIMER_LEFT_{i}_TM': '60.0', f'PRIMER_RIGHT_{i}_TM': '60.5',
                     f'PRIMER_PAIR_{i}_PRODUCT_SIZE': str(right - left + 1),
                     f'PRIMER_PAIR_{i}_PENALTY': str(i + 0.25),
                     f'PRIMER_LEFT_{i}_HAIRPIN_TH': '11.1'})
    return data


def test_snv_does_not_claim_gel_discrimination_and_reconstructs_products():
    item = template()
    result = evaluate_candidates(item, PrimerDesignSettings(), response(item))
    assert result['assay'] == 'PCR_SANGER'
    assert result['snv_gel_discriminates_alleles'] is False
    assert result['candidate_pair_count'] == 3
    assert result['specificity_status'] == 'NOT_ASSESSED'
    assert result['experimental_status'] == 'NOT_RUN'
    pair = result['candidates'][0]
    assert pair['reference_product_bp'] == pair['alternate_product_bp'] == 501
    assert pair['reference_product_sequence'][200] == item.ref
    assert pair['alternate_product_sequence'][200] == item.alt
    assert pair['reverse_sequence_5to3'] == reverse_complement(item.reference_sequence[581:601])
    assert pair['forward_start_1based'] == 1101
    assert pair['reverse_start_1based'] == 1582
    assert pair['reverse_end_1based'] == 1601


@pytest.mark.parametrize('ref,alt,delta,kind', [('ACGT','A',-3,'DEL'), ('A','ACGT',3,'INS')])
def test_indel_products_are_allele_reconstructed(ref,alt,delta,kind):
    item=template(ref,alt)
    result=evaluate_candidates(item,PrimerDesignSettings(),response(item))
    pair=result['candidates'][0]
    assert result['assay']=='FLANKING_PCR' and result['variant']['variant_type']==kind
    assert pair['alternate_product_bp']-pair['reference_product_bp']==delta
    assert pair['alternate_product_sequence']==pair['reference_product_sequence'][:200]+alt+pair['reference_product_sequence'][200+len(ref):]
    assert pair['alternate_product_sequence'].startswith(pair['forward_sequence_5to3'])
    assert pair['alternate_product_sequence'].endswith(reverse_complement(pair['reverse_sequence_5to3']))


@pytest.mark.parametrize('field,value', [('position',True),('position',0),('sequence_start',1500),
                                      ('ref','N'),('alt','<INS>'),('alt','A[chr2:1['),
                                      ('sample','x\ny'),('source_variant_id',''),('source_record_sha256','bad')])
def test_invalid_variant_identity_is_rejected(field,value):
    with pytest.raises(InputValidationError):replace(template(),**{field:value})


def test_ref_mismatch_complex_substitution_and_empty_flank_rejected():
    item=template()
    with pytest.raises(InputValidationError,match='match'):
        replace(item,ref=next(x for x in 'ACGT' if x not in (item.ref,item.alt)))
    with pytest.raises(InputValidationError,match='Only'):
        replace(item,ref='AC',alt='GT')
    with pytest.raises(InputValidationError,match='flanking'):
        replace(item,sequence_start=item.position)


@pytest.mark.parametrize('field,value', [('requested_pairs',True),('requested_pairs',0),('requested_pairs',21),
                                      ('timeout_seconds',float('nan')),('timeout_seconds',0),('primer_min_size',12),
                                      ('primer_min_tm',70),('primer_max_gc',10),('min_pair_shift',0),
                                      ('product_size_ranges',((800,200),)),('product_size_ranges',()),
                                      ('thermodynamic_parameters_path','p\nP3_FILE_FLAG=1')])
def test_invalid_settings_are_rejected(field,value):
    with pytest.raises(InputValidationError):PrimerDesignSettings(**{field:value})


def test_boulder_input_constrains_unchanged_binding_flanks():
    item=template('ACGT','A');text=build_boulder_input(item,PrimerDesignSettings())
    assert 'SEQUENCE_TARGET=300,4\n' in text
    assert 'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST=0,300,304,596\n' in text
    assert 'PRIMER_NUM_RETURN=15\n' in text
    assert 'P3_FILE_FLAG=0\n' in text
    assert text.endswith('\n=\n')


@pytest.mark.parametrize('text', ['A=1\n','A=1\nA=2\n=\n','A=1\n=\nB=2\n=\n',
                                 'garbage\n=\n','PRIMER_ERROR=bad template\n=\n'])
def test_malformed_or_error_boulder_output_fails(text):
    with pytest.raises(OutputValidationError):parse_boulder_record(text)


def test_output_values_with_equals_are_preserved():
    assert parse_boulder_record('PRIMER_WARNING=a=b\n=\n')['PRIMER_WARNING']=='a=b'


@pytest.mark.parametrize('change', ['right_sequence','right_position','left_overlap','product_size','nonfinite_tm','negative_penalty'])
def test_corrupt_native_coordinates_or_metrics_fail(change):
    item=template();tags=response(item)
    if change=='right_sequence':tags['PRIMER_RIGHT_0_SEQUENCE']='A'*20
    elif change=='right_position':tags['PRIMER_RIGHT_0']='0,20'
    elif change=='left_overlap':tags['PRIMER_LEFT_0']='290,20'
    elif change=='product_size':tags['PRIMER_PAIR_0_PRODUCT_SIZE']='500'
    elif change=='nonfinite_tm':tags['PRIMER_LEFT_0_TM']='nan'
    else:tags['PRIMER_PAIR_0_PENALTY']='-1'
    with pytest.raises(OutputValidationError):evaluate_candidates(item,PrimerDesignSettings(),tags)


def test_positional_diversity_and_shortfall_are_honest():
    item=template();tags=response(item,((100,600),(101,601),(100,600)))
    result=evaluate_candidates(item,PrimerDesignSettings(),tags)
    assert result['candidate_pair_count']==1 and result['requested_count_satisfied'] is False
    assert len(result['rejected_candidates'])==2


def test_no_candidates_is_valid_design_failure_not_fake_pairs():
    item=template();tags=response(item,())
    tags['PRIMER_LEFT_EXPLAIN']='no acceptable primers'
    result=evaluate_candidates(item,PrimerDesignSettings(),tags)
    assert result['status']=='NO_CANDIDATES' and result['candidates']==[]
    assert result['primer3_explanations']['PRIMER_LEFT_EXPLAIN']=='no acceptable primers'


def test_alternate_product_outside_user_range_is_rejected():
    item=template('A','A'+'C'*500)
    result=evaluate_candidates(item,PrimerDesignSettings(),response(item))
    assert result['candidate_pair_count']==0
    assert all(x['reason']=='ALTERNATE_PRODUCT_OUTSIDE_REQUESTED_RANGE' for x in result['rejected_candidates'])
