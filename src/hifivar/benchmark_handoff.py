"""Explicit source-to-truth benchmark consumer for the owned analysis launcher.

Receipts bind execution and context-specific metric artifacts, not universal
accuracy, biological identity, clinical utility or truth-resource authenticity.
Independent Phase13 rules keep their original interface.
"""
from __future__ import annotations

import copy
import gzip
from dataclasses import asdict
import os
from pathlib import Path
from typing import Mapping

import yaml

from hifivar.annotation import audit_annotation_vcf_identity
from hifivar.attempt_receipts import artifact_identity, _path, _read
from hifivar.benchmark_inputs import audit_benchmark_regions
from hifivar.downstream_plan import (_PlanLoader, _ROW_KEYS, _configured_catalog, _mapping,
    _plan_bytes, _row, _unknown_keys)
from hifivar.exceptions import ConfigurationError, InputValidationError, OutputValidationError
from hifivar.serialization import standardize_data, write_json_atomic, write_yaml_atomic

TRACKS=('small_variants','sv','assembly_sv','tr')
CATEGORY={'small_variants':'small','sv':'sv','assembly_sv':'sv','tr':'tr'}
CONTROL_NAMES=('effective.yaml','source_bindings.json')


def enabled(config: Mapping) -> bool:
    return config.get('benchmark',{}).get('enabled') is True


def without_benchmark(config: dict) -> dict:
    result=copy.deepcopy(config)
    result['benchmark'].update(enabled=False,source_plan=None)
    for track in TRACKS:result['benchmark'][track]['enabled']=False
    return result


def read_source_plan(path: Path, config: dict, context):
    """Read one explicit query per enabled track for one declared truth sample."""
    source,raw=_plan_bytes(path)
    try:
        doc=_mapping(yaml.load(raw.decode('utf-8-sig'),Loader=_PlanLoader),'Benchmark plan')
    except (UnicodeError,yaml.YAMLError,RecursionError) as error:
        raise InputValidationError('Benchmark source plan must be bounded strict UTF8 YAML.') from error
    _unknown_keys(doc,frozenset({'schema_version','benchmark'}),'Benchmark plan')
    if type(doc.get('schema_version')) is not int or doc['schema_version']!=1:
        raise InputValidationError('Benchmark source schema_version must be integer1.')
    rows=_mapping(doc.get('benchmark'),'Benchmark sources')
    selected={name for name in TRACKS if config['benchmark'][name].get('enabled')}
    if not selected or set(rows)!=selected:
        raise InputValidationError('Benchmark plan must cover each enabled track exactly once.')
    sample=config['benchmark']['sample_id']
    if sample not in context.sample_ids:
        raise InputValidationError('Benchmark sample_id is absent from the analysis context.')
    catalog=_configured_catalog(config,source.parent) if 'tr' in rows else None
    parsed={}
    for name in TRACKS:
        if name not in rows:continue
        fields=_mapping(rows[name],'Benchmark source')
        _unknown_keys(fields,_ROW_KEYS-{'state','variant_ids','bnd_bundle_manifest'},'Benchmark source')
        row=_row(fields,'benchmark',context.sample_ids,context.reference.build,config,source.parent,catalog)
        if row.sample_id!=sample or row.track!=CATEGORY[name] or row.output_role in {'gvcf','bnd_vcf'}:
            raise InputValidationError('Benchmark query must select the declared sample/category native VCF, not gVCF/BND-derived role.')
        if name=='assembly_sv' and row.origin=='GENERATED' and not row.producer.startswith('assembly_sv.'):
            raise InputValidationError('Generated assembly benchmark requires an assembly-SV producer.')
        if name=='tr' and row.catalog_id!=config['benchmark']['tr'].get('catalog_id'):
            raise InputValidationError('Benchmark TR truth/query catalog SHA256 must match explicitly.')
        parsed[name]=row
    return parsed


def _query(path: Path, index: Path, sample: str, context) -> dict:
    if not str(path).endswith('.vcf.gz') or index!=Path(str(path)+'.tbi'):
        raise InputValidationError('Owned benchmark requires BGZF VCF with conventional TBI.')
    from hifivar.launcher import _ready_vcf
    _ready_vcf(path,sample,context)
    audit=audit_annotation_vcf_identity(path,sample,context.reference)
    if audit['explicit_bnd_derived_header']:
        raise InputValidationError('Benchmark consumes native calls, not derived BND review representations.')
    # Renaming a gVCF does not make it an ordinary query. A bare '*' is a
    # legitimate spanning-deletion allele and is distinct from '<*>'.
    records=0
    try:
        with gzip.open(path,'rt',encoding='utf-8') as handle:
            for line in handle:
                if line.startswith('#'):continue
                fields=line.rstrip('\r\n').split('\t',8)
                if len(fields)<8:
                    raise InputValidationError('Benchmark VCF record requires at least eight columns.')
                if {'<NON_REF>','<*>'}.intersection(fields[4].split(',')):
                    raise InputValidationError('Benchmark requires native variant VCF, not gVCF reference-block alleles.')
                records+=1
    except (OSError,EOFError,UnicodeError) as error:
        raise InputValidationError('Unable to stream benchmark VCF record alleles.') from error
    audit.update(gvcf_record_check='PASS_NO_RESERVED_REFERENCE_BLOCK_ALLELES',
                 records_scanned_for_gvcf=records,record_validation='ALT_COLUMN_GVCF_BLOCK_CHECK_ONLY')
    return audit


def prepare_intent(config: dict, context, base: Path, root: Path):
    """Normalize/freeze query intent, truth and region resources before mkdir."""
    if not enabled(config):return {},[],{}
    section=config['benchmark']
    from hifivar.sample import validate_sample_id
    validate_sample_id(section['benchmark_id'])
    if not section.get('source_plan') or section.get('overwrite'):
        raise ConfigurationError('hifivar run benchmark requires source_plan and overwrite=false.')
    from hifivar.launcher import _absolute
    if config.get('tr',{}).get('catalog'):
        config['tr']['catalog']=_absolute(config['tr']['catalog'],base)
    path=Path(_absolute(section['source_plan'],base))
    rows=read_source_plan(path,config,context)
    inputs=[path];audits={}
    if config.get('tr',{}).get('catalog'):inputs.append(Path(config['tr']['catalog']))
    document={}
    for name,row in rows.items():
        track=section[name]
        if track.get('query_vcf') is not None:
            raise ConfigurationError('source_plan and legacy query_vcf are mutually exclusive.')
        for key in ('truth_vcf','confident_bed'):
            if track.get(key):track[key]=_absolute(track[key],base)
        for item in track.get('stratifications',[]):item['path']=_absolute(item['path'],base)
        truth=Path(track['truth_vcf']);index=Path(str(truth)+'.tbi')
        audit=_query(truth,index,row.sample_id,context)
        inputs.extend((truth,index));regions=[]
        if name in {'sv','assembly_sv'} and not track.get('confident_bed'):
            raise ConfigurationError('Owned SV benchmark requires an explicit confident_bed and its version.')
        if track.get('confident_bed'):
            if not track.get('confident_bed_version'):
                raise ConfigurationError('Owned benchmark BED version must be explicit.')
            regions.append({'name':'confident','path':track['confident_bed'],'version':track['confident_bed_version']})
        regions+=track.get('stratifications',[])
        for item in regions:
            if str(item['version']).strip().casefold()=='latest':raise ConfigurationError('Benchmark BED version must be pinned.')
            audit_benchmark_regions(Path(item['path']),context.reference);inputs.append(Path(item['path']))
        if track['truth_version'].strip().casefold()=='latest':raise ConfigurationError('Truth version must be pinned.')
        if row.origin=='EXTERNAL':
            _query(row.path,row.index_path,row.sample_id,context)
            if os.path.samefile(row.path,truth):raise InputValidationError('Query and truth must not alias the same file.')
            inputs.extend((row.path,row.index_path))
        fields=standardize_data(asdict(row));fields.pop('state');fields.pop('variant_ids');fields.pop('bnd_bundle_manifest')
        document[name]=fields
        audits[name]={'truth_header':audit,'regions':regions,'truth_authenticity':'USER_DECLARED_NOT_PROVIDER_AUTHENTICATED',
                      'REF_sequence_check':'NOT_RUN_HEADER_COMPATIBILITY_ONLY'}
    section['source_plan']=str(root/'benchmark_source_plan.yaml')
    return {'benchmark_source_plan.yaml':yaml.safe_dump({'schema_version':1,'benchmark':document},sort_keys=False)},inputs,audits


def resolve_sources(config,context,root,attempt):
    rows=read_source_plan(root/'benchmark_source_plan.yaml',config,context);result={}
    frozen={x['path']:x for x in _read(root/'launch_manifest.json')['identity']['inputs']}
    for name,row in rows.items():
        path,index,tool,version,bundle=row.path,row.index_path,row.source_tool,row.source_version,None
        link=None
        if row.origin=='GENERATED':
            from hifivar.downstream_handoff import _generated_source
            path,index,tool,version,link,bundle=_generated_source(root,attempt,row,context.reference)
        if bundle is not None:raise InputValidationError('BND-derived benchmark source is not a native query.')
        query_audit=_query(path,index,row.sample_id,context)
        truth=Path(config['benchmark'][name]['truth_vcf']);truth_index=Path(str(truth)+'.tbi')
        truth_audit=_query(truth,truth_index,row.sample_id,context)
        if os.path.samefile(path,truth):raise InputValidationError('Query/truth source aliases are forbidden.')
        roles={'reference_fasta':artifact_identity(context.reference.fasta),'reference_fai':artifact_identity(context.reference.fai),
               'vcf':artifact_identity(path),'index':artifact_identity(index),
               'truth_vcf':artifact_identity(truth),'truth_index':artifact_identity(truth_index)}
        track=config['benchmark'][name]
        regions=[]
        if track.get('confident_bed'):regions.append(Path(track['confident_bed']))
        regions.extend(Path(x['path']) for x in track.get('stratifications',[]))
        for i,region in enumerate(regions):
            audit_benchmark_regions(region,context.reference);roles['region_'+str(i)]=artifact_identity(region)
        for role,info in roles.items():
            if (row.origin=='EXTERNAL' or role not in {'vcf','index'}) and frozen.get(info['path'])!=info:
                raise InputValidationError('Benchmark query/truth/region differs from frozen launch input.')
        result[name]={'source_inputs':roles,'origin':row.origin,'producer_receipt':link,'catalog_id':row.catalog_id,
          'tool':tool,'source_version':version,'source_version_scope':'SELECTED_ATTEMPT_RECEIPT' if link else 'USER_DECLARED_NOT_REPROBED',
          'query_header':query_audit,'truth_header':truth_audit,'truth_metadata':{k:track.get(k) for k in ('truth_version','truth_source','confident_bed_version','region_class')},
          'REF_sequence_check':'NOT_RUN_HEADER_COMPATIBILITY_ONLY'}
    return standardize_data(result)


def materialize(config,context,root,attempt,*,calling_config):
    from hifivar.downstream_handoff import consumer_config
    result=resolve_sources(config,context,root,attempt)
    output=_path(attempt/'benchmark',root=root)
    if output.exists():raise OutputValidationError('Benchmark materialization exists; use a fresh attempt.')
    effective=consumer_config(calling_config)
    for group in ('annotation','cohort','review'):effective[group]['enabled']=False
    effective['benchmark']=copy.deepcopy(config['benchmark']);effective['benchmark']['source_plan']=None
    for name,source in result.items():effective['benchmark'][name]['query_vcf']=source['source_inputs']['vcf']['path']
    from hifivar.config import validate_config
    validate_config(effective,require_complete=True)
    output.mkdir(parents=True)
    write_yaml_atomic(effective,output/'effective.yaml');write_json_atomic(result,output/'source_bindings.json')
    record={'kind':'hifivar_benchmark_materialization','schema_version':1,'attempt_id':attempt.name,
      'context_sha256':artifact_identity(attempt/'execution_context.json')['sha256'],
      'plan':artifact_identity(root/'benchmark_source_plan.yaml'),
      'controls':{n:artifact_identity(output/n) for n in CONTROL_NAMES},
      'sample_id':config['benchmark']['sample_id'],'source_receipts':[v['producer_receipt'] for v in result.values() if v['producer_receipt']]}
    write_json_atomic(record,output/'materialization.json')
    return output/'effective.yaml',record


def validate_materialization(root,attempt):
    from hifivar.attempt_receipts import _context
    execution,identity,_=_context(attempt/'execution_context.json',live=False)
    if execution['dry_run'] or not identity['config']['benchmark'].get('source_plan'):
        raise InputValidationError('Benchmark execution requires a native source-plan attempt.')
    path=attempt/'benchmark';record=_read(path/'materialization.json');status=_read(attempt/'status.json')
    if (record.get('kind')!='hifivar_benchmark_materialization' or record.get('schema_version')!=1 or record.get('attempt_id')!=attempt.name
        or record.get('context_sha256')!=artifact_identity(attempt/'execution_context.json')['sha256']
        or record.get('plan')!=artifact_identity(root/'benchmark_source_plan.yaml')
        or status.get('benchmark_materialization_sha256')!=artifact_identity(path/'materialization.json')['sha256']):
        raise InputValidationError('Benchmark controls do not belong to the selected attempt.')
    if set(record.get('controls',{}))!=set(CONTROL_NAMES):raise InputValidationError('Incomplete benchmark controls.')
    for name,info in record['controls'].items():
        if artifact_identity(path/name)!=info:raise InputValidationError('Benchmark controls changed.')
    effective=yaml.safe_load((path/'effective.yaml').read_text(encoding='utf-8'))
    from hifivar.context import AnalysisContext
    context=AnalysisContext.from_config(effective)
    if effective['reference']!=identity['config']['reference']:
        raise InputValidationError('Benchmark reference differs from frozen analysis context.')
    fresh=resolve_sources(identity['config'],context,root,attempt)
    if _read(path/'source_bindings.json')!=fresh:raise InputValidationError('Benchmark sources changed after materialization.')
    from hifivar.downstream_handoff import consumer_config
    # Reconstruct the only permissible consumer config from its actual BAM context.
    expected=consumer_config(effective)
    expected['benchmark']=copy.deepcopy(identity['config']['benchmark']);expected['benchmark']['source_plan']=None
    for group in ('annotation','cohort','review'):expected[group]['enabled']=False
    for name,source in fresh.items():expected['benchmark'][name]['query_vcf']=source['source_inputs']['vcf']['path']
    if effective!=expected:raise InputValidationError('Benchmark consumer reactivates producers or changes its policy.')
    if record['source_receipts']!=[v['producer_receipt'] for v in fresh.values() if v['producer_receipt']]:
        raise InputValidationError('Benchmark producer receipt links changed.')
    return record,effective,fresh


def outputs(config,root):
    b=Path(root)/'benchmark'/config['benchmark']['benchmark_id'];result={'manifest':b/'benchmark_manifest.json','metrics':b/'benchmark_metrics.tsv'}
    for name in TRACKS:
        track=config['benchmark'][name]
        if not track.get('enabled'):continue
        sub=b/name;result[name+'_result']=sub/'benchmark_result.json'
        if name=='small_variants':
            result['small_summary']=sub/'happy.summary.csv'
            result['small_metrics']=sub/('happy.metrics.json.gz' if track['metrics_compression']=='gzip' else 'happy.metrics.json')
        elif name=='tr':result['tr_loci']=sub/'locus_comparison.tsv'
        else:
            for role,file in {'summary':'summary.json','tp_base':'tp-base.vcf.gz','tp_comp':'tp-comp.vcf.gz','fp':'fp.vcf.gz','fn':'fn.vcf.gz'}.items():
                result[name+'_'+role]=sub/'truvari'/file
    return result


def input_paths(config):
    result=[config['reference']['fasta'],str(config['reference']['fasta'])+'.fai']
    for name in TRACKS:
        t=config['benchmark'][name]
        if not t.get('enabled'):continue
        for key in ('query_vcf','truth_vcf'):
            result.extend((t[key],str(t[key])+'.tbi'))
        if t.get('confident_bed'):result.append(t['confident_bed'])
        result.extend(x['path'] for x in t.get('stratifications',[]))
    return list(dict.fromkeys(result))


def receipt_binding(context_path,context,identity,root,samples,inputs,input_manifest,*,live):
    record,effective,sources=validate_materialization(root,context_path.parent)
    selected=(identity['config']['benchmark']['sample_id'],)
    expected={name:{role:Path(info['path']) for role,info in row['source_inputs'].items()} for name,row in sources.items()}
    if tuple(samples)!=selected or Path(input_manifest or '')!=context_path.parent/'benchmark/source_bindings.json' or inputs!=expected:
        raise InputValidationError('Benchmark actual source/selected sample differs from its owned materialization.')
    return {'mode':'EXPLICIT_BENCHMARK_SOURCE_PLAN','source_inputs':{name:row['source_inputs'] for name,row in sources.items()},
      'input_manifest':artifact_identity(Path(input_manifest)),'controls':record['controls'],
      'materialization':artifact_identity(context_path.parent/'benchmark/materialization.json'),
      'source_receipts':record['source_receipts'],'sample_id':selected[0],
      'expected_execution_kind':'IN_PROCESS' if set(sources)=={'tr'} else 'NATIVE',
      'expected_outputs':{role:str(path) for role,path in outputs(effective,root/'results').items()}}


def execute(config,context,*,happy_wrapper=None,truvari_wrapper=None):
    """Execute existing comparison engines serially, retain metrics and sources."""
    from hifivar.attempt_receipts import ATTEMPT_CONTEXT_ENV, attempt_step
    from hifivar.benchmark import BenchmarkRegion,TruthSet,compare_tr_vcfs,BenchmarkVariantClass
    from hifivar.happy import HappyRequest,HappyWrapper
    from hifivar.truvari import TruvariRequest,TruvariThresholds,TruvariWrapper,parse_truvari_summary,stratify_truvari_outputs
    from hifivar.serialization import write_json_atomic
    import csv
    context_path=Path(os.environ.get(ATTEMPT_CONTEXT_ENV,''))
    if not context_path.is_file():raise InputValidationError('Owned benchmark requires the launcher execution context.')
    root=context_path.parent.parent.parent;record,effective,sources=validate_materialization(root,context_path.parent)
    if config!=effective:raise InputValidationError('Actual benchmark script config differs from frozen consumer controls.')
    out=outputs(effective,root/'results');section=effective['benchmark'];sample=section['sample_id']
    source_inputs={n:{r:Path(i['path']) for r,i in x['source_inputs'].items()} for n,x in sources.items()}
    versions={};results=[]
    with attempt_step('benchmark',(sample,),outputs=out,inputs=source_inputs,input_manifest=context_path.parent/'benchmark/source_bindings.json') as step:
        for name in TRACKS:
            if name not in sources:continue
            t=section[name];sub=out[name+'_result'].parent
            truth=TruthSet(Path(t['truth_vcf']),t['truth_version'],context.reference.build,t['truth_source'],t.get('catalog_id'))
            if name=='small_variants':
                region=BenchmarkRegion('confident',Path(t['confident_bed']),t['confident_bed_version'])
                strata=tuple(BenchmarkRegion(x['name'],Path(x['path']),x['version'],x['region_class']) for x in t['stratifications'])
                request=HappyRequest(section['benchmark_id'],sample,context.reference,Path(t['query_vcf']),truth,region,sub/'happy',t['engine'],t['threads'],t['summary_filter'],strata,False)
                result=(happy_wrapper or HappyWrapper(executable=t['happy_executable'],configured_version=t.get('happy_version'))).run(request,stderr_path=root/'logs/benchmark/happy.stderr.log')
                if result.status.value!='completed' or result.metrics_artifact!=out['small_metrics']:raise OutputValidationError('hap.py did not complete its expected metric contract.')
                payload=result.to_dict();payload['execution_status']=payload.pop('status');payload.update(status='PASS',variant_class='small_variant',query_path=t['query_vcf'],tool='hap.py',summary_filter=t['summary_filter'])
                versions['happy']=result.version
            elif name=='tr':
                result=compare_tr_vcfs(benchmark_id=section['benchmark_id'],sample_id=sample,query_vcf=Path(t['query_vcf']),truth_set=truth,query_catalog_id=sources[name]['catalog_id'],output_tsv=out['tr_loci'],overwrite=False)
                payload=result.to_dict();versions['tr_exact']='1'
            else:
                policy=TruvariThresholds(**{k:t.get(k) for k in ('refdist','pctseq','pctsize','pctovl','sizemin','sizemax','bnddist')},pass_only=t['pass_only'])
                request=TruvariRequest(sample,context.reference,Path(t['truth_vcf']),Path(t['query_vcf']),sub/'truvari',False,Path(t['confident_bed']),policy)
                result=(truvari_wrapper or TruvariWrapper(executable=t['truvari_executable'])).run(request,stderr_path=root/('logs/benchmark/'+name+'.stderr.log'))
                if result.status.value!='completed':raise OutputValidationError('Truvari did not complete.')
                kind=BenchmarkVariantClass.ASSEMBLY_SV if name=='assembly_sv' else BenchmarkVariantClass.SV
                metrics=parse_truvari_summary(result.summary_path,variant_class=kind)
                stratified,unsupported=stratify_truvari_outputs(request.output_directory,size_bins=tuple(t['size_bins']),variant_class=kind)
                payload=result.to_dict();payload['execution_status']=payload.pop('status');payload.update(benchmark_id=section['benchmark_id'],sample_id=sample,variant_class=name,
                  status='PARTIAL' if unsupported else 'PASS',query_path=t['query_vcf'],truth_set=truth.to_dict(),tool='truvari',metrics=[m.to_dict() for m in metrics+stratified],unsupported_length_stratification_types=list(unsupported))
                versions[name+'_truvari']=result.version
            payload['source_binding']=sources[name];payload['benchmark_policy']=t
            write_json_atomic(payload,out[name+'_result']);results.append(payload)
        write_json_atomic({'benchmark_id':section['benchmark_id'],'sample_id':sample,'reference_build':context.reference.build,'results':results,
          'scope':'EXPLICIT_QUERY_TRUTH_REGION_POLICY_NOT_UNIVERSAL_ACCURACY','raw_query_modified':False},out['manifest'])
        with out['metrics'].open('x',encoding='utf-8',newline='') as f:
            writer=csv.writer(f,delimiter='\t',lineterminator='\n');writer.writerow(('benchmark_id','sample','variant_class','status','tool','stratum','metric','value','source_field'))
            for result in results:
                for m in result.get('metrics',[]):writer.writerow((section['benchmark_id'],sample,result['variant_class'],result['status'],result['tool'],m.get('stratum'),m['name'],m['value'],m.get('source_field')))
        step.complete(outputs=out,tool_versions=versions,validation={'validated':True,'scope':'EXISTING_ENGINE_METRIC_ARTIFACTS_CONTEXT_BOUND'},
          execution_kind='IN_PROCESS' if set(sources)=={'tr'} else 'NATIVE',metadata={'reference_build':context.reference.build,
           'benchmark_statuses':{r['variant_class']:r['status'] for r in results},'truth_authenticity':'USER_DECLARED_NOT_PROVIDER_AUTHENTICATED',
           'REF_sequence_check':'NOT_RUN_HEADER_COMPATIBILITY_ONLY','raw_queries_modified':False,'universal_accuracy':False})
    return results
