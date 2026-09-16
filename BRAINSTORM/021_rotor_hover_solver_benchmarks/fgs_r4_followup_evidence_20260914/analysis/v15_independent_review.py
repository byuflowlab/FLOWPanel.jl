import csv, pathlib, statistics as S, math, tomllib, hashlib, json
base=pathlib.Path('BRAINSTORM/021_rotor_hover_solver_benchmarks/fgs_r4_followup_evidence_20260914'); root=base/'diag-v15-13694724'; out=base/'analysis'
read=lambda p:list(csv.DictReader(p.open()))
num=lambda r,k:float(r[k])
def save(name, rows):
 with (out/name).open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
stages=['initialization','fmm','influence_mapping','residual','leaf_solve','nonself_product','scatter','remaining_iteration','final_update']
summ=[]; batches=[]; budget=[]; checks=[]; baseline_hashes=None
for j in [1,4,8,16,32,64]:
 p=root/f'j{j}-b1/results'; r=read(p/'diagnostic_trials.csv'); eq=read(p/'instrumentation_equivalence.csv'); pv=read(p/'cpu_thread_complete_validation.csv')
 assert len(r)==40 and len(eq)==2 and len(pv)==1
 for x in r+eq+pv:
  assert x['accepted']=='true' and x['authoritative_evaluator']=='certified_fmm'
  assert 0<=num(x,'authoritative_rel_l2')<=1e-6
  for k in ['relative_solution_delta','solution_delta','counterpart_solution_rel_l2']:
   if k in x: assert math.isfinite(num(x,k)) and 0<=num(x,k)<=1e-8
  if 'solved' in x: assert x['solved']=='true'
 for x in eq:
  assert 0<=num(x,'direct_rel_l2')<=1e-6 and 0<=num(x,'direct_fmm_delta')<=1e-7
  assert x['counterpart_history_identical']=='true' and x['history_length']=='28' and x['iterations']=='27'
 assert pv[0]['eligible']=='true' and pv[0]['fmm_certified']=='true'
 assert num(pv[0],'retained_bytes')<500*2**30 and num(pv[0],'process_peak_rss_bytes')<500*2**30
 u=[x for x in r if x['instrumented']=='false']; ins=[x for x in r if x['instrumented']=='true']
 med=lambda rr,k:S.median(num(x,k) for x in rr)
 um=med(u,'solve_seconds'); im=med(ins,'solve_seconds')
 summ.append(dict(j=j,uninstrumented_median_s=um,uninstrumented_min_s=min(num(x,'solve_seconds') for x in u),uninstrumented_max_s=max(num(x,'solve_seconds') for x in u),instrumented_median_s=im,overhead_percent=100*(im/um-1),speedup_j1=37.5593957/um,instrumented_median_excluding_b2t1_s=med([x for x in ins if not(x['batch']=='2' and x['trial']=='1')],'solve_seconds'),first_instrumented_outer_gap_s=num(ins[0],'solve_seconds')-num(ins[0],'total_stage_seconds')))
 for b in range(1,5):
  rr=[x for x in r if int(x['batch'])==b];assert len(rr)==10 and sorted(int(x['trial']) for x in rr)==list(range(1,11));assert all(x['instrumented']==('true' if b%2==0 else 'false') for x in rr)
  batches.append(dict(j=j,batch=b,instrumented=b%2==0,median_s=med(rr,'solve_seconds'),min_s=min(num(x,'solve_seconds') for x in rr),max_s=max(num(x,'solve_seconds') for x in rr)))
 errors=[]
 for x in ins:
  assert x['outer_count']=='28' and x['sweep_count']=='81' and x['leaf_visit_count']=='86508'
  vals=[num(x,k+'_seconds') for k in stages];assert all(math.isfinite(v) and v>=0 for v in vals)
  err=max(abs(sum(vals)-num(x,'exclusive_sum_seconds')),abs(num(x,'total_stage_seconds')-num(x,'exclusive_sum_seconds')-num(x,'unaccounted_seconds')));errors.append(err);assert err<1e-6
  assert 0<=num(x,'unaccounted_seconds')<.1 and 0<=num(x,'solve_seconds')-num(x,'total_stage_seconds')<1
 for k in stages+['exclusive_sum','unaccounted','total_stage']:
  budget.append(dict(j=j,stage=k,median_seconds=med(ins,k+'_seconds'),median_fraction_of_instrumented_solve=S.median(num(x,k+'_seconds')/num(x,'solve_seconds') for x in ins)))
 hashes={n:hashlib.sha256((p/n).read_bytes()).hexdigest() for n in ['gemv_census.csv','gemv_census.toml','dependency_edges.csv']}
 if baseline_hashes is None:baseline_hashes=hashes
 assert hashes==baseline_hashes
 checks.append(dict(j=j,timed_rows=len(r),direct_controls=len(eq),max_stage_reconciliation_error_s=max(errors),bc_rel_l2=max(num(x,'authoritative_rel_l2') for x in r+eq+pv),direct_bc_rel_l2=max(num(x,'direct_rel_l2') for x in eq),direct_fmm_delta=max(num(x,'direct_fmm_delta') for x in eq),max_repeat_delta=max(num(x,'relative_solution_delta') for x in r),history_length=28,outer_checks=28,update_iterations=27,sweeps=81,census_identical=True))
p=root/'j1-b1/results';c=read(p/'gemv_census.csv'); e=read(p/'dependency_edges.csv');t=tomllib.loads((p/'gemv_census.toml').read_text());edges={(int(x['source_leaf']),int(x['dependent_leaf'])) for x in e}; assert len(edges)==len(e)==95390
leaves={int(x['leaf']):x for x in c};adj={k:set() for k in leaves};outdeg={k:0 for k in leaves}
for a,b in edges: assert a!=b;adj[a].add(b);adj[b].add(a);outdeg[a]+=1
for k,x in leaves.items():
 assert len(adj[k])==int(x['conflict_degree']) and outdeg[k]==int(x['dependent_leaves'])
 assert int(x['matrix_elements'])==int(x['m'])*int(x['n']) and int(x['matrix_bytes'])==8*int(x['matrix_elements']) and int(x['scatter_entries'])==int(x['m'])
 assert all(x['potential_color']!=leaves[v]['potential_color'] for v in adj[k])
assert sum(len(v) for v in adj.values())//2==48627
assert sum(int(x['matrix_bytes']) for x in c)==t['matrix_bytes']==2862850032
assert sum(int(x['scatter_entries']) for x in c)==t['scatter_entries']==5431340
colors={int(x['potential_color']) for x in c};assert len(colors)==79
for color in colors:
 cc=[x for x in c if int(x['potential_color'])==color];assert len(cc)==t['potential_color_sizes'][color-1];assert sum(int(x['matrix_bytes']) for x in cc)==t['potential_color_matrix_bytes'][color-1]
def quant(v,q):
 a=sorted(v); f=(len(a)-1)*q; lo=int(f);return a[lo]+(a[min(lo+1,len(a)-1)]-a[lo])*(f-lo)
census=[]
for k in ['m','n','matrix_bytes','scatter_entries','target_interactions','dependent_leaves','conflict_degree']:
 v=[int(x[k]) for x in c];census.append(dict(metric=k,count=len(v),total=sum(v),minimum=min(v),p25=quant(v,.25),median=quant(v,.5),p75=quant(v,.75),p95=quant(v,.95),maximum=max(v)))
v=t['potential_color_sizes'];census.append(dict(metric='potential_color_size',count=len(v),total=sum(v),minimum=min(v),p25=quant(v,.25),median=quant(v,.5),p75=quant(v,.75),p95=quant(v,.95),maximum=max(v)))
for name,rows in [('v15_independent_scaling.csv',summ),('v15_independent_batches.csv',batches),('v15_independent_stage_budget.csv',budget),('v15_independent_checks.csv',checks),('v15_independent_census.csv',census)]:save(name,rows)
print(json.dumps({'scaling':summ,'max_reconciliation_error':max(x['max_stage_reconciliation_error_s'] for x in checks),'census':census},indent=2))
