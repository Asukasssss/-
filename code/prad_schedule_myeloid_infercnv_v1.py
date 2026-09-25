"""Run official inferCNV in isolated per-donor directories on server165."""
import argparse, os, subprocess, json, time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import pandas as pd
p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--pilot',action='store_true');a=p.parse_args()
r=a.root;plans=pd.read_csv(r/'private/donor_plan.tsv',sep='\t')
runtime=Path((r/'source/env_path.txt').read_text().strip())
jobs=[(d,m) for d in plans.donor_alias for m in ['primary','reference_split']]
if a.pilot:jobs=jobs[:1]
env=os.environ.copy()
env.update(OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',R_LIBS_USER=str(runtime/'lib/R/library'))
env['PATH']=str(runtime/'bin')+':'+env['PATH']
def run(job):
 d,m=job;out=r/'private'/d/m;out.mkdir(exist_ok=True)
 if (out/'DONE').exists():return dict(donor_alias=d,mode=m,status='DONE',reused=True)
 start=time.time()
 if (out/'run.log').exists():(out/'run.log').rename(out/('run.previous.'+str(time.time_ns())+'.log'))
 with (out/'run.log').open('w') as f:
  proc=subprocess.run([str(runtime/'bin/Rscript'),str(r/'prad_run_myeloid_infercnv_v1.R'),str(r),d,m],stdout=f,stderr=subprocess.STDOUT,env=env)
 result=dict(donor_alias=d,mode=m,status='DONE' if proc.returncode==0 and (out/'DONE').exists() else 'FAILED',returncode=proc.returncode,elapsed_seconds=round(time.time()-start,1))
 (out/'job_status.json').write_text(json.dumps(result,indent=2));return result
results=[]
with ThreadPoolExecutor(max_workers=1 if a.pilot else 3) as ex:
 for f in as_completed([ex.submit(run,j) for j in jobs]):
  z=f.result();results.append(z);print(json.dumps(z),flush=True)
pd.DataFrame(results).to_csv(r/'private'/('pilot_status.tsv' if a.pilot else 'job_status.tsv'),sep='\t',index=False)
if any(z['status']!='DONE' for z in results):raise SystemExit(1)
