"""Bounded server runner: one pilot then at most three independent R processes."""
from pathlib import Path
import concurrent.futures,json,os,subprocess,time
import pandas as pd
R=Path(__file__).resolve().parent
RSCRIPT=R/'env/bin/Rscript'
def run_unit(p):
    p=Path(p)
    if (p/'INFERENCE_DONE').exists():return {'unit_path':str(p),'status':'DONE_REUSED'}
    env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='2',MKL_NUM_THREADS='1')
    with (p/'infercnv.log').open('w') as log:
        z=subprocess.run([str(RSCRIPT),str(R/'myeloid_infercnv_run.R'),str(p),str(R)],env=env,stdout=log,stderr=subprocess.STDOUT)
    status='DONE' if z.returncode==0 and (p/'INFERENCE_DONE').exists() else 'FAILED'
    return {'unit_path':str(p),'status':status,'exit_code':z.returncode}
if __name__=='__main__':
    assert RSCRIPT.exists()
    chk=subprocess.run([str(RSCRIPT),'-e','stopifnot(requireNamespace("infercnv",quietly=TRUE));cat(as.character(packageVersion("infercnv")))'],capture_output=True,text=True)
    (R/'public/infercnv_environment.json').write_text(json.dumps({'returncode':chk.returncode,'version':chk.stdout,'stderr':chk.stderr},indent=2))
    if chk.returncode:raise RuntimeError('infercnv load failed')
    assert (R/'PREPARED.json').exists()
    d=pd.read_csv(R/'private/units_index.tsv',sep='\t');d=d[d.ready].sort_values(['cohort','unit'])
    first=d[d.cohort=='GSE242230'].iloc[0]['path']
    result=[run_unit(first)]
    (R/'private/run_progress.json').write_text(json.dumps(result,indent=2))
    if result[0]['status']=='FAILED':raise RuntimeError('Pilot failed; other units not launched')
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:
        fs=[ex.submit(run_unit,p) for p in d.path if p!=first]
        for f in concurrent.futures.as_completed(fs):
            result.append(f.result());(R/'private/run_progress.json').write_text(json.dumps(result,indent=2))
            print('FINISHED_UNITS',len(result),'OF',len(d),flush=True)
    (R/'INFERENCE_BATCH_FINISHED.json').write_text(json.dumps({'total':len(result),'failed':sum(r['status']=='FAILED' for r in result)}))
