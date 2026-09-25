"""Wait for this run's dependency setup, then infer and summarize on server."""
from pathlib import Path
import json,subprocess,time,traceback
R=Path(__file__).resolve().parent
try:
    start=time.monotonic()
    while True:
        log=(R/'setup_retry.log').read_text(errors='replace') if (R/'setup_retry.log').exists() else ''
        if 'SETUP_EXIT 0' in log and (R/'PREPARED.json').exists():break
        if 'SETUP_EXIT ' in log and 'SETUP_EXIT 0' not in log:raise RuntimeError('Dependency setup failed')
        if time.monotonic()-start>1800:raise TimeoutError('Dependency readiness timeout;inference not started')
        time.sleep(15)
    (R/'.running').write_text('INFERENCE')
    subprocess.run(['python3',str(R/'myeloid_infercnv_schedule.py')],cwd=R,check=True)
    subprocess.run(['python3',str(R/'myeloid_infercnv_summarize.py')],cwd=R,check=True)
    subprocess.run(['python3',str(R/'myeloid_infercnv_manifest.py')],cwd=R,check=True)
    (R/'PIPELINE_FINISHED.json').write_text(json.dumps({'status':'COMPLETED_EXECUTION','check_failed_units_and_unevaluable_tests':True}))
    (R/'.running').unlink()
except Exception:
    (R/'FAILED.txt').write_text(traceback.format_exc())
    raise
