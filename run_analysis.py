from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import os,sys,subprocess,shutil,argparse
root=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('--start-at');a=ap.parse_args()
r=os.environ.get('RSCRIPT',shutil.which('Rscript') or '')
if not r or not os.environ.get('NHATS_RAW_DIR'):sys.exit('Set RSCRIPT and NHATS_RAW_DIR before running.')
for d in ['baseline/data_derived','baseline/tables','data','tables','figures','logs','documentation']:(root/d).mkdir(parents=True,exist_ok=True)
steps=['baseline_02_build.R','baseline_04_main_v2.R','01_scoped_data.R','02_models_repair.R','03_MI_repair.R','04_risks_scoped.R','05_ipcw_scoped.R','06_descriptives_scoped.R','07_figures_scoped.R','08_flow_scoped.R','10_MI_diagnostics.R','11_check_MI_pool.R','12_integration_models.R','13_integration_figures.R']
if a.start_at:
 if a.start_at not in steps:sys.exit('Unknown resume step')
 steps=steps[steps.index(a.start_at):]
def run(s,extra=(),tag=''):
 print(s,tag,flush=True)
 with (root/'logs'/(s+tag+'.log')).open('wb') as f:
  res=subprocess.run([r,'--vanilla',str(root/'code'/s),str(root),*map(str,extra)],stdout=f,stderr=subprocess.STDOUT)
 if res.returncode:raise RuntimeError('Failed: '+s+tag)
workers=max(1,min(8,int(os.environ.get('MI_WORKERS','4'))))
for s in steps:
 if s=='03_MI_repair.R' and workers>1:
  bounds=[]
  for i in range(workers):bounds.append((1+i*20//workers,(i+1)*20//workers))
  with ThreadPoolExecutor(max_workers=workers) as ex:
   futures=[ex.submit(run,'03_MI_chunk.R',b,'_'+str(b[0])) for b in bounds]
   for future in futures:future.result()
  run('03_MI_combine.R')
 else:run(s)
print('ALL_ANALYSIS_STEPS_COMPLETED',flush=True)
