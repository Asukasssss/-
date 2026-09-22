"""Independent synthetic fixtures; never written into scientific report outputs."""
import sys,unittest,tempfile,copy,importlib.util,json,hashlib,subprocess
from pathlib import Path
sys.dont_write_bytecode=True
R=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('renderer',R/'code/camp_results_report.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
import pandas as pd
def fixture():
 c=dict(cancer='SYNTHETIC',input_commit='0'*40,columns={'relation_id':'relation_key'},families={'metabolite_primary':'MP','metabolite_sensitivity':'MS','association_primary':'AP','association_sensitivity':'AS'},expected={'metabolites':2,'workpool':1,'relations':2,'genes':3,'current_genes':2},examples=[dict(gene='G1',relation_id='M1|G1',role='synthetic',reason='test',limitation='synthetic')],sc=dict(profile_sources=['profiles'],stability_sources=['stability'],partition='ALL',cohorts=['StudyA'],partition_field='partition',cohort_field='cohort',type_field='celltype',stability_top_field='top_lineage',detection_field='mean_detection_fraction'),legacy_source_display_columns=['old_top'])
 rows=[]
 for i in range(2):rows.append(dict(gene='G'+str(i+1),metabolite_key='M'+str(i+1),metabolite_name='met'+str(i+1),relation_id='M'+str(i+1)+'|G'+str(i+1),cohort='Synthetic',analysis_type='MP',effect=str(i),p_value='.01' if i==0 else '.8',q_value='.08' if i==0 else '.8',n='7',status='DONE'))
 d=dict(metabolites=pd.DataFrame(rows),association=pd.DataFrame([dict(x,analysis_type='AP') for x in rows]),genes=pd.DataFrame({'gene':['G1','G2','G3']}),reader=pd.DataFrame({'gene':['G1','G2','G3'],'old_top':['Type1']*3}),relations=pd.DataFrame([dict(relation_key=x['relation_id'],gene=x['gene'],metabolite_key=x['metabolite_key'],metabolite_name=x['metabolite_name']) for x in rows]),mapping=pd.DataFrame({'metabolite_key':['M1'],'mapping_state':['DIRECT_MAPPED']}),profiles=pd.DataFrame([dict(gene=g,cohort='StudyA',partition='ALL',celltype='Type1',status='DONE',n='3',effect='0',mean_detection_fraction='0') for g in ['G1','G2','G3']]),stability=pd.DataFrame([dict(gene=g,cohort='StudyA',partition='ALL',status='NOT_EVALUABLE',reason='low_detection',top_lineage='Type1',bootstrap_top_frequency='NA') for g in ['G1','G2','G3']]))
 return c,d
class DisplayTests(unittest.TestCase):
 def test_different_sizes_missing_optional_modules_and_one_study(self):
  c,d=fixture();z=m.prepare(c,d);v=m.validate(c,d,z)
  self.assertEqual(v['counts']['genes'],3);self.assertTrue(z['rna'].empty);self.assertTrue(z['sensitivity'].empty);self.assertNotIn('subtypes',d)
  self.assertEqual(z['source_status'].cohort.nunique(),1)
 def test_unevaluable_rank_never_becomes_localization(self):
  c,d=fixture();z=m.prepare(c,d)
  self.assertTrue(z['source_status'].display_top.eq('暂不可定位').all());self.assertNotIn('old_top',z['reader']);self.assertTrue(z['source_status'].mechanical_top.eq('Type1').all())
 def test_zero_effect_is_not_missing(self):
  c,d=fixture();self.assertEqual(m.tiers(d['metabolites'])[0],'仅P<0.05')
  t=d['metabolites'].copy();t.loc[0,'status']='NOT_EVALUABLE';self.assertEqual(m.tiers(t)[0],'不可评估')
 def test_duplicate_and_wrong_relationship_rejected(self):
  c,d=fixture();c['examples'][0]['relation_id']='M2|G2'
  with self.assertRaises(ValueError):m.prepare(c,d)
 def test_render_refuses_nonempty_output(self):
  c,d=fixture();z=m.prepare(c,d)
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp);(p/'sentinel').write_text('keep')
   with self.assertRaises(FileExistsError):m.build(c,d,z,{},pd.DataFrame(),p,'unused')
   self.assertEqual((p/'sentinel').read_text(),'keep')
 def test_synthetic_render_missing_rna_subtypes_mask_single_study(self):
  c,d=fixture();real=json.loads((R/'configs/BRCA_report_v1.yaml').read_text(encoding='utf-8'))
  c.update(font_candidates=real['font_candidates'],narratives={k:['合成测试问题','合成数据，仅用于程序测试。','不是任何癌种实际结果。'] for k in real['narratives']},design=[dict(label='合成样本',value=7)],assets=[],files={k:'SYNTHETIC/'+k for k in d},modules={'rna':False,'subtypes':False,'single_cell':True},subtypes=[])
  c['sc'].update(celltypes=['Type1'],labels={'Type1':'合成类别'},heatmap_page_rows=2)
  for key in ['metabolites','association']:
   d[key]['ci_lower']='-0.5';d[key]['ci_upper']='0.5'
  d['metabolites']['pairs_higher']='4';d['metabolites']['pairs_lower']='3';d['metabolites']['pairs_equal']='0'
  z=m.prepare(c,d)
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp);cfg=p/'synthetic.yaml';cfg.write_text(json.dumps(c),encoding='utf-8');out=p/'render'
   result=m.build(c,d,z,{},pd.DataFrame([{'source':'synthetic'}]),out,cfg)
   self.assertEqual(result['modules']['rna'],'NOT_EVALUABLE');self.assertEqual(result['modules']['subtypes'],'NOT_RUN')
   self.assertTrue((out/'SYNTHETIC_主报告.pdf').exists());self.assertFalse(any('subtype' in x.name for x in (out/'figures').glob('*')))
 def test_validate_only_does_not_create_output_or_change_tracked_inputs(self):
  tracked=subprocess.check_output(['git','-c','core.quotepath=false','ls-files'],cwd=R,encoding='utf-8').splitlines()
  before={p:hashlib.sha256((R/p).read_bytes()).hexdigest() for p in tracked if (R/p).is_file()}
  with tempfile.TemporaryDirectory() as tmp:
   out=Path(tmp)/'must_not_exist'
   subprocess.run([sys.executable,str(R/'code/camp_results_report.py'),'--config',str(R/'configs/BRCA_report_v1.yaml'),'--out',str(out),'--validate-only'],cwd=R,check=True,capture_output=True)
   self.assertFalse(out.exists())
  after={p:hashlib.sha256((R/p).read_bytes()).hexdigest() for p in before}
  self.assertEqual(before,after)
if __name__=='__main__':unittest.main()
