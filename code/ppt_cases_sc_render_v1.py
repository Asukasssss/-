"""Reuse the accepted LYPLA1 renderer with explicit gene and page substitutions."""
from pathlib import Path
import sys,shutil
R=Path(sys.argv[1]);AS=R.parent/'20260929T071500Z_four_sc_v3'
for n in ['msyh.ttc','msyhbd.ttc','sysu_logo.png']:shutil.copy2(AS/n,R/n)
template=(R/'ppt_four_sc_render_v3.py').read_text()
for c,g,p in [('BRCA','ASNS',2),('COAD','UCKL1',9),('PDAC','SLC6A6',5),('PRAD','SLC6A6',6)]:
 code=template.replace("enumerate(['BRCA','COAD','PDAC','PRAD'],11)",f"enumerate(['{c}'],{p})").replace('LYPLA1',g)
 code=code.replace("f'{c}：从全细胞表达定位到同患者证据'",f"f'{{c}}：{g} 的全细胞定位与同患者证据'")
 code=code.replace("name=f'{page}_{c}_", "name=f'{page:02d}_{c}_")
 # Pair scope and conclusions remain supplied by the audited per-case metadata.
 exec(compile(code,str(R/'ppt_four_sc_render_v3.py'),'exec'),{'__name__':'__main__'})
