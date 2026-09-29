from pathlib import Path
import subprocess, tarfile, io, shutil, json, html, hashlib

BASE=Path('D:/CodexData/visualizations/2026/09/25/01a0d8cd-4545-75b0-ba43-922588deea1a')
OUT=BASE/'LYPLA1_spatial_selection_all'
OUT.mkdir(exist_ok=True)
ROOT=Path('C:/Users/Administrator/Documents/New project 5')
REMOTE='/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/results/collaborative/'
RUNS=[('COAD','20260928T084322Z_lypla1_spatial_v1/public'),('COAD','20260928T094500Z_lypla1_pathology_v1/public'),('COAD','20260928T103900Z_lypla1_stereoseq_v1/final_delivery'),('PDAC','20260928T094758Z_geomx_bell2025_lypla1_v1/public'),('PDAC','20260928T102605Z_gse235315_lypla1_spatial_v1/public'),('PDAC','20260928T104922Z_moncada_lypla1_slice_v1/public')]
items=[]
def add(cancer,group,p,note,source=None):
    dest=OUT/'images'/cancer/group/p.name
    dest.parent.mkdir(parents=True,exist_ok=True)
    if p.resolve()!=dest.resolve(): shutil.copy2(p,dest)
    items.append(dict(cancer=cancer,group=group,name=p.stem,note=note,src=dest.relative_to(OUT).as_posix(),source=source or str(p),sha256=hashlib.sha256(dest.read_bytes()).hexdigest()))

notes={'HER2ST':'作者病理区域；乳腺腺体不等于纯正常上皮','CTA':'作者癌区/免疫区标注；未标注区不是正常上皮','Wu':'空间表达及细胞背景，按原图注释解释','GSE210616':'癌区标注限制见原分析；不以表达亮度推断癌区','HD3':'高分辨率；缺少可核验癌区/正常导管标签','normal_duct_audit':'正常导管可获取性审计图，非癌区富集证明'}
for group,note in notes.items():
    d=BASE/('LYPLA1_'+group+('_spatial' if group!='normal_duct_audit' else ''))
    for p in sorted(d.glob('*.png')): add('BRCA',group,p,note)
repo=ROOT/'camp-ppt-unified-20260928'
for p in sorted((repo/'results/BRCA/07_INTEGRATION/20260928T150000Z_ppt_redraw_v1/server_figures').glob('spatial_PRAD*.png')):
    add('PRAD','Erickson2022',p,'7 张切片来自同一患者；作者病理共识标签，5 张可作癌区—良性腺体比较')
prad=ROOT/'camp-prad-analysis/results/PRAD/06_EXTERNAL'
for run,group in [('20260928T093400Z_lypla1_spatial_v1','Hirz2023'),('20260928T101700Z_lypla1_erickson_v1','Erickson2022_summary')]:
    for p in sorted((prad/run).glob('*.png')):add('PRAD',group,p,'Hirz：作者 RCTD 推断标签；Erickson：作者病理标签。汇总图不另算切片。')

def render():
    counts={c:sum(x['cancer']==c for x in items) for c in ['BRCA','COAD','PDAC','PRAD']}
    for c in counts:
        selected=[x for x in items if x['cancer']==c]
        for i,x in enumerate(selected,1):x['id']=f'{c}-{i:02d}'
    payload=json.dumps(items,ensure_ascii=False).replace('</','<\\/')
    template='''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>LYPLA1 空间图全量选片</title><style>
body{margin:0;background:#f3f6f5;color:#173d32;font-family:"Microsoft YaHei",sans-serif}header{position:sticky;top:0;background:#fff;padding:18px 28px;z-index:2;border-bottom:1px solid #d5e2db}h1{font-size:25px;margin:0 0 8px}p{margin:8px 0;color:#60746c}button,input{padding:9px 15px;border:1px solid #b8cfc2;background:white;border-radius:6px;margin:4px;color:#174535;cursor:pointer}button.active{background:#00583e;color:white}main{padding:20px;display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:20px}.card{background:white;border-radius:10px;padding:15px;border:1px solid #dce6df}.card img{width:100%;height:285px;object-fit:contain;cursor:zoom-in}.name{font-weight:bold;font-size:17px}.note{font-size:12px;color:#6b756f}.picked{outline:3px solid #91b94c}dialog{width:94vw;max-width:none;height:90vh;padding:10px;border:0;border-radius:8px}dialog img{width:100%;height:83%;object-fit:contain}dialog::backdrop{background:#000b}.toolbar{display:flex;justify-content:space-between;align-items:center}textarea{width:95%;height:100px}a{color:#00734f}@media(max-width:900px){main{grid-template-columns:1fr}} </style>
<header><h1>LYPLA1 · 四癌空间图选片库</h1><p>已有结果原图完整列出；点图放大，勾选后复制编号。汇总图、不同分辨率和重复排版不等于独立切片。</p><div id="tabs"></div><input id="search" placeholder="搜索切片 / 数据集 / 编号" oninput="draw()"><button onclick="copySelection()">复制所选编号</button><span id="count"></span></header><main id="grid"></main><dialog id="viewer"><div class="toolbar"><b id="caption"></b><button onclick="document.getElementById('viewer').close()">关闭 ×</button></div><img id="large"><p id="detail"></p><a id="original" target="_blank">打开原图</a></dialog>
<script>const data=PAYLOAD;let cancer='BRCA';const chosen=new Set(JSON.parse(localStorage.getItem('lypla1_spatial_choices')||'[]'));const esc=s=>s.replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));function draw(){document.getElementById('tabs').innerHTML=['BRCA','COAD','PDAC','PRAD','全部'].map(c=>`<button class="${c===cancer?'active':''}" onclick="cancer='${c}';draw()">${c} (${c==='全部'?data.length:data.filter(x=>x.cancer===c).length})</button>`).join('');const q=document.getElementById('search').value.toLowerCase();const list=data.filter(x=>(cancer==='全部'||x.cancer===cancer)&&JSON.stringify(x).toLowerCase().includes(q));document.getElementById('grid').innerHTML=list.map(x=>`<article class="card ${chosen.has(x.id)?'picked':''}"><div class="toolbar"><span class="name">${esc(x.id)} · ${esc(x.name)}</span><label><input type="checkbox" ${chosen.has(x.id)?'checked':''} onchange="pick('${x.id}',this.checked)">选用</label></div><p>${esc(x.group)}</p><img loading="lazy" src="${x.src}" onclick="show('${x.id}')"><p class="note">${esc(x.note)}</p></article>`).join('');document.getElementById('count').textContent=`当前 ${list.length} 幅 · 已选 ${chosen.size} 幅`;}function pick(id,v){v?chosen.add(id):chosen.delete(id);localStorage.setItem('lypla1_spatial_choices',JSON.stringify([...chosen]));draw()}function show(id){const x=data.find(x=>x.id===id);document.getElementById('caption').textContent=x.id+' · '+x.group+' · '+x.name;document.getElementById('large').src=x.src;document.getElementById('detail').textContent=x.note;document.getElementById('original').href=x.src;document.getElementById('viewer').showModal()}async function copySelection(){const text=data.filter(x=>chosen.has(x.id)).map(x=>x.id+' '+x.group+' '+x.name).join('\\n');try{await navigator.clipboard.writeText(text);alert('已复制 '+chosen.size+' 幅图的编号')}catch(e){prompt('复制这些编号',text)}}draw();</script></html>'''
    (OUT/'index.html').write_text(template.replace('PAYLOAD',payload),encoding='utf-8')
    (OUT/'inventory.json').write_text(json.dumps(items,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(counts),flush=True)

render()
import concurrent.futures
def fetch(entry):
    c,run=entry; directory=REMOTE+c+'/B/'+run
    cmd="cd "+directory+" && find . -maxdepth 3 -type f \\( -name '*.png' -o -name 'README_CN.md' \\) -print0 | tar --null -T - -cf -"
    for attempt in range(2):
        proc=subprocess.run(['ssh','-o','BatchMode=yes','-o','ConnectTimeout=10','server165',cmd],capture_output=True,timeout=110)
        if proc.returncode==0:break
    if proc.returncode: return c,run,[],proc.stderr.decode(errors='replace')
    group=run.split('/')[0]; dest=OUT/'remote'/c/group;dest.mkdir(parents=True,exist_ok=True);paths=[]
    with tarfile.open(fileobj=io.BytesIO(proc.stdout),mode='r:') as tf:
        for member in tf.getmembers():
            if not member.isfile():continue
            p=dest/Path(member.name).name
            if p.suffix.lower()!='.png' and p.name!='README_CN.md':continue
            p.write_bytes(tf.extractfile(member).read())
            if p.suffix=='.png':paths.append((p,member.name))
    return c,run,paths,None
errors=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
    for c,run,paths,error in pool.map(fetch,RUNS):
        if error: errors.append(dict(cancer=c,run=run,error=error));continue
        group=run.split('/')[0]
        note='已有空间分析原图；注释类型及证据范围见该数据集原报告。'
        if 'geomx' in group:note='GeoMx 区域/患者表达比较汇总，不是完整切片表达地图。'
        if c=='COAD':note='CRC（含结肠与直肠），不等于纯 COAD。'
        if '084322' in group:note+=' 本套未取得逐点病理癌区标签，仅展示表达位置。'
        if 'pathology' in group:note+=' 作者病理区域标签；完整保留非肿瘤上皮及反向结果。'
        if 'stereoseq' in group:note+=' Stereo-seq；复用作者癌相关上皮/正常上皮注释，约 50×50 微米空间单元。'
        if 'moncada' in group:note='PDAC-A ST1：同一患者一张切片；中右点位对应，H&E 尚未完成像素级配准。LYPLA1 检出稀疏。'
        for p,member in sorted(paths):add(c,group,p,note,REMOTE+c+'/B/'+run+'/'+member.removeprefix('./'))
        render()
(OUT/'fetch_errors.json').write_text(json.dumps(errors,ensure_ascii=False,indent=2),encoding='utf-8')
print('DONE',len(items),'errors',len(errors),flush=True)
