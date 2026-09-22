import fs from 'node:fs/promises';
import path from 'node:path';
import { Workbook, SpreadsheetFile } from '@oai/artifact-tool';
const out=path.resolve(process.argv[2]);
const payload=JSON.parse(await fs.readFile(path.join(out,'repro','workbook_data.json'),'utf8'));
const wb=Workbook.create();
function letter(n){let s='';while(n){n--;s=String.fromCharCode(65+n%26)+s;n=Math.floor(n/26);}return s;}
const labels={gene:'基因',human_gene_id:'稳定基因ID',current_direct:'当前直接池',current_conditional:'当前条件池',current_unresolved:'当前未解决',historical_retained:'历史保留',identity_status:'身份状态',RNA_effect:'RNA配对均值差',RNA_p_value:'RNA P值',RNA_q_value:'RNA q值',RNA_status:'RNA状态',Lee_display_top:'Lee来源',Uhlitz_display_top:'Uhlitz来源',Lee_shared_top:'Lee共同类别来源',Uhlitz_shared_top:'Uhlitz共同类别来源',descriptive_stable_concordance:'宽类别描述性稳定一致',counterevidence_or_limit_cn:'反证与限制',metabolite_name:'代谢物',metabolite_key:'代谢物稳定键',effect:'效应',p_value:'P值',q_value:'q值',status:'状态',reason:'原因',n:'实际样本数',cohort:'队列',run_id:'运行版本',analysis_type:'分析类型',celltype:'细胞类别',mean_expression:'供者等权平均表达',mean_detection_fraction:'供者平均检出比例',source_status:'来源状态',display_top:'可解释来源',mechanical_top:'机械最高类别'};
let index=0;const qa=[];
for(const [name,data] of Object.entries(payload)){
 if(!data.columns.length)continue;
 const s=wb.worksheets.add(name);s.showGridLines=false;
 const rows=data.rows.map(row=>row.map((v,j)=>{
  if(v===null)return 'NA';
  if(typeof v==='string'&&/^-?(?:\d+\.?\d*|\.\d+)(?:e[+-]?\d+)?$/i.test(v)&&!/(?:_id|sha|gene|symbol|seed)/i.test(data.columns[j]))return Number(v);
  if(typeof v==='string'&&v.startsWith('='))return "'"+v;
  return v;
 }));
 const heads=data.columns.map(c=>labels[c]||c);
 const range=s.getRangeByIndexes(0,0,rows.length+1,heads.length);range.values=[heads,...rows];
 range.format.font={name:'Microsoft YaHei',size:10};range.format.rowHeight=22;range.format.columnWidth=20;
 const header=s.getRangeByIndexes(0,0,1,heads.length);header.format={fill:'#087F8C',font:{name:'Microsoft YaHei',bold:true,color:'#FFFFFF'},wrapText:true,rowHeight:46};
 for(let j=0;j<heads.length;j++){
  const col=s.getRangeByIndexes(0,j,rows.length+1,1);const field=data.columns[j];
  col.format.columnWidth=Math.min(48,Math.max(16,heads[j].length*1.25));
  if(/(?:p_value|q_value|_P$|_q$)/i.test(field))s.getRangeByIndexes(1,j,Math.max(rows.length,1),1).setNumberFormat('0.000E+00');
  else if(/(?:effect|difference|frequency|fraction|ci_|detection|expression)/i.test(field))s.getRangeByIndexes(1,j,Math.max(rows.length,1),1).setNumberFormat('0.0000');
  else if(/^(n|n_|planned|evaluable)|pairs$/.test(field))s.getRangeByIndexes(1,j,Math.max(rows.length,1),1).setNumberFormat('#,##0');
 }
 if(name==='先读我'){s.getRangeByIndexes(0,0,rows.length+1,1).format.columnWidth=100;s.getRangeByIndexes(1,0,rows.length,1).format.wrapText=true;}
 else{const t=s.tables.add(`A1:${letter(heads.length)}${rows.length+1}`,true,'COADTable'+(++index));t.showFilterButton=true;s.freezePanes.freezeRows(1);}
 qa.push({sheet:name,rows:rows.length,columns:heads.length});
}
wb.recalculate();
console.log((await wb.inspect({kind:'table',range:'基因易读表!A1:H5',include:'values',tableMaxRows:5,tableMaxCols:8,maxChars:1800})).ndjson);
await fs.mkdir(path.join(out,'workbook_previews'),{recursive:true});
for(const x of qa){const blob=await wb.render({sheetName:x.sheet,range:`A1:${x.columns===1?'A':'F'}${Math.min(x.rows+1,7)}`,scale:1.3});await fs.writeFile(path.join(out,'workbook_previews',x.sheet+'.png'),new Uint8Array(await blob.arrayBuffer()));}
const file=await SpreadsheetFile.exportXlsx(wb);await file.save(path.join(out,'COAD_完整结果册.xlsx'));
await fs.writeFile(path.join(out,'workbook_validation.json'),JSON.stringify({status:'EXPORTED',sheets:qa,numeric_types:true,filters:true,frozen_headers:true,visual_review:'PENDING'},null,2));
console.log('WORKBOOK_EXPORTED');
