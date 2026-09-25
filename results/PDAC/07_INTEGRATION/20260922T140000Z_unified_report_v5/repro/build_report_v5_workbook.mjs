import fs from 'node:fs/promises';
import {Workbook,SpreadsheetFile} from '@oai/artifact-tool';
const out=process.argv[2],d=JSON.parse(await fs.readFile(out+'/workbook_data.json','utf8'));
const wb=Workbook.create();
const cn={gene:'基因',cohort:'队列',stable_gene_id:'稳定基因ID',in_current_pool:'当前直接池',history_only:'仅历史保留',RNA_effect:'RNA均值差',RNA_p:'RNA P',RNA_q:'RNA q',RNA_n:'RNA配对数',RNA_background:'RNA支持层级',metabolite_name:'代谢物',metabolite_key:'代谢物键',relation_id:'关系编号',effect:'效应（见尺度）',p_value:'P',q_value:'q',n:'有效单位数',celltype:'细胞类别',status:'状态',reason:'原因',mean_detection_fraction:'等权检出比例',effect_type:'效应尺度',normalization:'表达处理尺度',module:'模块',step:'阶段',planned:'计划数',evaluable:'可评估数',P_lt005:'P<0.05数',q_lt005:'q<0.05数',test_family:'检验家族'};
const col=n=>{let s='';for(n++;n;n=Math.floor((n-1)/26))s=String.fromCharCode(65+(n-1)%26)+s;return s;};
let t=0;
for(const x of d.sheets){
 const sh=wb.worksheets.add(x.name);sh.showGridLines=false;
 const end=col(x.headers.length-1),last=x.rows.length+1;
 sh.getRange(`A1:${end}${last}`).values=[x.headers.map(h=>cn[h]||h),...x.rows.map(r=>r.map(v=>typeof v==='string'&&/^[=+@]/.test(v)?"'"+v:v))];
 sh.getRange(`A1:${end}${last}`).format={font:{name:'Microsoft YaHei',size:10},rowHeight:40,verticalAlignment:'center'};
 sh.getRange(`A1:${end}1`).format={fill:'#2C526B',font:{color:'#FFFFFF',bold:true},wrapText:true,rowHeight:72,horizontalAlignment:'center'};
 sh.tables.add(`A1:${end}${last}`,true,'PDACV5_'+(++t));sh.freezePanes.freezeRows(1);sh.freezePanes.freezeColumns(1);
 const widths=[];
 x.headers.forEach((h,i)=>{const letter=col(i),numeric=x.numeric.includes(h);const width=h==='reason'||h.includes('limitation')?90:h.includes('name')||h.includes('path')||h==='metabolite_key'?50:h==='cohort'?34:numeric?18:28;widths.push(width);sh.getRange(`${letter}1:${letter}${last}`).format.columnWidth=width;
 if(!numeric)sh.getRange(`${letter}2:${letter}${last}`).format.wrapText=true;
 if(numeric)sh.getRange(`${letter}2:${letter}${last}`).setNumberFormat(/p_value|q_value|RNA_p|RNA_q/.test(h)?'0.000000':/fraction|frequency/.test(h)?'0.0%':/^(n|planned|evaluable)$|pairs|n_cells|P_lt|q_lt/.test(h)?'0':'0.0000');
 });
 x.rows.forEach((r,i)=>{let height=40;for(let j=0;j<r.length;j++){if(typeof r[j]!=='string')continue;const weight=[...r[j]].reduce((n,ch)=>n+(ch.charCodeAt(0)>255?2:1),0);height=Math.max(height,Math.ceil(weight/(widths[j]*.82))*15+12);}sh.getRange(`A${i+2}:${end}${i+2}`).format.rowHeight=Math.min(409,height);});
 console.log(JSON.stringify({sheet:x.name,rows:x.rows.length}));
}
wb.recalculate();
console.log((await wb.inspect({kind:'table',range:'基因易读表!A1:I4',include:'values',tableMaxRows:4,tableMaxCols:9,maxChars:1800})).ndjson);
console.log((await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#NUM!|#SPILL!',options:{useRegex:true,maxResults:5},maxChars:600})).ndjson);
await fs.mkdir(out+'/workbook_previews',{recursive:true});
for(let i=0;i<d.sheets.length;i++){
 const x=d.sheets[i],end=col(Math.min(x.headers.length,7)-1);
 const png=await wb.render({sheetName:x.name,range:`A1:${end}${Math.min(x.rows.length+1,7)}`,scale:1.2,format:'png'});
 await fs.writeFile(out+`/workbook_previews/sheet_${String(i+1).padStart(2,'0')}.png`,new Uint8Array(await png.arrayBuffer()));
}
const file=await SpreadsheetFile.exportXlsx(wb);await file.save(out+'/PDAC_完整结果册.xlsx');
console.log(JSON.stringify({exported:true,sheets:d.sheets.length}));
