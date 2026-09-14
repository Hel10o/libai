// Run a copy of this file in an OS temporary directory whose node_modules
// symlink points to CODEX_PRIMARY_RUNTIME_NODE_MODULES, using CODEX_PRIMARY_RUNTIME_NODE.
// All XLSX authoring uses the documented @oai/artifact-tool API.
import fs from 'node:fs/promises';
import path from 'node:path';
import { FileBlob, SpreadsheetFile, Workbook } from '@oai/artifact-tool';

const args=Object.fromEntries(process.argv.slice(2).map(s=>{const p=s.indexOf('=');return [s.slice(0,p).replace(/^--/,''),s.slice(p+1)];}));
const root=args.root||process.cwd();
const mode=args.mode||'probe';
function report(stage, extra={}) { console.log(JSON.stringify({stage,elapsed_s:(Date.now()-start)/1000,rss_MB:process.memoryUsage().rss/2**20,...extra})); }
const start=Date.now();
const font=args.font||'Arial';
async function templateFile(name) {
  const bundled=name.replace('.xlsx','_template.xlsx');
  for(const candidate of [path.join(root,'inputs',bundled),path.join(root,'q4_complete_delivery/v1/inputs',bundled),path.join(root,'A题/附件/附件3',name)]) {
    try { await fs.access(candidate); return candidate; } catch {}
  }
  throw new Error(`Missing supplied template ${bundled}`);
}
// Decimal ROUND_HALF_UP of the canonical decimal representation of each
// double. Values far from a half-way case use an equivalent fast path.
function halfUp4(value) {
  if(value===null) return null;
  if(!Number.isFinite(value)) throw new Error('Non-finite in-domain value');
  const sign=value<0?-1:1,v=Math.abs(value),scaled=v*1e4;
  if(Math.abs(scaled-Math.floor(scaled)-0.5)>1e-8) return sign*Math.floor(scaled+0.5)/1e4;
  const [mant,expText='0']=v.toString().split('e');
  const [whole,fraction='']=mant.split('.');
  const integer=BigInt(whole+fraction);
  const places=fraction.length-Number(expText)-4;
  const rounded=places<=0?integer*10n**BigInt(-places):(integer+5n*10n**BigInt(places-1))/(10n**BigInt(places));
  return sign*Number(rounded)/1e4;
}

async function saveCheck(wb, outfile, label, sampleSheets=[]) {
  wb.recalculate();
  const checks=[];
  for(const [sn,rg] of sampleSheets) {
    const chk=await wb.inspect({kind:'table',range:`'${sn}'!${rg}`,include:'values,formulas',tableMaxRows:3,tableMaxCols:24,maxChars:5000});
    checks.push(JSON.parse(chk.ndjson.split('\n').filter(Boolean)[0]));
  }
  const errors=await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A|#NUM!|#NULL!|#SPILL!|#CALC!',options:{useRegex:true,maxResults:10},maxChars:2000});
  report('artifact_checked',{label,checks,errors:errors.ndjson});
  await fs.mkdir(path.dirname(outfile),{recursive:true});
  const output=await SpreadsheetFile.exportXlsx(wb);
  await output.save(outfile);
  report('artifact_saved',{label,outfile});
}

function formatSheet(sh,firstRow,lastRow,lastCol,header) {
  sh.getRange('A1:F5').clear({applyTo:'all'});
  sh.getRangeByIndexes(0,0,1,header.length).values=[header];
  sh.getRangeByIndexes(0,0,1,header.length).format={font:{name:font,size:10,bold:true},verticalAlignment:'center',horizontalAlignment:'center',wrapText:true,rowHeight:31};
  sh.getRange('A1').format.columnWidth=28;
  sh.getRangeByIndexes(0,1,1,header.length-1).format.columnWidth=10;
  sh.getRangeByIndexes(firstRow-1,0,lastRow-firstRow+1,lastCol).format={font:{name:font,size:10},verticalAlignment:'center',horizontalAlignment:'right'};
  sh.getRangeByIndexes(firstRow-1,0,lastRow-firstRow+1,1).setNumberFormat('0.########');
  sh.getRangeByIndexes(firstRow-1,1,lastRow-firstRow+1,lastCol-1).setNumberFormat('0.0000');
  sh.freezePanes.freezeRows(1);
  sh.freezePanes.freezeColumns(1);
}

if(mode==='q2chunk') {
  const dir=args.inputDir;
  const meta=JSON.parse(await fs.readFile(path.join(dir,'q2_meta.json'),'utf8'));
  const times=await fs.readFile(path.join(dir,'q2_time.f64'));
  const profiles=await fs.readFile(path.join(dir,'q2_profile.f64'));
  const lo=Number(args.lo),hi=Number(args.hi);
  if(!(lo>=1&&hi<=meta.shape[0]&&hi>lo)) throw new Error('Invalid source index interval');
  const wb=await SpreadsheetFile.importXlsx(await FileBlob.load(await templateFile('result2.xlsx')));
  for(const [sn,component] of [['温度',0],['水分浓度',1]]) {
    const sh=wb.worksheets.getItem(sn);
    formatSheet(sh,lo+1,hi,22,['时间\\到药材中心的距离',...meta.radius_cm]);
    for(let b=lo;b<hi;b+=1000) {
      const e=Math.min(hi,b+1000),rows=[];
      for(let i=b;i<e;i++) {
        const row=[times.readDoubleLE(i*8)];
        for(let j=0;j<21;j++) row.push(halfUp4(profiles.readDoubleLE(((i*21+j)*2+component)*8)));
        rows.push(row);
      }
      sh.getRangeByIndexes(b,0,e-b,22).values=rows;
    }
  }
  report('q2_chunk_populated',{lo,hi});
  if(args.preview==='true') for(const sn of ['温度','水分浓度']) await inspectAndRender(wb,sn,`result2_${sn}_head`,'A1:V8');
  await saveCheck(wb,args.out,`q2_${lo}_${hi}`,[['温度',`A${lo+1}:V${Math.min(lo+2,hi)}`],['水分浓度',`A${hi}:V${hi}`]]);
}

if(mode==='q4') {
  const data=JSON.parse(await fs.readFile(args.input,'utf8'));
  const indexes=data.time_s.map((t,i)=>t>0?i:-1).filter(i=>i>=0);
  const wb=await SpreadsheetFile.importXlsx(await FileBlob.load(await templateFile('result4.xlsx')));
  const sh=wb.worksheets.getItem('Sheet1');
  formatSheet(sh,2,indexes.length+1,23,['时间\\到药材中心的距离',...Array.from({length:21},(_,j)=>j/10),'药材表面']);
  const rows=indexes.map(i=>[data.time_s[i],...data.C_fixed[i].map(halfUp4),halfUp4(data.C_surface[i])]);
  sh.getRangeByIndexes(1,0,rows.length,23).values=rows;
  sh.getRange('W1').format.columnWidth=13;
  sh.getRange('Y1').values=[['实际半径R(t) / cm']];
  sh.getRange('Y1').format={font:{name:font,size:10,bold:true},wrapText:true,columnWidth:18,verticalAlignment:'center',horizontalAlignment:'center'};
  sh.getRangeByIndexes(1,24,rows.length,1).values=indexes.map(i=>[data.radius_m[i]*100]);
  sh.getRangeByIndexes(1,24,rows.length,1).setNumberFormat('0.000000');
  sh.getRange('AA1').values=[['输出说明']];
  sh.getRange('AA2:AA6').values=[['时间单位：秒；固定径向坐标单位：厘米。'],['水分浓度为干基含水率，单位kg水/kg干物质。'],['固定径向点超出当前半径时留空，不代表零。'],['表面值取r=R(t)；与固定点重合时两列数值一致。'],['含水率按十进制四位ROUND_HALF_UP舍入；未舍入轨迹另附，严格终判使用未舍入全域最大值。']];
  sh.getRange('AA1:AA6').format={font:{name:font,size:10},columnWidth:94};
  await inspectAndRender(wb,'Sheet1','result4_head','A1:Y8');
  await inspectAndRender(wb,'Sheet1','result4_tail',`A${Math.max(2,rows.length-5)}:Y${rows.length+1}`);
  await saveCheck(wb,args.out,'q4',[['Sheet1','A1:Y3'],['Sheet1',`A${rows.length+1}:Y${rows.length+1}`]]);
}

async function inspectAndRender(wb,sheetName,tag,range='A1:W8') {
  const inspect=await wb.inspect({kind:'table',range:`'${sheetName}'!${range}`,include:'values,formulas',tableMaxRows:8,tableMaxCols:24,maxChars:5000});
  console.log(inspect.ndjson);
  const img=await wb.render({sheetName,range,scale:1.3,format:'png'});
  await fs.writeFile(path.join(args.previewDir||process.env.DRYING_WORKBOOK_TEMP||path.join(process.cwd(),'workbook_temp'),`${tag}.png`),new Uint8Array(await img.arrayBuffer()));
}

if(mode==='preview') {
  const wb=await SpreadsheetFile.importXlsx(await FileBlob.load(args.input));
  for(const sheetName of (args.sheets||'Sheet1').split(',')) {
    await inspectAndRender(wb,sheetName,`${args.tag||'preview'}_${sheetName}`,args.range||'A1:W8');
  }
}

if(mode==='probe') {
  const n=Number(args.rows||20000);
  for(const name of ['result2.xlsx','result4.xlsx']) {
    const original=await SpreadsheetFile.importXlsx(await FileBlob.load(await templateFile(name)));
    report('template_import',{name});
    const names=name==='result2.xlsx'?['温度','水分浓度']:['Sheet1'];
    for(const sn of names) await inspectAndRender(original,sn,`template_${name}_${sn}`,'A1:F5');
  }
  const wb=Workbook.create();
  for(const sn of ['capacity_probe_T','capacity_probe_C']) {
    const sh=wb.worksheets.add(sn);
    sh.getRange('A1:V1').values=[['BENCHMARK ONLY — synthetic data',...Array.from({length:21},(_,i)=>i/10)]];
    for(let i=0;i<n;i+=2000) {
      const count=Math.min(2000,n-i);
      const rows=Array.from({length:count},(_,j)=>[i+j+1,...Array.from({length:21},(_,c)=>(i+j+c)/10000)]);
      sh.getRangeByIndexes(i+1,0,count,22).values=rows;
    }
    sh.getRange(`B2:V${n+1}`).setNumberFormat('0.0000');
    sh.freezePanes.freezeRows(1);
    report('probe_populated',{sheet:sn,rows:n});
  }
  wb.recalculate(); report('probe_recalculate');
  const out=await SpreadsheetFile.exportXlsx(wb); report('probe_export');
  await out.save(path.join(process.env.DRYING_WORKBOOK_TEMP||path.join(process.cwd(),'workbook_temp'),`capacity_probe_${n}.xlsx`));
  report('probe_saved');
}
