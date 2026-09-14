#!/usr/bin/env python3
"""Fill the original result1 template and independently audit every result cell.
The delivered workbook was created with artifact_tool. A standard-library
OOXML fallback is provided for a local Python without that optional package;
that fallback was not used in the delivered workbook or runtime verification.
"""
from __future__ import annotations
import argparse,hashlib,json,math,zipfile
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np
NS='http://schemas.openxmlformats.org/spreadsheetml/2006/main'
Q=lambda name:'{'+NS+'}'+name
COLS=[chr(65+i) for i in range(22)]

def rounded(a):
    return [[float(f'{v:.4f}') for v in row] for row in a]

def export_excel(template,raw,out,workbook=None):
    template,raw,out=map(Path,(template,raw,out));out.parent.mkdir(parents=True,exist_ok=True)
    if template.resolve()==out.resolve():raise ValueError('Never overwrite the original template')
    z=np.load(raw);assert z['temperature_degC'].shape==(1801,21)
    blocks=[('温度','temperature_degC','温度/℃'),('水分浓度','moisture_dry_basis','干基含水率/(kg/kg)')]
    try:
        from artifact_tool import Blob,SpreadsheetFile
    except ImportError:
        return export_stdlib_fallback(template,z,out,blocks)
    wb=workbook if workbook is not None else SpreadsheetFile.import_xlsx(Blob.load(str(template)))
    for name,key,unit in blocks:
        sh=wb.worksheets.get_item(name)
        header=['时间/s \\ 径向距离/cm\n'+unit]+[float(f'{r:.1f}') for r in z['radius_cm']]
        sh.get_range('A1:V1').values=[header]
        data=rounded(z[key][1:])
        sh.get_range('A2:V1801').values=[[i+1]+row for i,row in enumerate(data)]
        sh.get_range('B2:V1801').set_number_format('0.0000')
        sh.get_range('B1:V1').set_number_format('0.0')
        sh.get_range('A2:A1801').set_number_format('0')
        sh.get_range('A1:V1801').format.horizontal_alignment='center'
        sh.get_range('A1:V1801').format.vertical_alignment='center'
        sh.get_range('A1:V1801').format.row_height=18
        sh.get_range('B1:V1801').format.column_width=11
        sh.get_range('A1:A1801').format.column_width=29
        sh.get_range('A1:V1').format.row_height=38
        sh.get_range('A1').format.wrap_text=True
        sh.get_range('A1:V1').format.font.bold=True
        sh.freeze_panes.freeze_rows(1);sh.freeze_panes.freeze_columns(1)
    SpreadsheetFile.export_xlsx(wb).save(str(out))
    return {'engine':'artifact_tool','path':str(out)}

def export_stdlib_fallback(template,z,out,blocks):
    """Portable fallback, supplied but NOT executed in this delivery environment.
    Retains template ZIP parts, replaces the two result sheets and appends one
    numeric style. Uses public OOXML/ZIP formats, no spreadsheet dependencies.
    """
    ET.register_namespace('',NS)
    with zipfile.ZipFile(template) as zz: parts={n:zz.read(n) for n in zz.namelist()}
    styles=ET.fromstring(parts['xl/styles.xml']);nf=styles.find(Q('numFmts'))
    if nf is None:nf=ET.Element(Q('numFmts'),{'count':'0'});styles.insert(0,nf)
    ids=[int(e.attrib['numFmtId']) for e in nf];fmtid=max([163]+ids)+1
    ET.SubElement(nf,Q('numFmt'),{'numFmtId':str(fmtid),'formatCode':'0.0000'});nf.set('count',str(len(nf)))
    xfs=styles.find(Q('cellXfs'));sidx=len(xfs)
    ET.SubElement(xfs,Q('xf'),{'numFmtId':str(fmtid),'fontId':'0','fillId':'0','borderId':'0','xfId':'0','applyNumberFormat':'1'})
    xfs.set('count',str(len(xfs)));parts['xl/styles.xml']=ET.tostring(styles,encoding='utf-8',xml_declaration=True)
    wb=ET.fromstring(parts['xl/workbook.xml']);rels=ET.fromstring(parts['xl/_rels/workbook.xml.rels'])
    targets={e.attrib['Id']:e.attrib['Target'] for e in rels}
    sheets=wb.find(Q('sheets'))
    for name,key,unit in blocks:
        record=next(s for s in sheets if s.attrib['name']==name)
        rid=record.attrib['{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id']
        p=targets[rid];p=p.lstrip('/') if p.startswith('/') else 'xl/'+p
        root=ET.Element(Q('worksheet'));ET.SubElement(root,Q('dimension'),{'ref':'A1:V1801'})
        views=ET.SubElement(root,Q('sheetViews'));view=ET.SubElement(views,Q('sheetView'),{'workbookViewId':'0'})
        ET.SubElement(view,Q('pane'),{'xSplit':'1','ySplit':'1','topLeftCell':'B2','activePane':'bottomRight','state':'frozen'})
        cols=ET.SubElement(root,Q('cols'))
        ET.SubElement(cols,Q('col'),{'min':'1','max':'1','width':'29','customWidth':'1'})
        ET.SubElement(cols,Q('col'),{'min':'2','max':'22','width':'11','customWidth':'1'})
        sd=ET.SubElement(root,Q('sheetData'))
        values=[['时间/s \\ 径向距离/cm; '+unit]+z['radius_cm'].tolist()]+[[i+1]+v for i,v in enumerate(rounded(z[key][1:]))]
        for i,row in enumerate(values,1):
            rr=ET.SubElement(sd,Q('row'),{'r':str(i)})
            for j,val in enumerate(row):
                a={'r':COLS[j]+str(i)}
                if i>1 and j>0:a['s']=str(sidx)
                if isinstance(val,str):
                    a['t']='inlineStr';cc=ET.SubElement(rr,Q('c'),a);ss=ET.SubElement(cc,Q('is'));ET.SubElement(ss,Q('t')).text=val
                else:
                    cc=ET.SubElement(rr,Q('c'),a);ET.SubElement(cc,Q('v')).text=format(float(val),'.17g')
        parts[p]=ET.tostring(root,encoding='utf-8',xml_declaration=True)
    with zipfile.ZipFile(out,'w',zipfile.ZIP_DEFLATED) as zz:
        for name,blob in parts.items():zz.writestr(name,blob)
    return {'engine':'stdlib OOXML fallback (not tested in the delivery runtime)','path':str(out)}

def audit_excel(path,raw,table_dir=None):
    """Independent ZIP/XML reader: does NOT trust the exporting library."""
    path,raw=Path(path),Path(raw);z=np.load(raw);rep={'file':str(path),'checks':{},'sheets':{}}
    with zipfile.ZipFile(path) as zz:
        assert zz.testzip() is None
        strings=[]
        if 'xl/sharedStrings.xml' in zz.namelist():
            strings=[''.join(e.text or '' for e in si.iter(Q('t'))) for si in ET.fromstring(zz.read('xl/sharedStrings.xml'))]
        wb=ET.fromstring(zz.read('xl/workbook.xml'));sheets=wb.find(Q('sheets'))
        names=[s.attrib['name'] for s in sheets];assert names==['温度','水分浓度'],names
        rels=ET.fromstring(zz.read('xl/_rels/workbook.xml.rels'));targets={r.attrib['Id']:r.attrib['Target'] for r in rels}
        st=ET.fromstring(zz.read('xl/styles.xml'));formats={int(e.attrib['numFmtId']):e.attrib['formatCode'] for e in st.find(Q('numFmts'))}
        styles=list(st.find(Q('cellXfs')))
        for s,key,tab in zip(sheets,['temperature_degC','moisture_dry_basis'],['temperature','moisture']):
            target=targets[s.attrib['{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id']]
            target=target.lstrip('/') if target.startswith('/') else 'xl/'+target
            root=ET.fromstring(zz.read(target));values={};cell_styles={}
            for c in root.findall('.//'+Q('sheetData')+'/'+Q('row')+'/'+Q('c')):
                v=c.find(Q('v'));typ=c.get('t','n');addr=c.attrib['r']
                if typ=='inlineStr':val=''.join(e.text or '' for e in c.iter(Q('t')))
                elif typ=='s':val=strings[int(v.text)]
                elif typ=='str':val=v.text if v is not None else ''
                elif v is None:val=None
                else:
                    assert typ in ('n',''),(addr,typ)
                    val=float(v.text)
                values[addr]=val;cell_styles[addr]=int(c.get('s','0'))
                assert c.find(Q('f')) is None,'No formulas expected in result arrays'
            assert '时间/s' in values['A1'] and '距离/cm' in values['A1']
            hdr=np.array([values[COLS[j]+'1'] for j in range(1,22)]);assert np.allclose(hdr,np.arange(21)/10,atol=1e-14)
            tt=np.array([values[f'A{i}'] for i in range(2,1802)]);assert np.array_equal(tt,np.arange(1,1801))
            a=np.array([[values[COLS[j]+str(i)] for j in range(1,22)] for i in range(2,1802)],float)
            assert a.shape==(1800,21) and np.isfinite(a).all()
            assert np.array_equal(a,np.array(rounded(z[key][1:]))),'Excel differs from declared four-decimal rounding'
            for i in range(2,1802):
                for j in range(1,22):
                    xf=styles[cell_styles[COLS[j]+str(i)]]
                    assert formats[int(xf.attrib['numFmtId'])]=='0.0000'
            assert not any(isinstance(v,str) and ('…' in v or '...' in v) for v in values.values())
            assert max(int(''.join(filter(str.isdigit,k))) for k in values)==1801
            if table_dir is not None:
                table=np.loadtxt(Path(table_dir)/f'table_{tab}.csv',delimiter=',',skiprows=1)
                assert np.array_equal(table[:,1:],a[table[:,0].astype(int)-1][:,::5])
            rep['sheets'][s.attrib['name']]={'time_count':1800,'radius_count':21,'result_count':37800,
              't_first_last':[int(tt[0]),int(tt[-1])],'r_first_last_cm':[float(hdr[0]),float(hdr[-1])],
              'first_center_surface':a[0,[0,-1]].tolist(),'last_center_surface':a[-1,[0,-1]].tolist(),
              'max_absolute_rounding_difference':float(np.max(np.abs(a-z[key][1:]))),
              'all_result_values_numeric_finite':True,'all_result_number_formats':'0.0000',
              'no_gaps_or_ellipses':True,'paper_tables_match':table_dir is not None}
    rep['checks']={'all_passed':True,'total_result_values':75600,'sheet_names_exact':True,'ZIP_CRC_passed':True}
    rep['sha256']=hashlib.sha256(path.read_bytes()).hexdigest();return rep

def main():
    ap=argparse.ArgumentParser();b=Path(__file__).resolve().parents[1]
    ap.add_argument('--template',type=Path,default=b/'input'/'result1_template.xlsx')
    ap.add_argument('--raw',type=Path,default=b/'output'/'q1_unrounded.npz')
    ap.add_argument('--out',type=Path,default=b/'output'/'result1.xlsx')
    args=ap.parse_args();info=export_excel(args.template,args.raw,args.out)
    audit=audit_excel(args.out,args.raw,args.out.parent)
    audit['export']=info;(args.out.parent/'excel_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(audit,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
