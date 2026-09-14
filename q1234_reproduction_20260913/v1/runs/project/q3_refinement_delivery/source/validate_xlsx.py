"""Strict OOXML readback validator for the frozen Q3 workbook.

This validator is read-only.  It checks the exact allowed Sheet1 rectangle,
header semantics, numeric types, four-decimal moisture format, row times, radii,
all moisture values against output/solution.npz, and the final strict endpoint.
"""
from __future__ import annotations
import argparse, json, re, zipfile
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np

NS={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
RID='{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id'
CELL_RE=re.compile(r'^([A-Z]+)([1-9][0-9]*)$')

def colnum(s):
    n=0
    for ch in s:n=n*26+ord(ch)-64
    return n

def read_book(path:Path):
    with zipfile.ZipFile(path) as z:
        if z.testzip(): raise ValueError('Corrupt XLSX')
        wb=ET.fromstring(z.read('xl/workbook.xml'))
        relroot=ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))
        rel={e.attrib['Id']:e.attrib['Target'] for e in relroot}
        sheets=wb.findall('s:sheets/s:sheet',NS)
        names=[s.attrib['name'] for s in sheets]
        if names!=['Sheet1']: raise ValueError(f'Unexpected sheets: {names}')
        target=rel[sheets[0].attrib[RID]]
        target=target[1:] if target.startswith('/') else 'xl/'+target
        root=ET.fromstring(z.read(target))
        strings=[]
        if 'xl/sharedStrings.xml' in z.namelist():
            sr=ET.fromstring(z.read('xl/sharedStrings.xml'))
            for si in sr.findall('s:si',NS):
                strings.append(''.join(t.text or '' for t in si.findall('.//s:t',NS)))
        # numFmt map by style id.
        styles=ET.fromstring(z.read('xl/styles.xml'))
        custom={int(e.attrib['numFmtId']):e.attrib['formatCode'] for e in styles.findall('s:numFmts/s:numFmt',NS)}
        builtin={0:'General',1:'0',2:'0.00',3:'#,##0',4:'#,##0.00'}
        xfs=styles.findall('s:cellXfs/s:xf',NS)
        style_fmt=[custom.get(int(x.attrib.get('numFmtId','0')),builtin.get(int(x.attrib.get('numFmtId','0')),str(x.attrib.get('numFmtId','0')))) for x in xfs]
        cells={}
        max_row=max_col=0
        for c in root.findall('.//s:sheetData/s:row/s:c',NS):
            ref=c.attrib['r'];m=CELL_RE.match(ref)
            if not m:raise ValueError(ref)
            col,row=m.group(1),int(m.group(2));max_row=max(max_row,row);max_col=max(max_col,colnum(col))
            if c.find('s:f',NS) is not None:raise ValueError(f'Formula not allowed: {ref}')
            typ=c.attrib.get('t','n');v=c.find('s:v',NS)
            if typ=='s': value=strings[int(v.text)]
            elif typ=='inlineStr': value=''.join(t.text or '' for t in c.findall('.//s:t',NS))
            elif v is None:value=None
            elif typ in ('str','e'):value=v.text
            else:value=float(v.text)
            sid=int(c.attrib.get('s','0'));fmt=style_fmt[sid] if sid<len(style_fmt) else None
            cells[ref]=(value,typ,fmt)
        return cells,max_row,max_col

def validate(xlsx:Path, solution:Path):
    cells,max_row,max_col=read_book(xlsx)
    z=np.load(solution);t=z['time_s'][1:];C=z['profile_TC'][1:,:,1];r=z['radius_cm']
    n=len(t)+1
    if (max_row,max_col)!=(n,22):raise ValueError(f'Exact extent must be A1:V{n}, got row={max_row}, col={max_col}')
    expected_refs=set()
    for row in range(1,n+1):
        for col in range(1,23):
            x=col;letters=''
            while x:
                x,rem=divmod(x-1,26);letters=chr(65+rem)+letters
            expected_refs.add(f'{letters}{row}')
    if set(cells)!=expected_refs:
        extra=sorted(set(cells)-expected_refs)[:10];missing=sorted(expected_refs-set(cells))[:10]
        raise ValueError(f'Unexpected cell set; extra={extra}, missing={missing}')
    if cells['A1'][0] != '时间/s\\到药材中心的距离/cm':raise ValueError('Wrong A1 title/unit')
    for j in range(21):
        ref=f'{chr(66+j)}1'
        if abs(float(cells[ref][0])-r[j])>1e-12:raise ValueError(f'Wrong radius {ref}')
    max_t=max_c=0.0;bad_format=[]
    for i in range(len(t)):
        row=i+2
        v,typ,_=cells[f'A{row}']
        if typ not in ('n',None):raise ValueError(f'Non-numeric time A{row}')
        max_t=max(max_t,abs(float(v)-float(t[i])))
        for j in range(21):
            ref=f'{chr(66+j)}{row}';v,typ,fmt=cells[ref]
            if typ not in ('n',None):raise ValueError(f'Non-numeric moisture {ref}')
            max_c=max(max_c,abs(float(v)-float(C[i,j])))
            if fmt!='0.0000':bad_format.append((ref,fmt))
    if max_t!=0 or max_c!=0:raise ValueError(f'Numeric mismatch max_t={max_t}, max_C={max_c}')
    if bad_format:raise ValueError(f'Wrong moisture format: {bad_format[:5]}')
    if C[-1].max()>=0.15:raise ValueError('Final unrounded field not strictly qualified')
    return {'passed':True,'sheet':'Sheet1','rows_including_header':n,'columns':22,'data_rows':len(t),'moisture_cells':int(C.size),'max_time_difference_s':max_t,'max_moisture_difference_kgkg':max_c,'last_time_s':float(t[-1]),'last_max_C':float(C[-1].max())}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--xlsx',type=Path,required=True);p.add_argument('--solution',type=Path,required=True);p.add_argument('--json-out',type=Path);a=p.parse_args();obj=validate(a.xlsx,a.solution);txt=json.dumps(obj,ensure_ascii=False,indent=2);print(txt);a.json_out and a.json_out.write_text(txt+'\n',encoding='utf-8')
