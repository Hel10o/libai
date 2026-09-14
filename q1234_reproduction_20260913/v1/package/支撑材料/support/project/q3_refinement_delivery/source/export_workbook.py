"""Export the frozen Q3 refinement workbook from output/solution.npz.

The official workbook in this package was produced with artifact_tool.  This
module provides the same artifact_tool path plus a standard-library OOXML
fallback for ordinary local reproduction.  No formulas are used.
"""
from __future__ import annotations
import argparse, io, json, os, zipfile
from pathlib import Path
import numpy as np


def payload(root: Path):
    a=np.load(root/'output/solution.npz')
    t=a['time_s'][1:]
    c=a['profile_TC'][1:,:,1]
    if c.shape!=(len(t),21) or not np.all(np.diff(t)>0) or not np.isfinite(c).all():
        raise ValueError('Invalid final data')
    if abs(t[0]-60)>1e-10:
        raise ValueError('Workbook must start at 60 s')
    if not np.array_equal(t[:-1],np.arange(1,len(t))*60.0):
        raise ValueError('Missing 60 s rows')
    if c[-1].max()>=0.15:
        raise ValueError('Final row is not strictly qualified')
    return t,c


def export_artifact(root: Path, out: Path):
    if out.exists(): raise FileExistsError(out)
    os.environ.setdefault('ARTIFACT_TOOL_RPC_DAEMON_STARTUP_TIMEOUT_S','35')
    from artifact_tool import Blob, SpreadsheetFile
    root=Path(root);t,c=payload(root);n=len(t)+1
    wb=SpreadsheetFile.import_xlsx(Blob.load(str(root/'inputs/result3_blank.xlsx')))
    sh=wb.worksheets.get_item('Sheet1')
    sh.get_range(f'A1:V{n}').values=[['时间/s\\到药材中心的距离/cm']+np.linspace(0,2,21).tolist()]+np.column_stack([t,c]).tolist()
    sh.get_range(f'A2:A{n}').format.number_format='0.00'
    sh.get_range(f'B2:V{n}').format.number_format='0.0000'
    SpreadsheetFile.export_xlsx(wb).save(str(out))
    return {'engine':'artifact_tool','rows':len(t),'columns':22,'bytes':out.stat().st_size}


def portable_bytes(root: Path):
    t,c=payload(root);n=len(t)+1;ns='http://schemas.openxmlformats.org/spreadsheetml/2006/main'
    parts={
      '[Content_Types].xml':'''<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="xml" ContentType="application/xml"/><Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/><Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/><Override PartName="/xl/worksheets/sheet1.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/></Types>''',
      '_rels/.rels':'''<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/></Relationships>''',
      'xl/workbook.xml':f'''<workbook xmlns="{ns}" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets><sheet name="Sheet1" sheetId="1" r:id="rId1"/></sheets></workbook>''',
      'xl/_rels/workbook.xml.rels':'''<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet1.xml"/><Relationship Id="rId2" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/></Relationships>''',
      'xl/styles.xml':f'''<styleSheet xmlns="{ns}"><numFmts count="2"><numFmt numFmtId="164" formatCode="0.0000"/><numFmt numFmtId="165" formatCode="0.00"/></numFmts><fonts count="1"><font><sz val="10"/><name val="Arial"/></font></fonts><fills count="2"><fill><patternFill patternType="none"/></fill><fill><patternFill patternType="gray125"/></fill></fills><borders count="1"><border><left/><right/><top/><bottom/><diagonal/></border></borders><cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs><cellXfs count="3"><xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/><xf numFmtId="164" fontId="0" fillId="0" borderId="0" xfId="0" applyNumberFormat="1"/><xf numFmtId="165" fontId="0" fillId="0" borderId="0" xfId="0" applyNumberFormat="1"/></cellXfs><cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles></styleSheet>'''
    }
    b=io.BytesIO()
    with zipfile.ZipFile(b,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for name,text in parts.items():
            z.writestr(name,('<?xml version="1.0" encoding="UTF-8"?>'+text).encode())
        with z.open('xl/worksheets/sheet1.xml','w') as f:
            hdr='<worksheet xmlns="%s"><dimension ref="A1:V%d"/><sheetData><row r="1"><c r="A1" t="inlineStr"><is><t>时间/s\\到药材中心的距离/cm</t></is></c>'%(ns,n)
            hdr+=''.join(f'<c r="{chr(66+j)}1" t="n"><v>{j/10:.1f}</v></c>' for j in range(21))+'</row>'
            f.write(hdr.encode())
            for i,(tt,cc) in enumerate(zip(t,c),2):
                row=f'<row r="{i}"><c r="A{i}" t="n" s="2"><v>{float(tt):.17g}</v></c>'
                row+=''.join(f'<c r="{chr(66+j)}{i}" t="n" s="1"><v>{float(v):.17g}</v></c>' for j,v in enumerate(cc))+'</row>'
                f.write(row.encode())
            f.write(b'</sheetData></worksheet>')
    return b.getvalue()


def export_portable(root: Path, out: Path):
    if out.exists(): raise FileExistsError(out)
    out.write_bytes(portable_bytes(root))
    return {'engine':'standard-library OOXML','bytes':out.stat().st_size}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1]);p.add_argument('--out',type=Path,required=True);p.add_argument('--engine',choices=['artifact','portable'],default='portable');a=p.parse_args()
    print(json.dumps(export_artifact(a.root,a.out) if a.engine=='artifact' else export_portable(a.root,a.out)))
