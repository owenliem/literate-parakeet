#!/usr/bin/env python3
from __future__ import annotations
import pathlib, subprocess, urllib.request

OUT=pathlib.Path('candidate_fulltexts'); OUT.mkdir(exist_ok=True)
URLS={
 'ctMoniTR_Andrews_2025.pdf':'https://jitc.bmj.com/content/jitc/13/9/e012454.full.pdf',
 'Fei_2025.pdf':'https://advanced.onlinelibrary.wiley.com/doi/pdfdirect/10.1002/advs.202506565',
}
for name,url in URLS.items():
    p=OUT/name
    req=urllib.request.Request(url,headers={'User-Agent':'Mozilla/5.0 systematic-review-audit'})
    try:
        with urllib.request.urlopen(req,timeout=120) as r: p.write_bytes(r.read())
        print(name,p.stat().st_size)
    except Exception as e:
        (OUT/(name+'.error.txt')).write_text(repr(e))
        print(name,'ERROR',e)
for p in OUT.glob('*.pdf'):
    subprocess.run(['pdftotext','-layout',str(p),str(p.with_suffix('.txt'))],check=False)
