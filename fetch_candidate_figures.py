#!/usr/bin/env python3
from pathlib import Path
import urllib.request

OUT=Path('candidate_figures'); OUT.mkdir(exist_ok=True)
URLS={
 'Fei_2025_Figure5.jpg':'https://cdn.ncbi.nlm.nih.gov/pmc/blobs/5323/12752622/1f56db054923/ADVS-12-e06565-g001.jpg',
}
for name,url in URLS.items():
    try:
        req=urllib.request.Request(url,headers={'User-Agent':'Mozilla/5.0 ctDNA systematic review'})
        with urllib.request.urlopen(req,timeout=120) as r:data=r.read()
        (OUT/name).write_bytes(data)
        print(name,len(data))
    except Exception as e:
        (OUT/(name+'.error.txt')).write_text(repr(e))
        print(name,'ERROR',e)
