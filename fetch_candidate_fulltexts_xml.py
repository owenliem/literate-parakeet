#!/usr/bin/env python3
from pathlib import Path
import urllib.request

OUT=Path('candidate_fulltexts_xml'); OUT.mkdir(exist_ok=True)
for pmcid in ['PMC12458711','PMC12752622']:
    url=f'https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML'
    req=urllib.request.Request(url,headers={'User-Agent':'ctDNA-PRISMA-audit/1.0'})
    try:
        with urllib.request.urlopen(req,timeout=120) as r:
            data=r.read()
        (OUT/f'{pmcid}.xml').write_bytes(data)
        (OUT/f'{pmcid}.url.txt').write_text(url)
        print(pmcid,len(data))
    except Exception as e:
        (OUT/f'{pmcid}.error.txt').write_text(repr(e))
        print(pmcid,'ERROR',e)
