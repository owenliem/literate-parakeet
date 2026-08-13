#!/usr/bin/env python3
from pathlib import Path
import urllib.request, xml.etree.ElementTree as ET, tarfile, io

OUT=Path('candidate_pmc_oa'); OUT.mkdir(exist_ok=True)
UA='ctDNA-PRISMA-audit/1.0'
for pmcid in ['PMC12458711','PMC12752622']:
    try:
        meta_url=f'https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?id={pmcid}'
        with urllib.request.urlopen(urllib.request.Request(meta_url,headers={'User-Agent':UA}),timeout=120) as r:
            meta=r.read()
        (OUT/f'{pmcid}_oa.xml').write_bytes(meta)
        root=ET.fromstring(meta)
        links=root.findall('.//link')
        href=''
        for link in links:
            if link.attrib.get('format')=='tgz': href=link.attrib.get('href',''); break
        if not href: raise RuntimeError('No tgz OA link')
        if href.startswith('ftp://'): href='https://'+href[len('ftp://'):]
        with urllib.request.urlopen(urllib.request.Request(href,headers={'User-Agent':UA}),timeout=180) as r:
            data=r.read()
        tgz=OUT/f'{pmcid}.tar.gz'; tgz.write_bytes(data)
        dest=OUT/pmcid; dest.mkdir(exist_ok=True)
        with tarfile.open(fileobj=io.BytesIO(data),mode='r:gz') as tf:
            tf.extractall(dest)
        print(pmcid, len(data), href)
    except Exception as e:
        (OUT/f'{pmcid}.error.txt').write_text(repr(e))
        print(pmcid,'ERROR',e)
