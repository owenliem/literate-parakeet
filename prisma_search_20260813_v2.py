#!/usr/bin/env python3
"""Fielded, reproducible database search for early ctDNA clearance during ICI in advanced NSCLC."""
from __future__ import annotations
import csv, datetime as dt, html, json, re, time, urllib.parse, urllib.request, xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any

OUT=Path('prisma_search_outputs_v2'); OUT.mkdir(exist_ok=True)
RUN=dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat()
UA='ctDNA-NSCLC-PRISMA/2.0 (systematic review audit; wynne.wijaya@oncology.ox.ac.uk)'

PUBMED='''("Carcinoma, Non-Small-Cell Lung"[Mesh] OR NSCLC[tiab] OR "non-small cell lung"[tiab] OR "non-small-cell lung"[tiab]) AND ("Circulating Tumor DNA"[Mesh] OR ctDNA[tiab] OR "circulating tumor DNA"[tiab] OR "cell-free tumor DNA"[tiab]) AND ("Immune Checkpoint Inhibitors"[Mesh] OR immunotherapy[tiab] OR pembrolizumab[tiab] OR nivolumab[tiab] OR atezolizumab[tiab] OR durvalumab[tiab] OR sintilimab[tiab] OR "PD-1"[tiab] OR "PD-L1"[tiab]) AND (clearance[tiab] OR undetectable[tiab] OR "molecular response"[tiab] OR kinetics[tiab] OR dynamics[tiab])'''

EPMC='''(TITLE_ABS:NSCLC OR TITLE_ABS:"non-small cell lung" OR TITLE_ABS:"non-small-cell lung") AND (TITLE_ABS:ctDNA OR TITLE_ABS:"circulating tumor DNA" OR TITLE_ABS:"cell-free tumor DNA") AND (TITLE_ABS:immunotherapy OR TITLE_ABS:"immune checkpoint" OR TITLE_ABS:pembrolizumab OR TITLE_ABS:nivolumab OR TITLE_ABS:atezolizumab OR TITLE_ABS:durvalumab OR TITLE_ABS:sintilimab OR TITLE_ABS:"PD-1" OR TITLE_ABS:"PD-L1") AND (TITLE_ABS:clearance OR TITLE_ABS:undetectable OR TITLE_ABS:"molecular response" OR TITLE_ABS:kinetics OR TITLE_ABS:dynamics)'''

CT_COND='non-small cell lung cancer'
CT_INTR='immunotherapy OR immune checkpoint OR PD-1 OR PD-L1 OR pembrolizumab OR nivolumab OR atezolizumab OR durvalumab OR sintilimab'
CT_TERM='("circulating tumor DNA" OR ctDNA) AND (clearance OR undetectable OR "molecular response" OR kinetics OR dynamics)'

def get(url:str)->bytes:
    err=None
    for i in range(5):
        try:
            req=urllib.request.Request(url,headers={'User-Agent':UA,'Accept':'application/json, application/xml, text/xml'})
            with urllib.request.urlopen(req,timeout=90) as r:return r.read()
        except Exception as e:
            err=e; time.sleep(min(2**(i+1),20))
    raise RuntimeError(url) from err

def js(url):return json.loads(get(url).decode())
def txt(n):return ''.join(n.itertext()).strip() if n is not None else ''
def nd(s):
    s=(s or '').strip().lower(); s=re.sub(r'^https?://(?:dx\.)?doi\.org/','',s); return re.sub(r'^doi:\s*','',s).rstrip('. ')
def nt(s):
    s=html.unescape(s or '').lower(); s=re.sub(r'<[^>]+>',' ',s); return re.sub(r'\s+',' ',re.sub(r'[^a-z0-9]+',' ',s)).strip()
def outcsv(name,rows):
    fields=['database','source_id','pmid','pmcid','doi','title','abstract','authors','journal','year','record_url','nct_id','overall_status','databases_found','dedup_key']
    with (OUT/name).open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore'); w.writeheader()
        for r in rows:
            q={k:(json.dumps(v,ensure_ascii=False) if isinstance(v,(list,dict)) else v) for k,v in r.items()}; w.writerow(q)

def pubmed():
    base='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
    u=base+'esearch.fcgi?'+urllib.parse.urlencode({'db':'pubmed','term':PUBMED,'retmode':'json','retmax':'100000','sort':'pub date'})
    p=js(u)['esearchresult']; ids=p['idlist']; rows=[]
    for a in range(0,len(ids),200):
        fu=base+'efetch.fcgi?'+urllib.parse.urlencode({'db':'pubmed','id':','.join(ids[a:a+200]),'retmode':'xml'})
        root=ET.fromstring(get(fu))
        for x in root.findall('.//PubmedArticle'):
            c=x.find('MedlineCitation'); art=c.find('Article') if c is not None else None
            pmid=txt(c.find('PMID')) if c is not None else ''; doi=''; pmcid=''
            for z in x.findall('PubmedData/ArticleIdList/ArticleId'):
                if z.attrib.get('IdType')=='doi':doi=nd(txt(z))
                elif z.attrib.get('IdType')=='pmc':pmcid=txt(z)
            authors=[]
            if art is not None:
                for au in art.findall('AuthorList/Author'):
                    n=txt(au.find('CollectiveName')) or ' '.join(filter(None,[txt(au.find('LastName')),txt(au.find('Initials'))]));
                    if n:authors.append(n)
            rows.append({'database':'PubMed/MEDLINE','source_id':pmid,'pmid':pmid,'pmcid':pmcid,'doi':doi,'title':txt(art.find('ArticleTitle')) if art is not None else '', 'abstract':' '.join(txt(q) for q in art.findall('Abstract/AbstractText')) if art is not None else '', 'authors':'; '.join(authors),'journal':txt(art.find('Journal/Title')) if art is not None else '', 'year':txt(art.find('Journal/JournalIssue/PubDate/Year')) if art is not None else '', 'record_url':f'https://pubmed.ncbi.nlm.nih.gov/{pmid}/'})
        time.sleep(.35)
    return rows,{'query':PUBMED,'url':u,'reported_count':int(p['count']),'retrieved_count':len(rows)}

def epmc():
    endpoint='https://www.ebi.ac.uk/europepmc/webservices/rest/search'; cursor='*'; rows=[]; count=None; first=''; pages=0
    while True:
        pages+=1; u=endpoint+'?'+urllib.parse.urlencode({'query':EPMC,'format':'json','pageSize':'1000','cursorMark':cursor,'resultType':'core'})
        first=first or u; p=js(u); count=int(p['hitCount']) if count is None else count; rr=p.get('resultList',{}).get('result',[])
        for x in rr:
            src=str(x.get('source') or ''); sid=str(x.get('id') or x.get('pmid') or x.get('doi') or '')
            rows.append({'database':'Europe PMC','source_id':sid,'pmid':str(x.get('pmid') or ''),'pmcid':str(x.get('pmcid') or ''),'doi':nd(str(x.get('doi') or '')),'title':str(x.get('title') or ''),'abstract':str(x.get('abstractText') or ''),'authors':str(x.get('authorString') or ''),'journal':str(x.get('journalTitle') or ''),'year':str(x.get('pubYear') or ''),'record_url':f'https://europepmc.org/article/{src}/{sid}'})
        nxt=p.get('nextCursorMark');
        if not rr or not nxt or nxt==cursor:break
        cursor=nxt; time.sleep(.2)
    return rows,{'query':EPMC,'url':first,'reported_count':count or 0,'retrieved_count':len(rows),'pages':pages}

def ctgov():
    endpoint='https://clinicaltrials.gov/api/v2/studies'; token=''; rows=[]; count=None; first=''; pages=0
    while True:
        pages+=1; prm={'format':'json','pageSize':'100','countTotal':'true','query.cond':CT_COND,'query.intr':CT_INTR,'query.term':CT_TERM}
        if token:prm['pageToken']=token
        u=endpoint+'?'+urllib.parse.urlencode(prm); first=first or u; p=js(u); count=int(p.get('totalCount',0)) if count is None else count
        for s in p.get('studies',[]):
            q=s.get('protocolSection',{}); ident=q.get('identificationModule',{}); desc=q.get('descriptionModule',{}); stat=q.get('statusModule',{}); nct=str(ident.get('nctId') or '')
            rows.append({'database':'ClinicalTrials.gov','source_id':nct,'nct_id':nct,'title':str(ident.get('officialTitle') or ident.get('briefTitle') or ''),'abstract':' '.join(filter(None,[str(desc.get('briefSummary') or ''),str(desc.get('detailedDescription') or '')])),'authors':str(ident.get('organization',{}).get('fullName') or ''),'journal':'ClinicalTrials.gov','year':str(stat.get('studyFirstPostDateStruct',{}).get('date') or '')[:4],'overall_status':str(stat.get('overallStatus') or ''),'record_url':f'https://clinicaltrials.gov/study/{nct}'})
        token=str(p.get('nextPageToken') or '')
        if not token:break
    return rows,{'condition_query':CT_COND,'intervention_query':CT_INTR,'term_query':CT_TERM,'url':first,'reported_count':count or 0,'retrieved_count':len(rows),'pages':pages}

def key(r):
    if r['database']=='ClinicalTrials.gov':return 'nct:'+r.get('nct_id','').lower()
    if r.get('pmid'):return 'pmid:'+r['pmid']
    if r.get('doi'):return 'doi:'+nd(r['doi'])
    return 'title:'+nt(r.get('title',''))
def dedup(rows):
    g={}
    for r in rows:g.setdefault(key(r),[]).append(r)
    uq=[]; dm=[]
    for k,v in g.items():
        keep=sorted(v,key=lambda r:(r['database']!='PubMed/MEDLINE',not bool(r.get('abstract')),-len(r.get('abstract',''))))[0].copy(); keep['dedup_key']=k; keep['databases_found']=sorted({x['database'] for x in v}); uq.append(keep)
        for d in v:
            if d is not keep and (d['database']!=keep['database'] or d['source_id']!=keep['source_id']):dm.append({'dedup_key':k,'kept_database':keep['database'],'kept_source_id':keep['source_id'],'duplicate_database':d['database'],'duplicate_source_id':d['source_id'],'title':keep['title']})
    return uq,dm

def main():
    pm,pl=pubmed(); ep,el=epmc(); ct,cl=ctgov(); raw=pm+ep+ct; uq,dm=dedup(raw)
    outcsv('pubmed_records.csv',pm); outcsv('europe_pmc_records.csv',ep); outcsv('clinicaltrials_records.csv',ct); outcsv('deduplicated_records.csv',uq)
    with (OUT/'duplicate_map.csv').open('w',encoding='utf-8',newline='') as f:
        fs=['dedup_key','kept_database','kept_source_id','duplicate_database','duplicate_source_id','title']; w=csv.DictWriter(f,fieldnames=fs);w.writeheader();w.writerows(dm)
    s={'run_at_utc':RUN,'sources':{'PubMed_MEDLINE':pl,'Europe_PMC':el,'ClinicalTrials_gov':cl},'raw_records_total':len(raw),'unique_records_after_cross_database_deduplication':len(uq),'duplicate_records_removed':len(raw)-len(uq),'publication_records_after_deduplication':sum(r['database']!='ClinicalTrials.gov' for r in uq),'trial_register_records_after_deduplication':sum(r['database']=='ClinicalTrials.gov' for r in uq),'deduplication_hierarchy':['NCT ID for registry records','PMID','DOI','normalized title'],'notes':['Counts are records, not unique studies.','Europe PMC query is restricted to title/abstract fields.','ClinicalTrials.gov records remain separate from publications.']}
    (OUT/'search_counts.json').write_text(json.dumps(s,indent=2,ensure_ascii=False)); (OUT/'search_strategies.json').write_text(json.dumps({'run_at_utc':RUN,'PubMed_MEDLINE':PUBMED,'Europe_PMC':EPMC,'ClinicalTrials_gov':{'condition':CT_COND,'intervention':CT_INTR,'term':CT_TERM}},indent=2,ensure_ascii=False)); print(json.dumps(s,indent=2))
if __name__=='__main__':main()
