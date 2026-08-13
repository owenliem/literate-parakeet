#!/usr/bin/env python3
"""Reproducible database search for the ctDNA-ICI advanced NSCLC review.

Sources:
- PubMed/MEDLINE via NCBI E-utilities
- Europe PMC via the Europe PMC REST API
- ClinicalTrials.gov via API v2

The script writes raw source exports, a cross-database deduplicated export,
duplicate mappings, exact search strings, and a machine-readable count log.
"""

from __future__ import annotations

import csv
import datetime as dt
import html
import json
import re
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any, Iterable

OUT = Path("prisma_search_outputs")
OUT.mkdir(exist_ok=True)
RUN_AT_UTC = dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat()
USER_AGENT = "ctDNA-NSCLC-PRISMA-search/1.0 (systematic-review audit; contact: wynne.wijaya@oncology.ox.ac.uk)"

PUBMED_QUERY = (
    '("Carcinoma, Non-Small-Cell Lung"[Mesh] OR NSCLC[tiab] OR '
    '"non-small cell lung"[tiab] OR "non-small-cell lung"[tiab]) AND '
    '("Circulating Tumor DNA"[Mesh] OR ctDNA[tiab] OR '
    '"circulating tumor DNA"[tiab] OR "cell-free tumor DNA"[tiab]) AND '
    '("Immune Checkpoint Inhibitors"[Mesh] OR immunotherapy[tiab] OR '
    'pembrolizumab[tiab] OR nivolumab[tiab] OR atezolizumab[tiab] OR '
    'durvalumab[tiab] OR sintilimab[tiab] OR "PD-1"[tiab] OR "PD-L1"[tiab]) AND '
    '(clearance[tiab] OR undetectable[tiab] OR "molecular response"[tiab] OR '
    'kinetics[tiab] OR dynamics[tiab])'
)

EUROPE_PMC_QUERY = (
    '(NSCLC OR "non-small cell lung" OR "non-small-cell lung") AND '
    '(ctDNA OR "circulating tumor DNA" OR "cell-free tumor DNA") AND '
    '(immunotherapy OR "immune checkpoint" OR pembrolizumab OR nivolumab OR '
    'atezolizumab OR durvalumab OR sintilimab OR "PD-1" OR "PD-L1") AND '
    '(clearance OR undetectable OR "molecular response" OR kinetics OR dynamics)'
)

CTGOV_CONDITION = "non-small cell lung cancer"
CTGOV_TERM = (
    '("circulating tumor DNA" OR ctDNA) AND '
    '(immunotherapy OR "immune checkpoint" OR "PD-1" OR "PD-L1" OR '
    'pembrolizumab OR nivolumab OR atezolizumab OR durvalumab OR sintilimab) AND '
    '(clearance OR undetectable OR response OR "molecular response" OR kinetics OR dynamics)'
)


def request_bytes(url: str, *, attempts: int = 5) -> bytes:
    last: Exception | None = None
    for attempt in range(1, attempts + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT, "Accept": "application/json, application/xml, text/xml"})
            with urllib.request.urlopen(req, timeout=90) as response:
                return response.read()
        except Exception as exc:  # noqa: BLE001
            last = exc
            if attempt == attempts:
                break
            time.sleep(min(2**attempt, 20))
    raise RuntimeError(f"Request failed after {attempts} attempts: {url}") from last


def request_json(url: str) -> dict[str, Any]:
    return json.loads(request_bytes(url).decode("utf-8"))


def text_of(node: ET.Element | None) -> str:
    if node is None:
        return ""
    return "".join(node.itertext()).strip()


def normalize_doi(value: str) -> str:
    value = value.strip().lower()
    value = re.sub(r"^https?://(?:dx\.)?doi\.org/", "", value)
    value = re.sub(r"^doi:\s*", "", value)
    return value.rstrip(". ")


def normalize_title(value: str) -> str:
    value = html.unescape(value or "").lower()
    value = re.sub(r"<[^>]+>", " ", value)
    value = value.replace("–", "-").replace("—", "-")
    value = re.sub(r"[^a-z0-9]+", " ", value)
    return re.sub(r"\s+", " ", value).strip()


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            clean = {k: (json.dumps(v, ensure_ascii=False) if isinstance(v, (list, dict)) else v) for k, v in row.items()}
            writer.writerow(clean)


def pubmed_search() -> tuple[list[dict[str, Any]], dict[str, Any]]:
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    params = {
        "db": "pubmed",
        "term": PUBMED_QUERY,
        "retmode": "json",
        "retmax": "100000",
        "usehistory": "y",
        "sort": "pub+date",
    }
    search_url = base + "esearch.fcgi?" + urllib.parse.urlencode(params)
    payload = request_json(search_url)["esearchresult"]
    ids = payload.get("idlist", [])
    count = int(payload["count"])
    if len(ids) != count:
        raise RuntimeError(f"PubMed returned {len(ids)} IDs for count={count}")

    records: list[dict[str, Any]] = []
    for start in range(0, len(ids), 200):
        batch = ids[start : start + 200]
        fetch_url = base + "efetch.fcgi?" + urllib.parse.urlencode({"db": "pubmed", "id": ",".join(batch), "retmode": "xml"})
        root = ET.fromstring(request_bytes(fetch_url))
        for article in root.findall(".//PubmedArticle"):
            citation = article.find("MedlineCitation")
            art = citation.find("Article") if citation is not None else None
            pmid = text_of(citation.find("PMID")) if citation is not None else ""
            title = text_of(art.find("ArticleTitle")) if art is not None else ""
            abstract = " ".join(text_of(x) for x in art.findall("Abstract/AbstractText")) if art is not None else ""
            journal = text_of(art.find("Journal/Title")) if art is not None else ""
            year = ""
            if art is not None:
                year = text_of(art.find("Journal/JournalIssue/PubDate/Year")) or text_of(art.find("Journal/JournalIssue/PubDate/MedlineDate"))
            doi = ""
            pmcid = ""
            for aid in article.findall("PubmedData/ArticleIdList/ArticleId"):
                kind = aid.attrib.get("IdType", "")
                if kind == "doi":
                    doi = normalize_doi(text_of(aid))
                elif kind == "pmc":
                    pmcid = text_of(aid)
            pub_types = [text_of(x) for x in art.findall("PublicationTypeList/PublicationType")] if art is not None else []
            authors = []
            if art is not None:
                for author in art.findall("AuthorList/Author"):
                    collective = text_of(author.find("CollectiveName"))
                    if collective:
                        authors.append(collective)
                    else:
                        name = " ".join(filter(None, [text_of(author.find("LastName")), text_of(author.find("Initials"))]))
                        if name:
                            authors.append(name)
            records.append({
                "database": "PubMed/MEDLINE",
                "source_id": pmid,
                "pmid": pmid,
                "pmcid": pmcid,
                "doi": doi,
                "title": title,
                "abstract": abstract,
                "authors": "; ".join(authors),
                "journal": journal,
                "year": year,
                "publication_types": pub_types,
                "record_url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else "",
            })
        time.sleep(0.4)
    return records, {"query": PUBMED_QUERY, "url": search_url, "reported_count": count, "retrieved_count": len(records)}


def europe_pmc_search() -> tuple[list[dict[str, Any]], dict[str, Any]]:
    endpoint = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
    cursor = "*"
    records: list[dict[str, Any]] = []
    reported_count: int | None = None
    page = 0
    first_url = ""
    while cursor:
        page += 1
        params = {
            "query": EUROPE_PMC_QUERY,
            "format": "json",
            "pageSize": "1000",
            "cursorMark": cursor,
            "resultType": "core",
        }
        url = endpoint + "?" + urllib.parse.urlencode(params)
        if not first_url:
            first_url = url
        payload = request_json(url)
        if reported_count is None:
            reported_count = int(payload["hitCount"])
        results = payload.get("resultList", {}).get("result", [])
        for item in results:
            pmid = str(item.get("pmid") or "")
            doi = normalize_doi(str(item.get("doi") or ""))
            source = str(item.get("source") or "")
            source_id = str(item.get("id") or pmid or doi)
            records.append({
                "database": "Europe PMC",
                "source_id": source_id,
                "source": source,
                "pmid": pmid,
                "pmcid": str(item.get("pmcid") or ""),
                "doi": doi,
                "title": str(item.get("title") or ""),
                "abstract": str(item.get("abstractText") or ""),
                "authors": str(item.get("authorString") or ""),
                "journal": str(item.get("journalTitle") or ""),
                "year": str(item.get("pubYear") or ""),
                "publication_types": item.get("pubTypeList", {}).get("pubType", []),
                "record_url": f"https://europepmc.org/article/{source}/{source_id}" if source and source_id else "",
            })
        next_cursor = payload.get("nextCursorMark")
        if not results or not next_cursor or next_cursor == cursor:
            break
        cursor = next_cursor
        time.sleep(0.25)
    if reported_count is None:
        reported_count = 0
    return records, {"query": EUROPE_PMC_QUERY, "url": first_url, "reported_count": reported_count, "retrieved_count": len(records), "pages": page}


def ctgov_search() -> tuple[list[dict[str, Any]], dict[str, Any]]:
    endpoint = "https://clinicaltrials.gov/api/v2/studies"
    page_token = ""
    records: list[dict[str, Any]] = []
    reported_count: int | None = None
    first_url = ""
    page = 0
    while True:
        page += 1
        params = {
            "format": "json",
            "pageSize": "100",
            "countTotal": "true",
            "query.cond": CTGOV_CONDITION,
            "query.term": CTGOV_TERM,
        }
        if page_token:
            params["pageToken"] = page_token
        url = endpoint + "?" + urllib.parse.urlencode(params)
        if not first_url:
            first_url = url
        payload = request_json(url)
        if reported_count is None:
            reported_count = int(payload.get("totalCount", 0))
        studies = payload.get("studies", [])
        for study in studies:
            p = study.get("protocolSection", {})
            ident = p.get("identificationModule", {})
            desc = p.get("descriptionModule", {})
            cond = p.get("conditionsModule", {})
            arms = p.get("armsInterventionsModule", {})
            design = p.get("designModule", {})
            status = p.get("statusModule", {})
            elig = p.get("eligibilityModule", {})
            nct = str(ident.get("nctId") or "")
            interventions = []
            for inter in arms.get("interventions", []) or []:
                interventions.append({"type": inter.get("type"), "name": inter.get("name")})
            title = str(ident.get("officialTitle") or ident.get("briefTitle") or "")
            summary_parts = [str(desc.get("briefSummary") or ""), str(desc.get("detailedDescription") or ""), str(elig.get("eligibilityCriteria") or "")]
            records.append({
                "database": "ClinicalTrials.gov",
                "source_id": nct,
                "nct_id": nct,
                "pmid": "",
                "pmcid": "",
                "doi": "",
                "title": title,
                "abstract": " ".join(x for x in summary_parts if x),
                "authors": str(ident.get("organization", {}).get("fullName") or ""),
                "journal": "ClinicalTrials.gov",
                "year": str(status.get("studyFirstPostDateStruct", {}).get("date") or "")[:4],
                "publication_types": [str(design.get("studyType") or "Trial registration")],
                "conditions": cond.get("conditions", []) or [],
                "interventions": interventions,
                "overall_status": str(status.get("overallStatus") or ""),
                "record_url": f"https://clinicaltrials.gov/study/{nct}" if nct else "",
            })
        page_token = str(payload.get("nextPageToken") or "")
        if not page_token:
            break
        time.sleep(0.25)
    if reported_count is None:
        reported_count = 0
    return records, {
        "condition_query": CTGOV_CONDITION,
        "term_query": CTGOV_TERM,
        "url": first_url,
        "reported_count": reported_count,
        "retrieved_count": len(records),
        "pages": page,
    }


def record_key(row: dict[str, Any]) -> str:
    database = row.get("database", "")
    if database == "ClinicalTrials.gov":
        return "nct:" + str(row.get("nct_id") or row.get("source_id") or "").lower()
    pmid = str(row.get("pmid") or "").strip()
    if pmid:
        return "pmid:" + pmid
    doi = normalize_doi(str(row.get("doi") or ""))
    if doi:
        return "doi:" + doi
    title = normalize_title(str(row.get("title") or ""))
    return "title:" + title


def deduplicate(rows: Iterable[dict[str, Any]]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    groups: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        groups.setdefault(record_key(row), []).append(row)
    unique: list[dict[str, Any]] = []
    duplicate_map: list[dict[str, Any]] = []
    for key, group in groups.items():
        preferred = sorted(group, key=lambda r: (r.get("database") != "PubMed/MEDLINE", not bool(r.get("abstract")), -len(str(r.get("abstract") or ""))))[0].copy()
        preferred["dedup_key"] = key
        preferred["databases_found"] = sorted({str(g.get("database")) for g in group})
        preferred["source_ids"] = sorted({str(g.get("source_id") or "") for g in group if g.get("source_id")})
        unique.append(preferred)
        if len(group) > 1:
            for duplicate in group[1:]:
                duplicate_map.append({
                    "dedup_key": key,
                    "kept_database": preferred.get("database", ""),
                    "kept_source_id": preferred.get("source_id", ""),
                    "duplicate_database": duplicate.get("database", ""),
                    "duplicate_source_id": duplicate.get("source_id", ""),
                    "title": preferred.get("title", ""),
                })
    unique.sort(key=lambda r: (str(r.get("database")), str(r.get("year")), normalize_title(str(r.get("title") or ""))))
    return unique, duplicate_map


def main() -> None:
    pubmed, pubmed_log = pubmed_search()
    europe, europe_log = europe_pmc_search()
    ctgov, ctgov_log = ctgov_search()

    all_raw = pubmed + europe + ctgov
    unique, duplicate_map = deduplicate(all_raw)

    common_fields = [
        "database", "source_id", "source", "nct_id", "pmid", "pmcid", "doi", "title", "abstract",
        "authors", "journal", "year", "publication_types", "conditions", "interventions", "overall_status", "record_url",
    ]
    write_csv(OUT / "pubmed_records.csv", pubmed, common_fields)
    write_csv(OUT / "europe_pmc_records.csv", europe, common_fields)
    write_csv(OUT / "clinicaltrials_records.csv", ctgov, common_fields)
    write_csv(OUT / "deduplicated_records.csv", unique, common_fields + ["dedup_key", "databases_found", "source_ids"])
    write_csv(OUT / "duplicate_map.csv", duplicate_map, ["dedup_key", "kept_database", "kept_source_id", "duplicate_database", "duplicate_source_id", "title"])

    total_raw = len(all_raw)
    unique_count = len(unique)
    summary = {
        "run_at_utc": RUN_AT_UTC,
        "sources": {
            "PubMed_MEDLINE": pubmed_log,
            "Europe_PMC": europe_log,
            "ClinicalTrials_gov": ctgov_log,
        },
        "raw_records_total": total_raw,
        "unique_records_after_cross_database_deduplication": unique_count,
        "duplicate_records_removed": total_raw - unique_count,
        "publication_records_after_deduplication": sum(1 for r in unique if r.get("database") != "ClinicalTrials.gov"),
        "trial_register_records_after_deduplication": sum(1 for r in unique if r.get("database") == "ClinicalTrials.gov"),
        "deduplication_hierarchy": ["NCT identifier for trial-register records", "PMID", "DOI", "normalized title"],
        "notes": [
            "Counts are database records returned at the run timestamp, not numbers of unique studies.",
            "Europe PMC includes MED/PubMed records; cross-database duplicates are removed using PMID, then DOI, then normalized title.",
            "ClinicalTrials.gov registrations are retained as separate register records rather than merged into journal publications.",
        ],
    }
    (OUT / "search_counts.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    (OUT / "search_strategies.json").write_text(json.dumps({
        "run_at_utc": RUN_AT_UTC,
        "PubMed_MEDLINE": PUBMED_QUERY,
        "Europe_PMC": EUROPE_PMC_QUERY,
        "ClinicalTrials_gov": {"condition": CTGOV_CONDITION, "term": CTGOV_TERM},
    }, indent=2, ensure_ascii=False), encoding="utf-8")
    (OUT / "README.txt").write_text(
        "PRISMA search export generated at " + RUN_AT_UTC + "\n"
        "Files contain raw database exports, a deduplicated master export, duplicate mappings, exact strategies, and counts.\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
