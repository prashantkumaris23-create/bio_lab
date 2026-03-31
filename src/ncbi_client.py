from __future__ import annotations

import xml.etree.ElementTree as ET
from typing import Iterable, List

import requests

from src.alignment import SequenceRecord, parse_fasta


ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def parse_accession_input(raw_value: str) -> List[str]:
    separators_normalized = raw_value.replace("\n", ",").replace(";", ",")
    items = [item.strip() for item in separators_normalized.split(",")]
    return [item for item in items if item]


def _request_params(email: str = "", api_key: str = "") -> dict[str, str]:
    params: dict[str, str] = {}
    if email.strip():
        params["email"] = email.strip()
    if api_key.strip():
        params["api_key"] = api_key.strip()
    return params


def search_nucleotide(
    term: str,
    retmax: int = 5,
    email: str = "",
    api_key: str = "",
) -> List[str]:
    if not term.strip():
        raise ValueError("NCBI search term is empty.")

    params = {
        "db": "nucleotide",
        "term": term.strip(),
        "retmode": "xml",
        "retmax": str(retmax),
    }
    params.update(_request_params(email, api_key))

    response = requests.get(ESEARCH_URL, params=params, timeout=60)
    response.raise_for_status()

    xml_root = ET.fromstring(response.text)
    ids = [element.text.strip() for element in xml_root.findall(".//IdList/Id") if element.text]
    if not ids:
        raise ValueError("No nucleotide records were found for the given NCBI query.")
    return ids


def fetch_fasta_by_ids(
    identifiers: Iterable[str],
    email: str = "",
    api_key: str = "",
) -> List[SequenceRecord]:
    ids = [identifier.strip() for identifier in identifiers if identifier.strip()]
    if not ids:
        raise ValueError("No accession IDs were provided.")

    params = {
        "db": "nucleotide",
        "id": ",".join(ids),
        "rettype": "fasta",
        "retmode": "text",
    }
    params.update(_request_params(email, api_key))

    response = requests.get(EFETCH_URL, params=params, timeout=120)
    response.raise_for_status()
    if not response.text.strip().startswith(">"):
        raise ValueError("NCBI did not return FASTA data. Please check the provided accessions or query.")

    records = parse_fasta(response.text)
    if len(records) < 2:
        raise ValueError("NCBI returned fewer than 2 FASTA records.")
    return records
