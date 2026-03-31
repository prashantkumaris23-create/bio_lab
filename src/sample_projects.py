from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class SampleProject:
    key: str
    title: str
    subtitle: str
    description: str
    input_mode: str
    search_query: str = ""
    accessions_text: str = ""
    fasta_text: str = ""
    search_limit: int = 5


SAMPLE_PROJECTS = [
    SampleProject(
        key="primate_cox1",
        title="Primate COX1",
        subtitle="NCBI accession project",
        description="Compare 3 primate mitochondrial COX1-related nucleotide records fetched directly from NCBI.",
        input_mode="Accession IDs",
        accessions_text="PV276824.1\nPQ495639.1\nPQ495638.1",
    ),
    SampleProject(
        key="primates_search",
        title="Primates Search",
        subtitle="Live NCBI query project",
        description="Run a filtered NCBI nucleotide search and build a tree from returned mitochondrial gene matches.",
        input_mode="NCBI Search",
        search_query="COX1[Gene] AND mitochondrion[filter] AND primates[Organism] AND 500:3000[SLEN]",
        search_limit=5,
    ),
    SampleProject(
        key="toy_fasta",
        title="Toy FASTA Demo",
        subtitle="Instant local project",
        description="Use a built-in DNA FASTA dataset for a fast, offline-friendly example project.",
        input_mode="Paste FASTA",
        fasta_text=""">Species_A
ATGCTAGCTAGCTAAGCTTACGAT
>Species_B
ATGCTAGATAGCTAAGCTTATGAT
>Species_C
ATGCGAGATAGCTAAGCATATGAT
>Species_D
ATGCGAGATCGCTAAGCATATGAT
""",
    ),
]
