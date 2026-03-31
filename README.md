# Comparative Genomics Tree Builder

This project is a Streamlit-based bioinformatics app for comparative genomics and phylogenetic tree building.

It can:

- fetch DNA sequences from NCBI
- accept accession IDs or pasted FASTA data
- align sequences
- calculate similarity and distance matrices
- build a Neighbor Joining phylogenetic tree
- visualize results in a web interface

## Frontend Used

The frontend of this project is built with:

- Streamlit
  - used to create the web app interface
- Custom CSS
  - used to style the app, layout, colors, cards, buttons, and dark theme
- Streamlit charts and tables
  - used for monitoring graphs, data tables, and result views
- Custom SVG rendering
  - used to display the phylogenetic tree visually inside the app

## Backend Used

The backend logic of this project is built with:

- Python
  - used for all application logic
- Requests
  - used to call NCBI APIs
- NCBI Entrez E-utilities
  - used to search and fetch nucleotide sequence data
- FASTA parsing logic
  - used to read and validate DNA sequence records
- Needleman-Wunsch algorithm
  - used for global pairwise sequence alignment
- Star-style multiple alignment logic
  - used to build a simple multi-sequence alignment
- Similarity and distance matrix calculation
  - used for comparative analysis
- Neighbor Joining algorithm
  - used to generate the phylogenetic tree

## Main Files

- [app.py](C:\Users\PK\Documents\CODING\tree-builder\app.py)
  - main Streamlit application
- [src/alignment.py](C:\Users\PK\Documents\CODING\tree-builder\src\alignment.py)
  - FASTA parsing and sequence alignment
- [src/ncbi_client.py](C:\Users\PK\Documents\CODING\tree-builder\src\ncbi_client.py)
  - NCBI API fetching
- [src/phylo.py](C:\Users\PK\Documents\CODING\tree-builder\src\phylo.py)
  - similarity, distance matrix, and tree generation
- [src/visualization.py](C:\Users\PK\Documents\CODING\tree-builder\src\visualization.py)
  - heatmap and phylogenetic tree rendering
- [src/sample_projects.py](C:\Users\PK\Documents\CODING\tree-builder\src\sample_projects.py)
  - built-in example datasets

## Run The Project

Install dependencies:

```bash
pip install -r requirements.txt
```

Run the app:

```bash
python -m streamlit run app.py --server.address 127.0.0.1 --server.port 8501
```

Open in browser:

```text
http://127.0.0.1:8501
```
