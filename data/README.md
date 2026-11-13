# Data

Created by `code/analysis.Rmd` on its first run, then reused on every later run.

| File | Description |
|---|---|
| `PRSV_CP_raw.gb` | Raw *Papaya ringspot virus* coat-protein records from NCBI GenBank, in GenBank format (which includes the country of collection). Downloaded with `rentrez`. Query: `"Papaya ringspot virus"[Organism] AND "coat protein"[Title] AND 700:1100[SLEN]`, max. 1,000 records |
| `PRSV_CP_download_info.txt` | Query, number of records found and downloaded, and download date |
| `PRSV_CP_clean.fasta` | Sequences that passed quality control (the analysis dataset) |

The raw download is committed on purpose. GenBank changes every day, so keeping the exact file makes the results reproducible. To analyse the current GenBank data, delete `PRSV_CP_raw.gb` and re-run the notebook.

**Source and licence:** NCBI GenBank (https://www.ncbi.nlm.nih.gov/genbank/). NCBI places no restrictions on the use of GenBank data. Cite GenBank and the accession numbers listed in `results/sequences_used.csv`.
