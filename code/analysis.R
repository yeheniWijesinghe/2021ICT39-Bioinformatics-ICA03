# ======================================================================
# SETUP
# ======================================================================
# ---- libraries ----
# dplyr is loaded LAST so its functions are not masked by Bioconductor ones.
suppressPackageStartupMessages({
  library(rentrez)     # download sequences from NCBI GenBank
  library(Biostrings)  # DNA sequence objects, GC content
  library(DECIPHER)    # multiple sequence alignment
  library(ape)         # genetic distances, neighbour-joining tree, bootstrap
  library(phangorn)    # midpoint rooting of the tree
  library(ggplot2)
  library(tibble)
  library(tidyr)
  library(readr)
  library(stringr)
  library(dplyr)
})
set.seed(2026)  # makes random sampling and bootstrapping reproducible

# ---- paths-and-settings ----
# Project folders (works whether knitted from code/ or from the repo root)
proj_root <- if (basename(getwd()) == "code") normalizePath("..") else normalizePath(".")
data_dir <- file.path(proj_root, "data")
fig_dir  <- file.path(proj_root, "figures")
res_dir  <- file.path(proj_root, "results")
for (d in c(data_dir, fig_dir, res_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# Analysis settings - change here, not further down
max_records     <- 1000  # max GenBank records downloaded
min_len         <- 750   # the PRSV CP gene is about 850-870 bp
max_len         <- 1000
max_per_country <- 5     # cap so that heavily sampled countries do not dominate
min_per_region  <- 5     # a region needs at least this many sequences
n_boot          <- 200   # bootstrap replicates for the tree

# Countries grouped into world regions
regions <- list(
  "South Asia"     = c("India", "Sri Lanka", "Bangladesh", "Pakistan", "Nepal", "Bhutan", "Maldives"),
  "Southeast Asia" = c("Thailand", "Viet Nam", "Vietnam", "Philippines", "Indonesia", "Malaysia",
                       "Myanmar", "Cambodia", "Laos", "Lao People's Democratic Republic", "Singapore"),
  "East Asia"      = c("China", "Taiwan", "Japan", "South Korea", "Korea", "Hong Kong"),
  "Middle East"    = c("Iran", "Saudi Arabia", "Oman", "Yemen", "Turkey", "Israel", "Iraq",
                       "Jordan", "United Arab Emirates"),
  "Africa"         = c("Egypt", "Sudan", "Kenya", "Tanzania", "Uganda", "Nigeria", "Ghana",
                       "Cameroon", "South Africa", "Ethiopia", "Reunion", "Mauritius", "Madagascar",
                       "Benin", "Togo", "Cote d'Ivoire", "Rwanda", "Mozambique", "Zambia",
                       "Malawi", "Burkina Faso", "Senegal", "Tunisia", "Morocco"),
  "Americas"       = c("USA", "Mexico", "Brazil", "Venezuela", "Colombia", "Ecuador", "Peru",
                       "Cuba", "Puerto Rico", "Jamaica", "Guatemala", "Honduras", "Costa Rica",
                       "Panama", "Dominican Republic", "Trinidad and Tobago", "Guadeloupe",
                       "Martinique", "Argentina", "Bolivia", "Nicaragua", "El Salvador",
                       "Belize", "Haiti", "Barbados", "French Guiana", "Chile", "Canada"),
  "Oceania"        = c("Australia", "New Zealand", "Papua New Guinea", "Fiji",
                       "French Polynesia", "New Caledonia"),
  "Europe"         = c("Italy", "Spain", "France", "Germany", "Greece", "Netherlands",
                       "Portugal", "United Kingdom")
)
region_lookup <- setNames(rep(names(regions), lengths(regions)), unlist(regions))

# Helper: save a ggplot to figures/ (300 dpi PNG for the README) and show it
save_gg <- function(p, name, w = 8, h = 5.5) {
  ggsave(file.path(fig_dir, name), p, width = w, height = h, dpi = 300, bg = "white")
  print(p)        # show the plot in RStudio's Plots pane
  invisible(p)
}


# ======================================================================
# DATA IMPORT (NCBI GENBANK)
# ======================================================================
# Records are downloaded in **GenBank format** (not plain FASTA) because that
# format includes the **country** where each virus sample was collected. The
# download is saved in `data/`. If the file already exists it is reused, so the
# analysis always runs on the same data (GenBank grows every day).
# ---- download ----
raw_file <- file.path(data_dir, "PRSV_CP_raw.gb")
term <- '"Papaya ringspot virus"[Organism] AND "coat protein"[Title] AND 700:1100[SLEN]'

if (!file.exists(raw_file)) {
  res <- entrez_search(db = "nuccore", term = term, use_history = TRUE, retmax = 0)
  n   <- min(res$count, max_records)
  chunks <- character(0)
  for (start in seq(0, n - 1, by = 100)) {
    chunks <- c(chunks, entrez_fetch(db = "nuccore", web_history = res$web_history,
                                     rettype = "gb", retmode = "text",
                                     retstart = start, retmax = min(100, n - start)))
    Sys.sleep(0.5)  # be polite to the NCBI server
  }
  writeLines(chunks, raw_file)
  writeLines(c(paste("Query:", term), paste("Records found:", res$count),
               paste("Downloaded:", n), paste("Date:", Sys.time())),
             file.path(data_dir, "PRSV_CP_download_info.txt"))
}
gb_text <- paste(readLines(raw_file, warn = FALSE), collapse = "\n")

# ---- parse ----
# Split the file into one text block per GenBank record ("//" ends a record)
records <- strsplit(gb_text, "\n//")[[1]]
records <- records[grepl("LOCUS", records)]

# Pull out the fields we need from one record
parse_record <- function(r) {
  seq_txt <- str_match(r, "\nORIGIN[^\n]*\n([\\s\\S]*)$")[, 2]
  tibble(
    id          = str_match(r, "\nVERSION\\s+(\\S+)")[, 2],
    definition  = str_squish(str_match(r, "DEFINITION\\s+([\\s\\S]*?)\nACCESSION")[, 2]),
    # newer GenBank records use /geo_loc_name, older ones /country
    country_raw = str_squish(str_match(r, '/(?:geo_loc_name|country)="([^"]*)"')[, 2]),
    host        = str_squish(str_match(r, '/host="([^"]*)"')[, 2]),
    year        = str_extract(str_match(r, '/collection_date="([^"]*)"')[, 2], "\\d{4}"),
    seq         = toupper(gsub("[^A-Za-z]", "", ifelse(is.na(seq_txt), "", seq_txt)))
  )
}
info <- bind_rows(lapply(records, parse_record))
cat("GenBank records read:", nrow(info), "\n")