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


# DATA IMPORT (NCBI GENBANK)

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


# CLEANING AND QUALITY CONTROL

info <- info |>
  distinct(id, .keep_all = TRUE) |>
  mutate(length      = nchar(seq),
         acgt        = str_count(seq, "[ACGT]"),
         ambig_ok    = length > 0 & (1 - acgt / pmax(length, 1)) < 0.01,
         country     = ifelse(grepl("Hawaii", country_raw, ignore.case = TRUE), "USA (Hawaii)",
                              str_trim(sub(":.*$", "", country_raw))),
         region      = ifelse(country == "USA (Hawaii)", "Americas", unname(region_lookup[country])),
         has_cp      = grepl("coat protein", definition, ignore.case = TRUE),
         has_country = !is.na(country_raw),
         known_region = !is.na(region),
         length_ok   = length >= min_len & length <= max_len)

passed <- info |> filter(has_cp, has_country, known_region, length_ok, ambig_ok)

# Balance: at most max_per_country per country, then drop small regions
sampled <- passed |> group_by(country) |> slice_sample(n = max_per_country) |> ungroup()
big_regions <- count(sampled, region) |> filter(n >= min_per_region) |> pull(region)
kept <- sampled |>
  filter(region %in% big_regions) |>
  mutate(label = paste0(gsub("[^A-Za-z]", "", country), "_", id)) |>
  arrange(region, country, id)

steps <- tibble(
  step = c("GenBank records downloaded (unique)", "Title mentions coat protein",
           "+ country recorded", "+ country in a known region",
           paste0("+ length ", min_len, "-", max_len, " bp"), "+ <1% ambiguous bases",
           paste0("After sampling max ", max_per_country, " per country and keeping regions with >= ",
                  min_per_region, " (final)")),
  sequences = c(nrow(info), sum(info$has_cp), sum(info$has_cp & info$has_country),
                sum(info$has_cp & info$has_country & info$known_region),
                sum(info$has_cp & info$has_country & info$known_region & info$length_ok),
                nrow(passed), nrow(kept)))
knitr::kable(steps, caption = "Quality-control filtering")

# Stop with a clear explanation if too few regions are left to compare
if (n_distinct(kept$region) < 2) {
  stop("Fewer than 2 regions passed QC - nothing to compare.\n",
       paste(capture.output(print(as.data.frame(steps))), collapse = "\n"), "\n",
       "Sequences per region before dropping small regions:\n",
       paste(capture.output(print(as.data.frame(count(sampled, region)))), collapse = "\n"))
}

# Countries with no region in the lookup (add them above if important)
unknown <- info |> filter(has_cp, has_country, !known_region) |> count(country, sort = TRUE)
if (nrow(unknown) > 0) knitr::kable(unknown, caption = "Countries not assigned to a region (excluded)")

write_csv(dplyr::select(kept, label, id, country, region, host, year, length, definition),
          file.path(res_dir, "sequences_used.csv"))
seqs <- DNAStringSet(setNames(kept$seq, kept$label))
writeXStringSet(seqs, file.path(data_dir, "PRSV_CP_clean.fasta"))


# EXPLORATORY DATA ANALYSIS

knitr::kable(count(kept, region, country), caption = "Sequences analysed per region and country")

p_counts <- ggplot(count(kept, region), aes(reorder(region, n), n)) +
  geom_col(fill = "#4c72b0") + coord_flip() +
  labs(x = NULL, y = "Number of sequences", title = "PRSV CP sequences per region (after QC)") +
  theme_bw()
save_gg(p_counts, "fig1_sequences_per_region.png", h = 4.5)

# ---- eda-gc ----
kept <- kept |> mutate(gc = as.numeric(letterFrequency(seqs, "GC", as.prob = TRUE)))

p_gc <- ggplot(kept, aes(reorder(region, gc, FUN = median), 100 * gc)) +
  geom_boxplot(outlier.shape = NA, fill = "grey90") +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.7) + coord_flip() +
  labs(x = NULL, y = "GC content (%)", title = "GC content of the CP gene by region") +
  theme_bw()
save_gg(p_gc, "fig2_gc_by_region.png", h = 4.5)

# Does GC content differ between regions? (non-parametric: small groups)
kw_gc <- kruskal.test(gc ~ factor(region), data = kept)
kw_gc


# HAWAIIAN REFERENCE STRAIN

ref_file  <- file.path(data_dir, "PRSV_HA_X67673_genome.fasta")
ref_label <- "Hawaii_HA_X67673.1"
if (!file.exists(ref_file)) {
  writeLines(entrez_fetch(db = "nuccore", id = "X67673", rettype = "fasta"), ref_file)
}
genome <- readDNAStringSet(ref_file)

# Align the genome with the longest CP sequence; keep the genome part that
# lines up with it (from the first to the last aligned CP base)
query <- seqs[which.max(width(seqs))]
pair  <- AlignSeqs(OrientNucleotides(c(query, genome), verbose = FALSE),
                   verbose = FALSE, processors = 1)
pm    <- as.matrix(pair)
cols  <- which(!pm[1, ] %in% c("-", "."))
ha_cp <- gsub("[-.]", "", paste(pm[2, min(cols):max(cols)], collapse = ""))

has_ref <- nchar(ha_cp) >= 700 && nchar(ha_cp) <= 1100
cat("Hawaiian HA CP region extracted:", nchar(ha_cp), "bp\n")
if (!has_ref) warning("Could not extract the CP region from X67673 - Hawaii comparison skipped.")

seqs_aln <- if (has_ref) c(seqs, DNAStringSet(setNames(ha_cp, ref_label))) else seqs

# Lookup tables: label -> region / country (the reference gets its own group)
region_of  <- c(setNames(kept$region, kept$label), setNames("Hawaii (reference)", ref_label))
country_of <- c(setNames(kept$country, kept$label), setNames("USA (Hawaii)", ref_label))


# MULTIPLE SEQUENCE ALIGNMENT


seqs_o <- OrientNucleotides(seqs_aln, verbose = FALSE)
aln    <- AlignSeqs(seqs_o, verbose = FALSE, processors = 1)

aln_mat   <- as.matrix(aln)
keep_cols <- colMeans(aln_mat == "-" | aln_mat == ".") <= 0.5
aln_trim  <- DNAStringSet(apply(aln_mat[, keep_cols], 1, paste0, collapse = ""))
names(aln_trim) <- names(aln)
writeXStringSet(aln_trim, file.path(res_dir, "PRSV_CP_alignment.fasta"))
cat("Alignment length:", ncol(aln_mat), "columns; after trimming:", sum(keep_cols), "\n")

# Read the alignment into ape's format for distances and trees
dna <- read.dna(file.path(res_dir, "PRSV_CP_alignment.fasta"), format = "fasta")
reg <- unname(region_of[rownames(dna)])
stopifnot("Labels in alignment do not match the sequence table" = !anyNA(reg))


# GENETIC DISTANCES: ARE STRAINS GROUPED BY REGION?

Dm <- as.matrix(dist.dna(dna, model = "K80", pairwise.deletion = TRUE))
if (anyNA(Dm) || any(!is.finite(Dm))) stop("Some distances could not be computed.")
Dm[Dm < 0] <- 0  # identical sequences can give -0 from rounding

# H1 uses the regional sequences only (the Hawaiian reference is left out)
use   <- reg != "Hawaii (reference)"
Dm_r  <- Dm[use, use]
reg_r <- reg[use]

# Every pair of sequences: same region or different regions?
idx   <- which(upper.tri(Dm_r), arr.ind = TRUE)
pairs <- tibble(distance = Dm_r[idx],
                type = ifelse(reg_r[idx[, 1]] == reg_r[idx[, 2]], "Same region", "Different regions"))

pairs |> group_by(type) |>
  summarise(n_pairs = n(), mean = mean(distance), median = median(distance),
            max = max(distance)) |>
  knitr::kable(digits = 4, caption = "K2P distances within and between regions")

# H1: are sequences from different regions more distant? (one-sided Wilcoxon test)
wilcox_reg <- wilcox.test(distance ~ type, data = pairs, alternative = "greater")
wilcox_reg

# ---- distance-histogram ----
p_hist <- ggplot(pairs, aes(100 * distance, fill = type)) +
  geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
  scale_fill_manual(values = c("Same region" = "#1b9e77", "Different regions" = "#d95f02")) +
  labs(x = "K2P distance (%)", y = "Number of sequence pairs", fill = NULL,
       title = "CP gene distances within vs between regions") +
  theme_bw() + theme(legend.position = "top")
save_gg(p_hist, "fig3_distance_histogram.png", h = 4.5)

# ---- region-heatmap ----
# Mean distance between every pair of regions (diagonal = within a region)
reg_levels <- sort(unique(reg_r))
reg_df <- expand.grid(region1 = reg_levels, region2 = reg_levels, stringsAsFactors = FALSE)
reg_df$mean_distance <- mapply(function(a, b) {
  m <- Dm_r[reg_r == a, reg_r == b, drop = FALSE]
  if (a == b) mean(m[upper.tri(m)]) else mean(m)
}, reg_df$region1, reg_df$region2)
write_csv(reg_df, file.path(res_dir, "mean_distance_between_regions.csv"))

p_heat <- ggplot(reg_df, aes(region1, region2, fill = 100 * mean_distance)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.1f", 100 * mean_distance),
                colour = mean_distance > median(mean_distance)), size = 3.2,
            show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black")) +
  scale_fill_gradient(low = "#deebf7", high = "#08306b") +
  labs(x = NULL, y = NULL, fill = "Mean K2P\ndistance (%)",
       title = "Mean CP distance within and between regions") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_gg(p_heat, "fig4_region_distance_heatmap.png", h = 5.5)


# PHYLOGENETIC TREE (NEIGHBOUR-JOINING WITH BOOTSTRAP)

# ---- nj-tree ----
# NJ tree on K2P distances; njs() copes with any missing distances
nj_fun <- function(x) {
  d <- dist.dna(x, model = "K80", pairwise.deletion = TRUE)
  if (anyNA(d)) njs(d) else nj(d)
}
tree <- nj_fun(dna)

# Bootstrap: resample alignment columns n_boot times and rebuild the tree
bs <- boot.phylo(tree, dna, nj_fun, B = n_boot, quiet = TRUE)
tree$node.label <- round(100 * bs / n_boot)
tree_rooted <- ladderize(midpoint(tree, node.labels = "support"))
write.tree(tree_rooted, file.path(res_dir, "tree_NJ_bootstrap.nwk"))

tree_levels <- sort(unique(reg))  # includes the Hawaiian reference, if present
region_cols <- setNames(hcl.colors(length(tree_levels), "Dark 3"), tree_levels)

draw_tree <- function() {
  tips <- tree_rooted$tip.label
  plot(tree_rooted, cex = 0.5, no.margin = TRUE, label.offset = 0.001,
       tip.color = region_cols[region_of[tips]],
       font = ifelse(country_of[tips] %in% c("Sri Lanka", "USA (Hawaii)"), 2, 1))  # bold = Sri Lanka / Hawaii
  s <- suppressWarnings(as.numeric(tree_rooted$node.label))
  nodelabels(ifelse(!is.na(s) & s >= 70, s, ""), frame = "none", cex = 0.45, adj = c(1.2, -0.4))
  add.scale.bar(cex = 0.6)
  legend("bottomright", legend = tree_levels, text.col = region_cols, bty = "n", cex = 0.7)
}
png(file.path(fig_dir, "fig5_NJ_tree.png"), width = 8, height = max(8, 0.12 * Ntip(tree_rooted)),
    units = "in", res = 300)
draw_tree(); dev.off()
draw_tree()


# SIMILARITY TO THE HAWAIIAN STRAIN USED IN RESISTANT PAPAYA

P <- as.matrix(dist.dna(dna, model = "raw", pairwise.deletion = TRUE))
is_hi <- unname(country_of[rownames(P)]) == "USA (Hawaii)"
has_hawaii <- sum(is_hi) > 0
cat("Hawaiian reference sequences:", sum(is_hi), "\n")
if (!has_hawaii) cat("No Hawaiian sequences passed QC - the next chunk is skipped.\n")

# ---- hawaii ----
if (has_hawaii) {
  ident <- tibble(label    = rownames(P)[!is_hi],
                  region   = reg[!is_hi],
                  country  = unname(country_of[rownames(P)[!is_hi]]),
                  identity = 100 * (1 - apply(P[!is_hi, is_hi, drop = FALSE], 1, min)))
  write_csv(ident, file.path(res_dir, "identity_to_hawaii.csv"))

  ident |> group_by(region) |>
    summarise(n = n(), mean_identity = mean(identity), min_identity = min(identity)) |>
    arrange(mean_identity) |>
    knitr::kable(digits = 1, caption = "CP identity (%) to the closest Hawaiian strain")

  if (any(ident$country == "Sri Lanka")) {
    knitr::kable(filter(ident, country == "Sri Lanka"), digits = 1,
                 caption = "Sri Lankan sequences: identity to the Hawaiian strain")
  }

  # H2: does identity to the Hawaiian strain differ between regions?
  kw_hi <- kruskal.test(identity ~ factor(region), data = ident)
  kw_hi

  p_hi <- ggplot(ident, aes(reorder(region, identity, FUN = median), identity)) +
    geom_boxplot(outlier.shape = NA, fill = "grey90") +
    geom_jitter(width = 0.15, size = 1.2, alpha = 0.7) + coord_flip() +
    labs(x = NULL, y = "CP identity to closest Hawaiian strain (%)",
         title = "How similar are PRSV strains to the Hawaiian strain?") +
    theme_bw()
  save_gg(p_hi, "fig6_identity_to_hawaii.png", h = 4.5)
}