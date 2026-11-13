# BiocManager installs Bioconductor packages (Biostrings, DECIPHER)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# CRAN packages: NCBI access, phylogenetics, data wrangling, plotting, reports
cran_pkgs <- c("rentrez", "ape", "phangorn",
               "dplyr", "tidyr", "tibble", "readr", "stringr", "ggplot2",
               "rmarkdown", "knitr")

# Bioconductor packages: DNA sequence handling and multiple alignment
bioc_pkgs <- c("Biostrings", "DECIPHER")

# Install only what is missing, so re-running the script is quick
missing_cran <- cran_pkgs[!vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran) > 0) install.packages(missing_cran)

missing_bioc <- bioc_pkgs[!vapply(bioc_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc) > 0) BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)

# Final check: report anything that still failed to install
all_pkgs <- c(cran_pkgs, bioc_pkgs)
ok <- vapply(all_pkgs, requireNamespace, logical(1), quietly = TRUE)
if (all(ok)) {
  message("All ", length(all_pkgs), " packages are installed.")
} else {
  warning("These packages did not install: ", paste(all_pkgs[!ok], collapse = ", "))
}
