# ######################################################################
# BMP15 genetic diversity in Muturu, Kuri and White Fulani cattle - MSc dissertation (Chapter 4)
# Complete analysis in one R script (all steps, in order).
#
# HOW TO USE
#   1. Put this file in the project folder that contains the data/ folder
#      (the same folder as BMP15_pipeline.Rproj) and open that .Rproj in RStudio.
#   2. Run the whole file:  source("BMP15_dissertation_analysis.R")   (or Ctrl + Shift + Enter)
#   Results are written to R_outputs/; R and package versions to logs/.
#
# Each step runs inside local({ ... }) so that its variables do not
# affect the next step.
# ######################################################################

# ======================================================================
# 0. PACKAGES (installed only if missing; needs internet the first time)
# ======================================================================
# ======================================================================
# install_packages.R  –  run ONCE on a new computer (needs internet)
# ======================================================================
cran <- c("ape", "pegas", "rentrez")
for (p in cran) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (p in c("Biostrings", "pwalign")) {
  if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, update = FALSE, ask = FALSE)
}
cat("\nPackage check:\n")
for (p in c(cran, "Biostrings", "pwalign")) cat(sprintf("  %-10s %s\n", p, ifelse(requireNamespace(p, quietly = TRUE), "OK", "MISSING")))

# ======================================================================
# 0. SETTINGS
# ======================================================================
# ======================================================================
# config.R  –  settings used by the BMP15 pipeline
# Change values here only; the scripts read them from this file.
# ======================================================================

cfg <- list(
  # Input data -----------------------------------------------------------
  ab1_dir        = "data/raw/ab1",                       # raw Sanger traces (.ab1), never edited
  alignment      = "data/raw/BMP15_alignment.fasta",     # Geneious export: reads mapped to the CDS
  review_table   = "data/curated/02_chromatogram_review.csv",  # curated genotype calls (manual review)
  reference_cds  = "data/reference/XM_027535004.1_BMP15_CDS.fasta",
  expected_qc    = "data/expected/expected_trace_qc.csv",

  # Amplicon and genome coordinates -------------------------------------
  amplicon_start = 606,          # first CDS base covered by the primers (c.606)
  amplicon_end   = 976,          # last CDS base (c.976)
  genome_chrom   = "X",          # chromosome name used in the VCF
  genome_offset  = 88809990,     # ARS-UCD2.0 position of c.606 (from BLASTn)
  assembly       = "GCA_002263795.4",  # ARS-UCD2.0

  # Peak-calling settings (same as Chapter 3) ---------------------------
  secondary_peak = 0.25,         # a second base counts if its peak is >= 25% of the highest
  c699_C_max     = 0.15,         # A/(A+C) below this  -> C only
  c699_A_min     = 0.90,         # A/(A+C) above this  -> A only (between -> C and A)
  del_support    = 0.60,         # share of positions 3' of c.787 explained by a 4-bp shift -> carrier
  anchor_c699    = "TGGCATGG",       # c.691-c.698; the next base is c.699
  anchor_c787    = "CCCTAAAGGCCTGA", # c.774-c.787; ends at c.787
  tolerance      = 0.002         # allowed difference when checking results against the dissertation
)

dir.create("R_outputs", showWarnings = FALSE)
dir.create("logs", showWarnings = FALSE)

# ======================================================================
# STEP 1: 00_trace_qc
# ======================================================================
local({
# ======================================================================
# 00_trace_qc.R  –  from raw .ab1 chromatograms to the measurements
#                   behind every genotype call
#
# For each trace this script:
#   1. reads the base calls, quality values and the four dye signals;
#   2. counts read length and bases with Phred quality >= 20;
#   3. at c.699 measures A / (A + C) from the peak heights;
#   4. 3' of c.787 measures how often the two highest peaks are the
#      normal base AND the base four positions downstream (the signature
#      of the c.788_791delAAGA deletion on one template);
#   5. turns these into automatic calls and compares them with the
#      curated calls in data/curated/02_chromatogram_review.csv.
# Uses base R only (no extra packages).
# ======================================================================
# (settings already loaded above)

# ---- 1. Minimal ABIF (.ab1) reader ---------------------------------
read_abif <- function(path) {
  b <- readBin(path, "raw", n = file.info(path)$size)
  if (rawToChar(b[1:4]) != "ABIF") stop(path, " is not an ABIF file")
  i32 <- function(off) readBin(b[(off + 1):(off + 4)], "integer", size = 4, endian = "big")
  i16 <- function(off) readBin(b[(off + 1):(off + 2)], "integer", size = 2, endian = "big", signed = TRUE)
  n_entries <- i32(6 + 12)
  dir_off   <- i32(6 + 20)
  tags <- list()
  for (k in 0:(n_entries - 1)) {
    e    <- dir_off + k * 28
    name <- rawToChar(b[(e + 1):(e + 4)])
    num  <- i32(e + 4)
    type <- i16(e + 8)
    cnt  <- i32(e + 12)
    size <- i32(e + 16)
    if (size <= 0) next
    off  <- if (size <= 4) e + 20 else i32(e + 20)
    tags[[paste0(name, ".", num)]] <- list(type = type, cnt = cnt, raw = b[(off + 1):(off + size)])
  }
  get_short <- function(key) readBin(tags[[key]]$raw, "integer", size = 2, n = tags[[key]]$cnt,
                                     endian = "big", signed = TRUE)
  order <- strsplit(rawToChar(tags[["FWO_.1"]]$raw[1:4]), "")[[1]]
  trace <- lapply(9:12, function(n) get_short(paste0("DATA.", n)))
  names(trace) <- order
  list(seq   = rawToChar(tags[["PBAS.2"]]$raw),
       qual  = as.integer(tags[["PCON.2"]]$raw),
       ploc  = get_short("PLOC.2") + 1L,            # +1: R vectors start at 1
       trace = trace)
}

# ---- 2. Reference coding sequence ----------------------------------
ref_lines <- readLines(cfg$reference_cds)
ref <- strsplit(paste(ref_lines[!grepl("^>", ref_lines)], collapse = ""), "")[[1]]

peak_heights <- function(r, p) {
  loc <- r$ploc[p]
  sapply(c("A", "C", "G", "T"), function(bse) r$trace[[bse]][loc])
}
peak_set <- function(h) names(h)[h >= cfg$secondary_peak * max(max(h), 1)]

# ---- 3. Measure every trace ----------------------------------------
files <- sort(list.files(cfg$ab1_dir, pattern = "\\.ab1$", full.names = TRUE))
cat("Trace files found:", length(files), "\n")
qc <- do.call(rbind, lapply(files, function(f) {
  r <- read_abif(f)
  s <- r$seq
  sample <- sub("\\.ab1$", "", basename(f))
  i699 <- regexpr(cfg$anchor_c699, s, fixed = TRUE)[1]
  i787 <- regexpr(cfg$anchor_c787, s, fixed = TRUE)[1]
  p787 <- if (i787 > 0) i787 + nchar(cfg$anchor_c787) - 1 else NA
  p699 <- if (i699 > 0) i699 + nchar(cfg$anchor_c699) else if (!is.na(p787)) p787 - 88 else NA
  a_frac <- NA
  if (!is.na(p699)) {
    h <- peak_heights(r, p699)
    a_frac <- round(h[["A"]] / max(h[["A"]] + h[["C"]], 1), 3)
  }
  support <- NA; n_pos <- 0
  if (!is.na(p787)) {
    hits <- 0
    for (cc in 790:899) {
      p <- p787 + (cc - 787)
      if (p > length(r$ploc)) break
      P <- peak_set(peak_heights(r, p))
      n_pos <- n_pos + 1
      hits <- hits + (ref[cc] %in% P && ref[cc + 4] %in% P)
    }
    if (n_pos > 0) support <- round(hits / n_pos, 3)
  }
  data.frame(Sample = sample, Read_length = nchar(s), Q20_bases = sum(r$qual >= 20),
             c699_A_fraction = a_frac, del4_shift_support = support, n_positions = n_pos)
}))

# ---- 4. Automatic calls and comparison with the curated calls ------
qc$auto_c699 <- ifelse(is.na(qc$c699_A_fraction), "not scored",
                ifelse(qc$c699_A_fraction < cfg$c699_C_max, "C",
                ifelse(qc$c699_A_fraction > cfg$c699_A_min, "A", "C/A")))
qc$auto_del <- ifelse(is.na(qc$del4_shift_support), "not scored",
               ifelse(qc$del4_shift_support >= cfg$del_support, "Carrier", "Non-carrier"))
review <- read.csv(cfg$review_table, stringsAsFactors = FALSE)
qc <- merge(qc, review[, c("Sample", "Use_in_analysis", "c699_call", "del4_status")], by = "Sample")
qc$c699_agrees <- ifelse(qc$Use_in_analysis != "Yes", "excluded",
                  ifelse(sub(" .*", "", qc$c699_call) == qc$auto_c699, "yes",
                  ifelse(qc$c699_call == "uncertain", "review (curated as uncertain)", "NO")))
qc$del_agrees <- ifelse(qc$Use_in_analysis != "Yes", "excluded",
                 ifelse(qc$del4_status == qc$auto_del, "yes", "NO"))

cat("\n===== Raw-trace QC and automatic calls =====\n")
print(qc[, c("Sample", "Read_length", "Q20_bases", "c699_A_fraction", "auto_c699", "c699_call",
             "del4_shift_support", "auto_del", "del4_status")], row.names = FALSE)
write.csv(qc, "R_outputs/00_trace_qc.csv", row.names = FALSE)

used <- qc[qc$Use_in_analysis == "Yes", ]
cat("\nAnimals analysed:", nrow(used),
    "| c.699 calls agreeing with curated table:", sum(used$c699_agrees == "yes"),
    "| deletion calls agreeing:", sum(used$del_agrees == "yes"), "\n")
if (any(used$c699_agrees == "NO" | used$del_agrees == "NO"))
  warning("Some automatic calls disagree with the curated table - see R_outputs/00_trace_qc.csv")

# ---- 5. Same numbers as when the analysis was first done? ----------
if (file.exists(cfg$expected_qc)) {
  ex <- read.csv(cfg$expected_qc, stringsAsFactors = FALSE)
  m <- merge(qc, ex, by = "Sample", suffixes = c("", ".expected"))
  same <- function(a, b) (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & abs(a - b) <= 0.002)
  diffs <- m[!(m$Read_length == m$Read_length.expected & m$Q20_bases == m$Q20_bases.expected &
               same(m$c699_A_fraction, m$c699_A_fraction.expected) &
               same(m$del4_shift_support, m$del4_shift_support.expected)), "Sample"]
  cat("Traces reproducing the original measurements exactly:", nrow(m) - length(diffs), "of", nrow(m), "\n")
  if (length(diffs)) cat("Differences in:", paste(diffs, collapse = ", "), "\n")
}
cat("\nDone. Results saved in R_outputs/00_trace_qc.csv\n")
})
while (!is.null(dev.list())) dev.off()

# ======================================================================
# STEP 2: 01_load_alignment
# ======================================================================
local({
# =====================================================================
# BMP15 dissertation - SCRIPT 01: Load alignment and build sample table
# =====================================================================

# ---- 1. Load the package we need ------------------------------------
library(ape)

# ---- 2. Tell R where your files are ---------------------------------
# R needs forward slashes ( / ), not backslashes ( \ )
dir.create("R_outputs", showWarnings = FALSE)   # folder for all results

# ---- 3. Read the alignment ------------------------------------------
# Rename your Geneious FASTA export to this simple name first
aln <- read.dna("data/raw/BMP15_alignment.fasta", format = "fasta", as.character = TRUE)
aln <- toupper(aln)
dim(aln)          # should print: 34  1213  (34 sequences incl. reference)

ref     <- aln[1, ]                         # row 1 = reference sequence
samples <- aln[-1, , drop = FALSE]          # rows 2-34 = your animals
rownames(samples) <- sub("\\.ab1$", "", rownames(samples))

# Convert alignment columns to CDS positions (c.1 = A of ATG)
cds_pos <- cumsum(ref != "-")

# ---- 4. Build the sample information table from the names -----------
parts <- do.call(rbind, strsplit(rownames(samples), "_"))
meta <- data.frame(
  Sample = rownames(samples),
  Breed  = c(MT = "Muturu", KR = "Kuri", WF = "White Fulani")[parts[, 1]],
  ID     = parts[, 2],
  Sex    = c(F = "Female", M = "Male", CM = "Castrated male")[parts[, 3]],
  Type   = c(T = "Taurine", Z = "Zebu")[parts[, 4]],
  row.names = NULL
)
# BMP15 is on the X chromosome: females carry 2 copies, males 1
meta$X_copies <- ifelse(meta$Sex == "Female", 2, 1)

# ---- 5. Where does each read start and end? -------------------------
span <- t(apply(samples, 1, function(s) range(which(s != "-"))))
meta$Start  <- cds_pos[span[, 1]]
meta$End    <- cds_pos[span[, 2]]
meta$Length <- meta$End - meta$Start + 1

# Base each animal carries at c.699 ("." = not covered)
col699 <- which(cds_pos == 699 & ref != "-")
meta$c699 <- ifelse(span[, 1] <= col699 & span[, 2] >= col699,
                    samples[, col699], ".")

print(meta)
write.csv(meta, "R_outputs/01_sample_table.csv", row.names = FALSE)

# ---- 6. Quick summaries ----------------------------------------------
cat("\nAnimals per breed and sex:\n")
print(table(meta$Breed, meta$Sex))

cat("\nGenotype at c.699 versus where the read ends:\n")
print(table(c699 = meta$c699, Read_ends_by_790 = meta$End <= 790))

# ---- 7. Coverage plot: how many animals have data at each position --
keep <- ref != "-"
cov_by_breed <- sapply(c("Muturu", "Kuri", "White Fulani"), function(b) {
  rows <- meta$Breed == b
  sapply(which(keep), function(j) sum(span[rows, 1] <= j & span[rows, 2] >= j))
})
x <- cds_pos[keep]
in_amp <- x >= 600 & x <= 980

png("R_outputs/01_coverage_plot.png", width = 1800, height = 900, res = 200)
matplot(x[in_amp], cov_by_breed[in_amp, ], type = "l", lty = 1, lwd = 2,
        col = c("darkgreen", "darkorange", "steelblue"),
        xlab = "BMP15 CDS position (c.)", ylab = "Number of animals with sequence",
        main = "Sequence coverage across the BMP15 exon 2 amplicon")
abline(v = c(668, 787), lty = 2, col = "grey40")
legend("topright", legend = colnames(cov_by_breed),
       col = c("darkgreen", "darkorange", "steelblue"), lwd = 2, bty = "n")
dev.off()

cat("\nDone. Results saved in R_outputs.\n")
})
while (!is.null(dev.list())) dev.off()

# ======================================================================
# STEP 3: 02_variant_tables
# ======================================================================
local({
# =====================================================================
# BMP15 dissertation - SCRIPT 02: Clean genotype table and first results
# =====================================================================

# ---- 1. Go to your folder and read the chromatogram review table ----
dir.create("R_outputs", showWarnings = FALSE)

review <- read.csv("data/curated/02_chromatogram_review.csv", stringsAsFactors = FALSE)
cat("Animals in review table:", nrow(review), "\n")      # should be 34

# Keep only the animals with usable sequence
dat <- review[review$Use_in_analysis == "Yes", ]
cat("Animals used in analysis:", nrow(dat), "\n")         # should be 29

# Put breeds in a fixed order for all tables
dat$Breed <- factor(dat$Breed, levels = c("Muturu", "Kuri", "White Fulani"))

# ---- 2. Simplify the c.699 call into three groups -------------------
dat$c699_group <- "Uncertain"
dat$c699_group[dat$c699_call == "C"] <- "C only"
dat$c699_group[dat$c699_call == "A"] <- "A only"
dat$c699_group[grepl("^C/A", dat$c699_call)] <- "C and A"
dat$c699_group <- factor(dat$c699_group,
                         levels = c("C only", "C and A", "A only", "Uncertain"))

# Deletion status as TRUE / FALSE
dat$del_carrier <- dat$del4_status == "Carrier"

# ---- TABLE 4.1: Animals analysed, by breed and sex -------------------
cat("\n===== TABLE 4.1: Animals analysed by breed and sex =====\n")
t41 <- table(dat$Breed, dat$Sex)
t41 <- addmargins(t41)
print(t41)
write.csv(as.data.frame.matrix(t41), "R_outputs/Table4.1_samples_by_breed_sex.csv")

# ---- TABLE 4.4: The variants identified -----------------------------
cat("\n===== TABLE 4.4: BMP15 variants identified (exon 2) =====\n")
t42 <- data.frame(
  Variant        = c("c.699C>A", "c.788_791delAAGA", "c.859G>A"),
  Type           = c("SNP (transversion)", "Deletion (4 bp)", "SNP (transition)"),
  Codon_change   = c("GGC>GGA (codon 233)", "Frameshift from codon 263", "GTT>ATT (codon 287)"),
  Protein_effect = c("p.Gly233= (synonymous)",
                     "p.Lys263Serfs*11 (premature stop; protein 272 aa instead of 403 aa)",
                     "p.Val287Ile (missense)"),
  Status         = c("Confirmed", "Putative (two-template pattern)", "Unconfirmed - not analysed"),
  stringsAsFactors = FALSE
)
print(t42, right = FALSE)
write.csv(t42, "R_outputs/Table4.2_variants_identified.csv", row.names = FALSE)

# ---- TABLE 4.5: c.699 calls by breed ---------------------------------
cat("\n===== TABLE 4.5: c.699C>A peak pattern by breed (counts) =====\n")
t43 <- table(dat$Breed, dat$c699_group)
print(addmargins(t43))
cat("\nSame table as row percentages:\n")
print(round(prop.table(t43, 1) * 100, 1))
write.csv(as.data.frame.matrix(addmargins(t43)), "R_outputs/Table4.3_c699_by_breed.csv")

# ---- TABLE 4.11: Deletion carrier frequency by breed, with 95% CI ----
cat("\n===== TABLE 4.11: c.788_791del carrier frequency by breed =====\n")
t44 <- do.call(rbind, lapply(levels(dat$Breed), function(b) {
  x <- dat$del_carrier[dat$Breed == b]
  ci <- binom.test(sum(x), length(x))$conf.int      # exact 95% CI
  data.frame(Breed = b,
             N = length(x),
             Carriers = sum(x),
             Non_carriers = sum(!x),
             Carrier_freq = round(mean(x), 3),
             CI95_lower = round(ci[1], 3),
             CI95_upper = round(ci[2], 3))
}))
# Add an "All breeds" row
x  <- dat$del_carrier
ci <- binom.test(sum(x), length(x))$conf.int
t44 <- rbind(t44, data.frame(Breed = "All breeds", N = length(x), Carriers = sum(x),
                             Non_carriers = sum(!x), Carrier_freq = round(mean(x), 3),
                             CI95_lower = round(ci[1], 3), CI95_upper = round(ci[2], 3)))
print(t44, row.names = FALSE)
write.csv(t44, "R_outputs/Table4.4_deletion_carrier_frequency.csv", row.names = FALSE)

# ---- LD CHECK (Section 4.6.1): Do the two variants go together? ---------------------
cat("\n===== LD CHECK (Section 4.6.1): c.699 pattern versus deletion status =====\n")
t45 <- table(c699 = dat$c699_group, Deletion = dat$del4_status)
print(t45)
write.csv(as.data.frame.matrix(t45), "R_outputs/Table4.5_c699_vs_deletion.csv")

# ---- TEST: Do the breeds differ in deletion carrier frequency? ------
cat("\n===== Fisher's exact tests: carrier frequency between breeds =====\n")
mk <- dat[dat$Breed %in% c("Muturu", "Kuri"), ]
mk$Breed <- droplevels(mk$Breed)
tab_mk <- table(mk$Breed, mk$del_carrier)
colnames(tab_mk) <- c("Non-carrier", "Carrier")
print(tab_mk)
f_mk <- fisher.test(tab_mk)
cat("Muturu vs Kuri: p-value =", signif(f_mk$p.value, 3), "\n")

tab_all <- table(dat$Breed, dat$del_carrier)
f_all <- fisher.test(tab_all)
cat("All three breeds: p-value =", signif(f_all$p.value, 3), "\n")
cat("(p < 0.05 means the breeds differ significantly)\n")

tests <- data.frame(Comparison = c("Muturu vs Kuri", "Muturu vs Kuri vs White Fulani"),
                    Test = "Fisher's exact test",
                    P_value = signif(c(f_mk$p.value, f_all$p.value), 3))
write.csv(tests, "R_outputs/Table4.6_breed_comparison_tests.csv", row.names = FALSE)

# ---- FIGURE 4.x: Carrier frequency bar chart --------------------------
png("R_outputs/Fig_deletion_carrier_frequency.png", width = 1600, height = 1200, res = 200)
bp <- barplot(t44$Carrier_freq[1:3], names.arg = t44$Breed[1:3],
              ylim = c(0, 1), col = c("darkgreen", "darkorange", "steelblue"),
              ylab = "Proportion of animals carrying c.788_791del",
              main = "BMP15 c.788_791del carrier frequency by breed")
arrows(bp, t44$CI95_lower[1:3], bp, t44$CI95_upper[1:3],
       angle = 90, code = 3, length = 0.08, lwd = 1.5)
text(bp, 0.05, paste0("n = ", t44$N[1:3]), col = "white", font = 2)
dev.off()

# ---- Save the clean genotype table for later scripts ---------------
write.csv(dat[, c("Sample", "Breed", "ID", "Sex", "Type",
                  "c699_call", "c699_group", "del4_status")],
          "R_outputs/02_clean_genotype_table.csv", row.names = FALSE)

cat("\nDone. Results saved in R_outputs.\n")
})
while (!is.null(dev.list())) dev.off()

# ======================================================================
# STEP 4: 03_protein_effect
# ======================================================================
local({
# =====================================================================
# BMP15 dissertation - SCRIPT 03: Effect of the variants on the protein
# =====================================================================

library(ape)
dir.create("R_outputs", showWarnings = FALSE)

# ---- 1. Get the reference BMP15 coding sequence (CDS) ---------------
aln <- read.dna("data/raw/BMP15_alignment.fasta", format = "fasta", as.character = TRUE)
ref <- toupper(aln[1, ])
ref <- paste(ref[ref != "-"], collapse = "")      # remove alignment gaps
cat("Reference CDS length:", nchar(ref), "bases\n")   # should be 1212

# ---- 2. A simple translation function (DNA -> protein) -------------
bases <- c("T", "C", "A", "G")
codon_names <- paste0(rep(bases, each = 16), rep(rep(bases, each = 4), 4), rep(bases, 16))
amino <- strsplit("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "")[[1]]
genetic_code <- setNames(amino, codon_names)

translate <- function(dna) {
  n <- nchar(dna)
  codons <- substring(dna, seq(1, n - 2, 3), seq(3, n, 3))
  protein <- paste(genetic_code[codons], collapse = "")
  sub("\\*.*$", "", protein)          # stop at the first stop codon (*)
}

# ---- 3. Make the three sequences: normal, c.699A, and the deletion ---
wt_dna <- ref
snp_dna <- ref; substr(snp_dna, 699, 699) <- "A"                      # c.699C>A
del_dna <- paste0(substr(ref, 1, 787), substr(ref, 792, nchar(ref)))  # c.788_791del

cat("Bases deleted at c.788-791:", substr(ref, 788, 791), "\n")      # should be AAGA

wt_prot  <- translate(wt_dna)
snp_prot <- translate(snp_dna)
del_prot <- translate(del_dna)

# ---- 4. Where is the cleavage site that releases the mature protein? --
# BMP15 is made as a long precursor. An enzyme (furin) cuts it after the
# motif R-R-A-R; the part after the cut is the active "mature" BMP15.
cut_site <- regexpr("RRAR", wt_prot)[1]
mature_start <- cut_site + 4
cat("Furin cleavage motif RRAR at amino acids", cut_site, "-", cut_site + 3, "\n")
cat("Mature BMP15 domain: amino acids", mature_start, "-", nchar(wt_prot), "\n")

# ---- TABLE 4.7: Protein consequences ---------------------------------
t47 <- data.frame(
  Sequence = c("Reference (wild type)", "c.699C>A", "c.788_791delAAGA"),
  Protein_length_aa = c(nchar(wt_prot), nchar(snp_prot), nchar(del_prot)),
  Identical_to_reference = c("-", ifelse(snp_prot == wt_prot, "Yes", "No"), "No"),
  Contains_cleavage_site_RRAR = c(grepl("RRAR", wt_prot), grepl("RRAR", snp_prot), grepl("RRAR", del_prot)),
  Mature_domain_present = c("Yes (aa 279-403)",
                            ifelse(snp_prot == wt_prot, "Yes (unchanged)", "Check"),
                            "No (protein ends before cleavage site)"),
  stringsAsFactors = FALSE
)
cat("\n===== TABLE 4.7: Predicted effect of variants on BMP15 protein =====\n")
print(t47, right = FALSE)
write.csv(t47, "R_outputs/Table4.7_protein_consequences.csv", row.names = FALSE)

# ---- 5. Show the wild-type vs deletion protein around the change -----
cat("\n===== Amino acids 250-300: wild type vs deletion =====\n")
cat("Wild type :", substr(wt_prot, 250, 300), "\n")
cat("Deletion  :", substr(del_prot, 250, 300), "*(stop)\n")
cat("            The two sequences differ from amino acid 263 onward.\n")

# ---- 6. Save protein sequences for the web tools (next step) --------
writeLines(c(">BMP15_wild_type", wt_prot,
             ">BMP15_c.788_791del", del_prot),
           "R_outputs/03_BMP15_proteins.fasta")

# ---- FIGURE: Protein diagram, wild type vs deletion ------------------
png("R_outputs/Fig_protein_diagram.png", width = 2000, height = 900, res = 200)
L <- nchar(wt_prot)
plot(NA, xlim = c(0, L + 20), ylim = c(0, 3), axes = FALSE,
     xlab = "Amino acid position", ylab = "",
     main = "Predicted BMP15 protein: wild type vs c.788_791del")
axis(1, at = seq(0, 400, 50))
# wild type
rect(1, 2, cut_site - 1, 2.5, col = "grey80")
rect(cut_site, 2, cut_site + 3, 2.5, col = "red")
rect(mature_start, 2, L, 2.5, col = "darkgreen")
text(0, 2.25, "Wild type", pos = 2, xpd = TRUE, cex = 0.8)
text(cut_site / 2, 2.7, "Signal peptide + pro-region", cex = 0.8)
text((mature_start + L) / 2, 2.7, "Mature domain (active)", cex = 0.8)
text(cut_site + 1, 1.85, "RRAR", col = "red", cex = 0.7)
# deletion
D <- nchar(del_prot)
rect(1, 0.7, 262, 1.2, col = "grey80")
rect(263, 0.7, D, 1.2, col = "orange")
text(0, 0.95, "Deletion", pos = 2, xpd = TRUE, cex = 0.8)
text(D + 5, 0.95, paste0("stop at ", D + 1), pos = 4, cex = 0.8)
legend("bottomright", fill = c("orange"), legend = "Altered amino acids (263-272)",
       bty = "n", cex = 0.8)
dev.off()

cat("\nDone. Results saved in R_outputs.\n")
})
while (!is.null(dev.list())) dev.off()

# ======================================================================
# STEP 5: 04_protparam_table
# ======================================================================
local({
# =====================================================================
# BMP15 dissertation - SCRIPT 04: ProtParam results table (Table 4.8)
# =====================================================================

dir.create("R_outputs", showWarnings = FALSE)

# ---- 1. Type in the ProtParam values ---------------------------------
t48 <- data.frame(
  Property = c("Number of amino acids",
               "Molecular weight (Da)",
               "Theoretical pI",
               "Instability index",
               "Aliphatic index",
               "GRAVY"),
  Wild_type   = c(403, 46114.50,  9.61, 48.96, 94.29, -0.281),
  c788_791del = c(272, 31447.91, 10.71, 43.13, 97.43, -0.277)
)

# ---- 2. How much did each value change? ------------------------------
t48$Difference     <- round(t48$c788_791del - t48$Wild_type, 3)
t48$Percent_change <- round(100 * t48$Difference / abs(t48$Wild_type), 1)

# ---- 3. Plain-language interpretation --------------------------------
t48$Interpretation <- c(
  "Deletion protein is 131 aa (32.5%) shorter",
  "Deletion protein is about 14.7 kDa lighter",
  "Both basic (pI > 7); deletion protein more basic",
  "Both > 40, so both predicted unstable in a test tube",
  "Similar thermal stability (aliphatic side chains)",
  "Both negative, so both hydrophilic (water-soluble)"
)

cat("\n===== TABLE 4.8: Physicochemical properties (ProtParam) =====\n")
print(t48, right = FALSE, row.names = FALSE)
write.csv(t48, "R_outputs/Table4.8_protparam.csv", row.names = FALSE)

cat("\nDone. Results saved in R_outputs.\n")
})
while (!is.null(dev.list())) dev.off()

# ======================================================================
# STEP 6: 05_structure_table
# ======================================================================
local({
# =====================================================================
# BMP15 dissertation - SCRIPT 05: SWISS-MODEL results (Table 4.9)
# =====================================================================

dir.create("R_outputs", showWarnings = FALSE)

# NA = value not shown by SWISS-MODEL for that model
t49 <- data.frame(
  Protein = c(rep("Wild type (403 aa)", 6), rep("c.788_791del (272 aa)", 5)),
  Model   = c("01", "02", "03", "04", "05", "06",
              "01", "02", "03", "04", "05"),
  Template = c("A0A3Q7W7G9.1.A", "5vqf.1.A", "5hly.1.A", "5ntu.1.A", "9r3s.1.B", "6umx.1.A",
               "A0A7J7UQ79.1.A", "6sf2.1.F", "9r3s.1.B", "5hlz.1.A", "5ntu.1.A"),
  Template_protein = c("BMP15 (AlphaFold model, organism not stated)",
                       "Pro-TGF-beta1", "Pro-activin A", "Pro-myostatin (GDF8)",
                       "Pro-TGF-beta1", "Pro/latent myostatin (GDF8)",
                       "BMP15, Myotis myotis (AlphaFold model)",
                       "BMP9 (GDF2) prodomain complex", "Pro-TGF-beta1",
                       "Pro-activin A", "Pro-myostatin (GDF8)"),
  Template_type = c("Predicted (AlphaFold)", rep("Experimental (X-ray)", 5),
                    "Predicted (AlphaFold)", rep("Experimental (X-ray)", 4)),
  Seq_identity_pct = c(83.04, 21.36, 20.63, 21.05, 21.02, 19.93,
                       79.04, 14.12, 17.65, 15.24, 12.00),
  Oligo_state = c("Monomer", "Homo-dimer", "Homo-dimer", "Homo-dimer", "Monomer", "Monomer",
                  "Monomer", "Monomer", "Monomer", "Monomer", "Monomer"),
  GMQE = c(0.70, 0.36, 0.41, 0.39, 0.35, 0.39,
           0.66, 0.33, 0.27, 0.24, 0.25),
  QMEANDisCo = c(NA, 0.42, 0.43, 0.41, 0.41, 0.41,
                 NA, 0.48, 0.44, 0.40, 0.40),
  Region_modelled = c("Full length (1-403)", rep("Prodomain + mature domain", 5),
                      "Full length (1-272)", rep("Prodomain only", 4)),
  stringsAsFactors = FALSE
)

cat("\n===== TABLE 4.9: SWISS-MODEL homology modelling of BMP15 =====\n")
print(t49, right = FALSE, row.names = FALSE)
write.csv(t49, "R_outputs/Table4.9_swissmodel.csv", row.names = FALSE)

# ---- Summary: wild type vs deletion -------------------------------
cat("\n===== Summary: wild type vs deletion =====\n")
summ <- do.call(rbind, lapply(split(t49, t49$Protein), function(d) {
  exp <- d[d$Template_type != "Predicted (AlphaFold)", ]
  data.frame(Protein = d$Protein[1],
             Models = nrow(d),
             Best_GMQE = max(d$GMQE),
             Dimer_models = sum(d$Oligo_state == "Homo-dimer"),
             Identity_experimental_templates = paste0(min(exp$Seq_identity_pct), "-",
                                                      max(exp$Seq_identity_pct), "%"),
             Mean_QMEANDisCo_experimental = round(mean(exp$QMEANDisCo), 3))
}))
print(summ, row.names = FALSE)
write.csv(summ, "R_outputs/Table4.9b_structure_summary.csv", row.names = FALSE)

cat("\nDone. Results saved in R_outputs.\n")
})
while (!is.null(dev.list())) dev.off()

# ======================================================================
# STEP 7: 06_diversity
# ======================================================================
local({
# =====================================================================
# BMP15 dissertation - SCRIPT 06: Genetic diversity within breeds
# =====================================================================
# BMP15 is on the X chromosome, so:
#  * Genotype-based indices (Ho, He, Fis, HWE) use FEMALES only (2 copies)
#  * Allele frequencies also use males, counted as 1 copy (hemizygous)
#  * Males showing two peaks (MT_64, KR_142) and WF_18 (uncertain call)
#    are left out of allele counting
# =====================================================================

dir.create("R_outputs", showWarnings = FALSE)

dat <- read.csv("R_outputs/02_clean_genotype_table.csv", stringsAsFactors = FALSE)
breeds <- c("Muturu", "Kuri", "White Fulani")

# ---- 1. Turn the c.699 peak pattern into a genotype -----------------
dat$Genotype <- NA
dat$Genotype[dat$c699_group == "C only"]  <- "CC"
dat$Genotype[dat$c699_group == "C and A"] <- "CA"
dat$Genotype[dat$c699_group == "A only"]  <- "AA"

# ---- 2. Hardy-Weinberg exact test (Wigginton et al., 2005) ----------
hwe_exact <- function(n_het, n_hom1, n_hom2) {
  n <- n_het + n_hom1 + n_hom2
  rare <- 2 * min(n_hom1, n_hom2) + n_het
  common <- 2 * n - rare
  if (rare == 0) return(1)
  probs <- numeric(rare + 1)                       # index = hets + 1
  mid <- floor(rare * common / (2 * n))
  if ((mid %% 2) != (rare %% 2)) mid <- mid + 1
  probs[mid + 1] <- 1
  hr <- (rare - mid) / 2; hc <- n - mid - hr; h <- mid
  while (h >= 2) {
    probs[h - 1] <- probs[h + 1] * h * (h - 1) / (4 * (hr + 1) * (hc + 1))
    h <- h - 2; hr <- hr + 1; hc <- hc + 1
  }
  hr <- (rare - mid) / 2; hc <- n - mid - hr; h <- mid
  while (h <= rare - 2) {
    probs[h + 3] <- probs[h + 1] * 4 * hr * hc / ((h + 2) * (h + 1))
    h <- h + 2; hr <- hr - 1; hc <- hc - 1
  }
  probs <- probs / sum(probs)
  min(1, sum(probs[probs <= probs[n_het + 1] + 1e-12]))
}

# ---- 3. Diversity indices from FEMALE genotypes ---------------------
diversity <- function(g) {
  g <- g[!is.na(g)]
  n <- length(g)
  CC <- sum(g == "CC"); CA <- sum(g == "CA"); AA <- sum(g == "AA")
  pA <- (CA + 2 * AA) / (2 * n); pC <- 1 - pA
  Ho  <- CA / n
  He  <- 1 - (pA^2 + pC^2)                         # expected heterozygosity
  uHe <- He * 2 * n / (2 * n - 1)                  # unbiased (Nei 1978)
  data.frame(N_females = n, CC = CC, CA = CA, AA = AA,
             Freq_C = round(pC, 3), Freq_A = round(pA, 3),
             Ho = round(Ho, 3), He = round(He, 3), uHe = round(uHe, 3),
             Fis = ifelse(uHe > 0, round(1 - Ho / uHe, 3), NA),
             PIC = round(1 - (pA^2 + pC^2) - 2 * pA^2 * pC^2, 3),
             Ne = round(1 / (pA^2 + pC^2), 3),
             HWE_p = round(hwe_exact(CA, CC, AA), 4))
}

fem <- dat[dat$Sex == "Female", ]
t410 <- do.call(rbind, lapply(breeds, function(b) cbind(Breed = b, diversity(fem$Genotype[fem$Breed == b]))))
t410 <- rbind(t410, cbind(Breed = "All breeds", diversity(fem$Genotype)))
t410$Note <- ifelse(t410$N_females < 5, "Too few animals - descriptive only", "")

cat("\n===== TABLE 4.12: c.699C>A genotype frequencies and diversity (females) =====\n")
print(t410, row.names = FALSE)
write.csv(t410, "R_outputs/Table4.10_diversity_females.csv", row.names = FALSE)

# ---- 4. Sex-aware allele frequencies (females x2, males x1) ---------
excl <- c("MT_64_M_T", "KR_142_CM_T", "WF_18_F_Z")
use  <- dat[!(dat$Sample %in% excl) & !is.na(dat$Genotype), ]
count_alleles <- function(d) {
  fe <- d[d$Sex == "Female", ]; ma <- d[d$Sex != "Female", ]
  nC <- 2 * sum(fe$Genotype == "CC") + sum(fe$Genotype == "CA") + sum(ma$Genotype == "CC")
  nA <- 2 * sum(fe$Genotype == "AA") + sum(fe$Genotype == "CA") + sum(ma$Genotype == "AA")
  data.frame(Animals = nrow(d), Females = nrow(fe), Males = nrow(ma),
             X_chromosomes = nC + nA, Count_C = nC, Count_A = nA,
             Freq_C = round(nC / (nC + nA), 3), Freq_A = round(nA / (nC + nA), 3))
}
t411 <- do.call(rbind, lapply(breeds, function(b) cbind(Breed = b, count_alleles(use[use$Breed == b, ]))))
t411 <- rbind(t411, cbind(Breed = "All breeds", count_alleles(use)))

cat("\n===== TABLE 4.13: c.699 allele frequencies, sex-aware (X-linked) =====\n")
print(t411, row.names = FALSE)
write.csv(t411, "R_outputs/Table4.11_allele_freq_sex_aware.csv", row.names = FALSE)

# ---- 5. Figure: allele frequencies by breed --------------------------
png("R_outputs/Fig_c699_allele_frequency.png", width = 1600, height = 1200, res = 200)
m <- t(as.matrix(t411[1:3, c("Freq_C", "Freq_A")]))
colnames(m) <- t411$Breed[1:3]
barplot(m, beside = FALSE, col = c("grey70", "darkorange"), ylim = c(0, 1.15),
        ylab = "Allele frequency", main = "BMP15 c.699C>A allele frequency by breed (X-linked counting)")
legend("top", legend = c("C (reference)", "A (variant)"), fill = c("grey70", "darkorange"),
       horiz = TRUE, bty = "n")
text(seq(0.7, by = 1.2, length.out = 3), 1.05, paste0("n = ", t411$X_chromosomes[1:3], " X"))
dev.off()

cat("\nDone. Results saved in R_outputs.\n")
})
while (!is.null(dev.list())) dev.off()

# ======================================================================
# STEP 8: 07_differentiation
# ======================================================================
local({
# =====================================================================
# BMP15 dissertation - SCRIPT 07: Genetic differentiation between breeds
# Uses the sex-aware c.699 allele counts from Script 06 (Table 4.13)
# =====================================================================

dir.create("R_outputs", showWarnings = FALSE)

af <- read.csv("R_outputs/Table4.11_allele_freq_sex_aware.csv", stringsAsFactors = FALSE)
af <- af[af$Breed != "All breeds", ]
rownames(af) <- af$Breed
print(af[, c("Breed", "X_chromosomes", "Count_C", "Count_A", "Freq_A")], row.names = FALSE)

# ---- 1. Pairwise statistics for two breeds -----------------------------
pair_stats <- function(b1, b2) {
  n1 <- af[b1, "X_chromosomes"]; n2 <- af[b2, "X_chromosomes"]
  a1 <- af[b1, "Count_A"];       a2 <- af[b2, "Count_A"]
  p1 <- a1 / n1; p2 <- a2 / n2

  # Nei's Gst (classic FST): (Ht - Hs) / Ht
  Hs <- (2 * p1 * (1 - p1) + 2 * p2 * (1 - p2)) / 2
  pbar <- (p1 + p2) / 2
  Ht <- 2 * pbar * (1 - pbar)
  Gst <- (Ht - Hs) / Ht

  # Hudson's FST (corrects for small sample size; Bhatia et al. 2013)
  num <- (p1 - p2)^2 - p1 * (1 - p1) / (n1 - 1) - p2 * (1 - p2) / (n2 - 1)
  den <- p1 * (1 - p2) + p2 * (1 - p1)
  Hudson <- num / den

  # Nei's (1972) standard genetic distance
  Jx <- p1^2 + (1 - p1)^2; Jy <- p2^2 + (1 - p2)^2
  Jxy <- p1 * p2 + (1 - p1) * (1 - p2)
  NeiD <- -log(Jxy / sqrt(Jx * Jy))

  # Fisher's exact test on allele counts (are allele frequencies different?)
  tab <- matrix(c(n1 - a1, a1, n2 - a2, a2), nrow = 2, byrow = TRUE)
  p <- fisher.test(tab)$p.value

  small <- min(n1, n2) < 10
  data.frame(Comparison = paste(b1, "vs", b2),
             Freq_A_1 = round(p1, 3), Freq_A_2 = round(p2, 3),
             Nei_Gst = round(Gst, 3),
             Hudson_FST = round(max(Hudson, 0), 3),      # negative values set to 0
             Nei_distance = round(NeiD, 3),
             Fisher_p = round(p, 4),
             Significant = ifelse(p < 0.05, "Yes", "No"),
             Note = ifelse(small, "White Fulani n = 2: descriptive only", ""))
}

t412 <- rbind(pair_stats("Muturu", "Kuri"),
              pair_stats("Muturu", "White Fulani"),
              pair_stats("Kuri", "White Fulani"))

cat("\n===== TABLE 4.15: Pairwise genetic differentiation (c.699C>A) =====\n")
print(t412, row.names = FALSE)
write.csv(t412, "R_outputs/Table4.12_pairwise_differentiation.csv", row.names = FALSE)

# ---- 2. Overall differentiation across all three breeds --------------
p <- af$Count_A / af$X_chromosomes
Hs <- mean(2 * p * (1 - p)); pbar <- mean(p); Ht <- 2 * pbar * (1 - pbar)
Gst_all <- (Ht - Hs) / Ht
p_all <- fisher.test(cbind(af$Count_C, af$Count_A))$p.value
cat("\nOverall Gst (3 breeds) =", round(Gst_all, 3),
    "  Fisher's exact test p =", round(p_all, 4), "\n")
cat("Guide (Wright 1978): FST 0-0.05 little, 0.05-0.15 moderate,",
    "0.15-0.25 great, > 0.25 very great differentiation\n")

write.csv(data.frame(Statistic = c("Overall Gst (3 breeds)", "Fisher exact p (3 breeds)"),
                     Value = round(c(Gst_all, p_all), 4)),
          "R_outputs/Table4.13_overall_differentiation.csv", row.names = FALSE)

# ---- 3. Figure: Nei distance matrix as a simple tree -----------------
D <- matrix(0, 3, 3, dimnames = list(af$Breed, af$Breed))
D["Muturu", "Kuri"] <- D["Kuri", "Muturu"] <- t412$Nei_distance[1]
D["Muturu", "White Fulani"] <- D["White Fulani", "Muturu"] <- t412$Nei_distance[2]
D["Kuri", "White Fulani"] <- D["White Fulani", "Kuri"] <- t412$Nei_distance[3]
write.csv(D, "R_outputs/Table4.14_Nei_distance_matrix.csv")

png("R_outputs/Fig_breed_dendrogram.png", width = 1400, height = 1100, res = 200)
plot(hclust(as.dist(D), method = "average"),
     main = "UPGMA dendrogram of breeds (Nei's distance, BMP15 c.699)",
     xlab = "", sub = "", ylab = "Nei's genetic distance")
dev.off()

cat("\nDone. Results saved in R_outputs.\n")
})
while (!is.null(dev.list())) dev.off()

# ======================================================================
# STEP 9: 08_phylogeny
# ======================================================================
local({
# =====================================================================
# BMP15 dissertation - SCRIPT 08: Haplotypes and phylogenetic relationship
# Needs internet (downloads BMP15 sequences of other species from NCBI)
# Needs packages: ape, rentrez, Biostrings, pwalign (install once - see chat)
# =====================================================================

for (p in c("ape", "rentrez", "Biostrings", "pwalign")) {
  if (!requireNamespace(p, quietly = TRUE)) stop("Package '", p, "' is not installed. See the install step in the chat.")
}
library(ape)

dir.create("R_outputs", showWarnings = FALSE)

# ---- 1. Our Nigerian cattle haplotypes (amplicon, in-frame c.607-975) --
aln <- read.dna("data/raw/BMP15_alignment.fasta", format = "fasta", as.character = TRUE)
ref <- toupper(aln[1, ]); ref <- ref[ref != "-"]
cds <- paste(ref, collapse = "")
H1 <- strsplit(substr(cds, 607, 975), "")[[1]]          # 369 bp = 123 codons
H3 <- H1; H3[699 - 606] <- "A"                            # c.699A, no deletion
H2 <- H3; H2[(788 - 606):(791 - 606)] <- "-"              # c.699A + c.788_791del

cat("Haplotype length:", length(H1), "bp\n")              # should be 369
cat("Deleted bases in H2:", paste(H1[182:185], collapse = ""), "\n")   # AAGA

seqs <- list(Nigeria_H1_ref_C699 = H1,
             Nigeria_H2_A699_del = H2,
             Nigeria_H3_A699_putative = H3)

# Haplotype description table
hap <- data.frame(
  Haplotype = c("H1", "H2", "H3"),
  c699 = c("C", "A", "A"),
  c788_791 = c("Present (AAGA)", "Deleted", "Present (AAGA)"),
  Protein = c("Normal (403 aa)", "Truncated (272 aa)", "Normal (403 aa)"),
  Evidence = c("Clean single trace in 12 animals",
               "Deletion template in all 15 carriers",
               "Inferred: 'A only' animals also show the non-deleted template"),
  Breeds_observed = c("Muturu, Kuri, White Fulani", "Muturu, Kuri, White Fulani", "Kuri"))
cat("\n===== TABLE 4.17: BMP15 haplotypes observed =====\n")
print(hap, row.names = FALSE, right = FALSE)
write.csv(hap, "R_outputs/Table4.15_haplotypes.csv", row.names = FALSE)

# ---- 2. Download BMP15 mRNA of related species from NCBI -------------
species <- c("Bos taurus", "Bos indicus", "Bos mutus", "Bison bison",
             "Bubalus bubalis", "Capra hircus", "Ovis aries", "Sus scrofa", "Homo sapiens")
info <- data.frame()
for (sp in species) {
  term <- paste0("BMP15[Gene Name] AND ", sp, "[Organism] AND refseq[filter] AND biomol_mrna[PROP]")
  hit <- tryCatch(rentrez::entrez_search(db = "nuccore", term = term, retmax = 1),
                  error = function(e) NULL)
  if (is.null(hit) || length(hit$ids) == 0) { cat("No sequence found for", sp, "\n"); next }
  txt <- rentrez::entrez_fetch(db = "nuccore", id = hit$ids[1], rettype = "fasta")
  tmp <- tempfile(fileext = ".fa"); writeLines(txt, tmp)
  s <- Biostrings::readDNAStringSet(tmp)[[1]]
  acc <- sub(" .*", "", names(Biostrings::readDNAStringSet(tmp))[1])

  # Align our 369-bp region to the mRNA and keep the matching part
  pa <- pwalign::pairwiseAlignment(pattern = Biostrings::DNAString(paste(H1, collapse = "")),
                                   subject = s, type = "global-local")
  p_chr <- strsplit(as.character(pwalign::alignedPattern(pa)), "")[[1]]
  s_chr <- strsplit(as.character(pwalign::alignedSubject(pa)), "")[[1]]
  region <- s_chr[p_chr != "-"]                  # keep cattle reference coordinates
  if (length(region) != length(H1)) { cat("Alignment problem for", sp, "- skipped\n"); next }

  label <- paste0(gsub(" ", "_", sp), "_", acc)
  seqs[[label]] <- region
  info <- rbind(info, data.frame(Species = sp, Accession = acc,
                                 Identity_to_H1_pct = round(100 * mean(region == H1), 2)))
  cat("Downloaded", sp, acc, "\n")
  Sys.sleep(0.4)                                 # be polite to NCBI
}
cat("\n===== SEQUENCES USED FOR TABLE 4.18: Comparison sequences from GenBank =====\n")
print(info, row.names = FALSE)
write.csv(info, "R_outputs/Table4.16_genbank_sequences.csv", row.names = FALSE)

# ---- 3. Build the alignment and save it (for MEGA / DnaSP if needed) --
x <- as.matrix(as.DNAbin(lapply(seqs, tolower)))
write.dna(x, "R_outputs/08_BMP15_haplotypes_species_alignment.fasta", format = "fasta")

# ---- 4. Neighbour-joining tree with 1000 bootstrap replicates ----------
make_tree <- function(xx) nj(dist.dna(xx, model = "K80", pairwise.deletion = TRUE))
tr <- make_tree(x)
set.seed(123)
bs <- boot.phylo(tr, x, make_tree, B = 1000, quiet = TRUE)
tr$node.label <- bs / 10                         # as percentages
og <- grep("Homo_sapiens", tr$tip.label, value = TRUE)
if (length(og) == 1) tr <- root(tr, outgroup = og, resolve.root = TRUE, edgelabel = TRUE)

png("R_outputs/Fig_BMP15_NJ_tree.png", width = 2000, height = 1400, res = 200)
plot(tr, cex = 0.8, no.margin = FALSE,
     main = "NJ tree of BMP15 exon 2 region (K80, 1000 bootstraps)")
nodelabels(round(as.numeric(tr$node.label)), frame = "none", adj = c(1.2, -0.4), cex = 0.7)
add.scale.bar()
dev.off()
write.tree(tr, "R_outputs/08_BMP15_NJ_tree.nwk")

# ---- 5. Pairwise distances ---------------------------------------------
d <- as.matrix(dist.dna(x, model = "K80", pairwise.deletion = TRUE))
write.csv(round(d, 4), "R_outputs/Table4.17_K80_distance_matrix.csv")

cat("\nDone. Results saved in R_outputs.\n")
})
while (!is.null(dev.list())) dev.off()

# ======================================================================
# STEP 10: 09_selection
# ======================================================================
local({
# =====================================================================
# BMP15 dissertation - SCRIPT 09: Signatures of selection
# Needs packages: ape, pegas (install pegas once - see chat)
# Uses outputs of Script 06 (Table 4.13) and Script 08 (alignment)
# =====================================================================

for (p in c("ape", "pegas")) {
  if (!requireNamespace(p, quietly = TRUE)) stop("Package '", p, "' is not installed. See the install step in the chat.")
}
library(ape)

x <- read.dna("R_outputs/08_BMP15_haplotypes_species_alignment.fasta", format = "fasta")
xc <- toupper(as.character(x))
h1 <- xc["Nigeria_H1_ref_C699", ]

# ---- 1. Corrected comparison table (identity over aligned bases only) ----
others <- rownames(xc)[!grepl("^Nigeria", rownames(xc))]
t416b <- do.call(rbind, lapply(others, function(r) {
  s <- xc[r, ]; ok <- s != "-"
  data.frame(Sequence = r,
             Aligned_bases = sum(ok),
             Coverage_pct = round(100 * mean(ok), 1),
             Identity_pct_aligned = round(100 * mean(s[ok] == h1[ok]), 2))
}))
cat("\n===== TABLE 4.18: Identity to Nigerian H1 over aligned bases =====\n")
print(t416b, row.names = FALSE)
write.csv(t416b, "R_outputs/Table4.16b_identity_aligned.csv", row.names = FALSE)

# ---- 2. Synonymous vs amino-acid-changing differences per species -------
bases <- c("T", "C", "A", "G")
codon_names <- paste0(rep(bases, each = 16), rep(rep(bases, each = 4), 4), rep(bases, 16))
code <- setNames(strsplit("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "")[[1]], codon_names)
codons <- function(s) sapply(seq(1, length(s), 3), function(i) paste(s[i:(i + 2)], collapse = ""))
c1 <- codons(h1)

full <- rownames(xc)[apply(xc, 1, function(r) !any(r == "-"))]   # gap-free sequences only
cmp <- setdiff(full[!grepl("^Nigeria_H1", full)], character(0))
t418 <- do.call(rbind, lapply(cmp, function(r) {
  c2 <- codons(xc[r, ])
  diff <- c1 != c2
  syn <- sum(diff & code[c1] == code[c2])
  non <- sum(diff & code[c1] != code[c2])
  data.frame(Sequence = r, Codons_differing = sum(diff),
             Synonymous = syn, Amino_acid_changing = non)
}))

# dN/dS (Li 1993 method, ape::dnds)
# Calculated one pair at a time (H1 against each sequence), so that one
# problem sequence - for example one identical to H1 - cannot stop the rest.
dn_pair <- function(r) {
  if (t418$Codons_differing[t418$Sequence == r] == 0) return(NA)   # identical to H1: nothing to estimate
  v <- tryCatch(as.matrix(dnds(x[c("Nigeria_H1_ref_C699", r), ], code = 1, codonstart = 1, quiet = TRUE))[1, 2],
                error = function(e) { message("dN/dS could not be calculated for ", r, ": ", conditionMessage(e)); NA })
  round(v, 3)
}
t418$dN_dS <- sapply(t418$Sequence, dn_pair)
t418$Interpretation <- ifelse(t418$Codons_differing == 0, "Identical to H1",
                       ifelse(is.na(t418$dN_dS) | !is.finite(t418$dN_dS), "Not estimable",
                       ifelse(t418$dN_dS < 1, "Purifying selection (dN/dS < 1)",
                       ifelse(t418$dN_dS > 1, "Positive selection (dN/dS > 1)", "Neutral"))))
cat("\n===== TABLE 4.19: Coding differences and dN/dS versus Nigerian H1 =====\n")
print(t418, row.names = FALSE)
write.csv(t418, "R_outputs/Table4.18_dNdS.csv", row.names = FALSE)

# ---- 3. Within-breed diversity and Tajima's D --------------------------
# Each X chromosome is H1 (c.699C) or H3 (c.699A). The 4-bp deletion is
# treated as missing data, as DnaSP does for gaps, so 1 segregating site.
af <- read.csv("R_outputs/Table4.11_allele_freq_sex_aware.csv", stringsAsFactors = FALSE)
hapC <- x["Nigeria_H1_ref_C699", ]
hapA <- x["Nigeria_H3_A699_putative", ]

pop_stats <- function(nC, nA, name) {
  n <- nC + nA
  s <- rbind(hapC[rep(1, nC), ], hapA[rep(1, nA), ])
  rownames(s) <- paste0(name, "_", seq_len(n))
  td <- tryCatch(pegas::tajima.test(s), error = function(e) NULL)
  data.frame(Population = name, X_chromosomes = n,
             Haplotypes = (nC > 0) + (nA > 0),
             Segregating_sites = length(seg.sites(s)),
             Hd = round(pegas::hap.div(s), 3),
             Pi = round(pegas::nuc.div(s), 5),
             Tajima_D = if (is.null(td) || n < 4) NA else round(td$D, 3),
             P_value = if (is.null(td) || n < 4) NA else round(td$Pval.beta, 3))
}
t419 <- do.call(rbind, lapply(seq_len(nrow(af)), function(i)
  pop_stats(af$Count_C[i], af$Count_A[i], af$Breed[i])))
t419$Note <- ifelse(t419$X_chromosomes < 4, "Too few - not tested", "")

cat("\n===== TABLE 4.14: Haplotype diversity, nucleotide diversity, Tajima's D =====\n")
print(t419, row.names = FALSE)
write.csv(t419, "R_outputs/Table4.19_diversity_TajimaD.csv", row.names = FALSE)

cat("\nGuide: Tajima's D near 0 = neutral; D > 0 (significant) = balancing selection",
    "\n       or population decline; D < 0 (significant) = purifying/positive selection",
    "\n       or population expansion. p < 0.05 = significant.\n")
cat("\nDone. Results saved in R_outputs.\n")
})
while (!is.null(dev.list())) dev.off()

# ======================================================================
# STEP 11: 10_LD_AMOVA_tables
# ======================================================================
local({
# =====================================================================
# BMP15 dissertation - SCRIPT 10: Linkage, AMOVA, sex and breed-type tables
# Extra analyses (items 5, 10, 13) - needs ape + pegas (already installed)
# =====================================================================

library(ape)
dat <- read.csv("R_outputs/02_clean_genotype_table.csv", stringsAsFactors = FALSE)
dat$Breed <- factor(dat$Breed, levels = c("Muturu", "Kuri", "White Fulani"))

# ---- 1. Linkage disequilibrium between c.699 and the deletion ----------
# Genotype-based ("composite") LD, because phase cannot be read from
# forward-only Sanger traces. WF_18 (uncertain c.699 call) excluded.
ld <- dat[dat$c699_group != "Uncertain", ]
A_present   <- as.numeric(ld$c699_group %in% c("C and A", "A only"))
Del_present <- as.numeric(ld$del4_status == "Carrier")

tab_ld <- table(c699_A = ifelse(A_present == 1, "A present", "A absent"),
                Deletion = ifelse(Del_present == 1, "Carrier", "Non-carrier"))
cat("\n===== 2x2 table: c.699 A allele vs deletion =====\n"); print(tab_ld)

r  <- cor(A_present, Del_present)
pA <- mean(A_present); pD <- mean(Del_present)
D  <- mean(A_present * Del_present) - pA * pD
Dmax <- if (D > 0) min(pA * (1 - pD), (1 - pA) * pD) else min(pA * pD, (1 - pA) * (1 - pD))
fp <- fisher.test(tab_ld)$p.value

t420 <- data.frame(Animals = length(A_present), r = round(r, 3), r_squared = round(r^2, 3),
                   D_prime = round(D / Dmax, 3), Fisher_p = signif(fp, 3))
cat("\n===== SECTION 4.6.1: Linkage disequilibrium (c.699 vs c.788_791del) =====\n")
print(t420, row.names = FALSE)
write.csv(t420, "R_outputs/Table4.20_LD.csv", row.names = FALSE)

# ---- 2. AMOVA on X-chromosome alleles at c.699 --------------------------
af <- read.csv("R_outputs/Table4.11_allele_freq_sex_aware.csv", stringsAsFactors = FALSE)
af <- af[af$Breed != "All breeds", ]

run_amova <- function(a) {
  allele <- unlist(mapply(function(c, A) c(rep(0, c), rep(1, A)), a$Count_C, a$Count_A))
  pop <- factor(rep(a$Breed, a$Count_C + a$Count_A))
  d <- dist(allele)
  set.seed(123)
  am <- pegas::amova(d ~ pop, nperm = 1000)
  vc <- am$varcomp
  pct <- 100 * vc$sigma2 / sum(vc$sigma2)
  data.frame(Source = c("Among breeds", "Within breeds"),
             df = am$tab$df[1:2],
             SSD = round(am$tab$SSD[1:2], 3),
             Variance_component = round(vc$sigma2, 4),
             Percent_variation = round(pct, 1),
             PhiST = c(round(vc$sigma2[1] / sum(vc$sigma2), 3), NA),
             P_value = c(vc$P.value[1], NA))
}

t421a <- run_amova(af)
t421b <- run_amova(af[af$Breed %in% c("Muturu", "Kuri"), ])
cat("\n===== TABLE 4.16 (a): AMOVA, three breeds =====\n"); print(t421a, row.names = FALSE)
cat("\n===== TABLE 4.16 (b): AMOVA, Muturu vs Kuri only =====\n"); print(t421b, row.names = FALSE)
write.csv(rbind(cbind(Analysis = "Three breeds", t421a),
                cbind(Analysis = "Muturu vs Kuri", t421b)),
          "R_outputs/Table4.21_AMOVA.csv", row.names = FALSE)

# ---- 3. Deletion carriers by sex and by breed type --------------------
cat("\n===== TABLE 4.6 (a): Deletion status by sex =====\n")
t422 <- table(dat$Sex, dat$del4_status); print(addmargins(t422))
cat("Fisher's exact test (sex vs deletion): p =", signif(fisher.test(t422)$p.value, 3), "\n")
write.csv(as.data.frame.matrix(addmargins(t422)), "R_outputs/Table4.22_deletion_by_sex.csv")

cat("\n===== TABLE 4.6 (b): Deletion status by breed type =====\n")
t423 <- table(dat$Type, dat$del4_status); print(addmargins(t423))
write.csv(as.data.frame.matrix(addmargins(t423)), "R_outputs/Table4.23_deletion_by_type.csv")

cat("\nDone. Results saved in R_outputs.\n")
})
while (!is.null(dev.list())) dev.off()

# ======================================================================
# STEP 12: 11_export_vcf
# ======================================================================
local({
# ======================================================================
# 11_export_vcf.R  –  VCF of the two BMP15 variants for submission to
#                     the European Variation Archive (EVA)
#
# Coordinates: ARS-UCD2.0 (GCA_002263795.4), chromosome X.
# The deletion is written left-aligned, as VCF requires; c.788_791delAAGA
# (HGVS, 3'-aligned) = X:88810169 TGAAA>T (same record as EVA PRJEB89089).
# BMP15 is X-linked: females are diploid, males haploid. Males in which
# both alleles were seen in the trace are written 0/1 (see Chapter 4).
# ======================================================================
# (settings already loaded above)
dat <- read.csv("R_outputs/02_clean_genotype_table.csv", stringsAsFactors = FALSE)
ord <- order(match(substr(dat$Sample, 1, 2), c("MT", "KR", "WF")), dat$ID)
dat <- dat[ord, ]
ids <- sub("^([A-Z]+_[0-9]+)_.*$", "\\1", dat$Sample)
male <- dat$Sex != "Female"

gt699 <- ifelse(dat$c699_group == "C only",  ifelse(male, "0", "0/0"),
         ifelse(dat$c699_group == "A only",  ifelse(male, "1", "1/1"),
         ifelse(dat$c699_group == "C and A", "0/1", ifelse(male, ".", "./."))))
gtdel <- ifelse(dat$del4_status == "Carrier", "0/1", ifelse(male, "0", "0/0"))

pos <- function(c) cfg$genome_offset + (c - cfg$amplicon_start)
header <- c(
  "##fileformat=VCFv4.2",
  paste0("##fileDate=", format(Sys.Date(), "%Y%m%d")),
  "##source=BMP15_pipeline (Sanger forward reads; peak-height review of chromatograms)",
  paste0("##reference=", cfg$assembly),
  paste0("##contig=<ID=", cfg$genome_chrom, ",assembly=", cfg$assembly, ">"),
  "##INFO=<ID=GENE,Number=1,Type=String,Description=\"Gene symbol\">",
  "##INFO=<ID=HGVS_C,Number=1,Type=String,Description=\"Coding change on XM_027535004.1 (HGVS)\">",
  "##INFO=<ID=EFFECT,Number=1,Type=String,Description=\"Predicted effect on the protein\">",
  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype. BMP15 is X-linked: females diploid, males haploid; males showing both alleles in the trace are given as 0/1\">",
  paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", ids), collapse = "\t"))
row699 <- c(cfg$genome_chrom, pos(699), ".", "C", "A", ".", "PASS",
            "GENE=BMP15;HGVS_C=c.699C>A;EFFECT=synonymous_variant_p.Gly233Gly", "GT", gt699)
rowdel <- c(cfg$genome_chrom, pos(785), ".", "TGAAA", "T", ".", "PASS",
            "GENE=BMP15;HGVS_C=c.788_791delAAGA;EFFECT=frameshift_variant_p.Lys263SerfsTer11", "GT", gtdel)
ref <- strsplit(paste(readLines(cfg$reference_cds)[-1], collapse = ""), "")[[1]]
stopifnot(ref[699] == "C", paste(ref[785:789], collapse = "") == "TGAAA")   # REF alleles must match
writeLines(c(header, paste(row699, collapse = "\t"), paste(rowdel, collapse = "\t")),
           "R_outputs/BMP15_variants_ARS-UCD2.0.vcf")
cat("VCF written for", length(ids), "animals: R_outputs/BMP15_variants_ARS-UCD2.0.vcf\n")
cat("c.699  ->", cfg$genome_chrom, pos(699), " C>A\n")
cat("deletion ->", cfg$genome_chrom, pos(785), " TGAAA>T\n")
})
while (!is.null(dev.list())) dev.off()

# ======================================================================
# STEP 13: 12_final_figures
# ======================================================================
local({
# ======================================================================
# 12_final_figures.R  -  final dissertation versions of Figure 4.7 and
#                        Figure 4.9, drawn from the pipeline outputs
#
# Inputs (made by earlier steps):
#   R_outputs/Table4.14_Nei_distance_matrix.csv   (step 07)
#   R_outputs/08_BMP15_NJ_tree.nwk                (step 08)
# Outputs:
#   R_outputs/Fig4_7_dendrogram.png
#   R_outputs/Fig4_9_NJ_tree_cladogram.png
#
# Run from the repository folder (open BMP15_pipeline.Rproj):  source("BMP15_dissertation_analysis.R")
# Package: ape
# ======================================================================
suppressPackageStartupMessages(library(ape))
dir.create("R_outputs", showWarnings = FALSE)

# ---- Figure 4.7: UPGMA dendrogram of the breeds (Nei's D) -------------
nei <- as.matrix(read.csv("R_outputs/Table4.14_Nei_distance_matrix.csv",
                          row.names = 1, check.names = FALSE))
cat("Nei's D matrix:\n"); print(nei)

upgma <- hclust(as.dist(nei), method = "average")   # UPGMA
upgma$height <- upgma$height / 2                     # branch height = D / 2 (ultrametric tree)

png("R_outputs/Fig4_7_dendrogram.png", width = 1900, height = 870, res = 300, family = "serif")
par(mar = c(4, 1, 1, 9))
dend <- as.dendrogram(upgma)
plot(dend, horiz = TRUE, axes = FALSE, edgePar = list(lwd = 2, col = "grey15"),
     nodePar = list(pch = NA), leaflab = "none",
     xlim = c(max(upgma$height) * 1.05, 0))
axis(1, at = seq(0, 0.125, by = 0.025), cex.axis = 0.8)
mtext("Genetic distance (Nei's D / 2)", side = 1, line = 2.5, cex = 0.9)
lab <- labels(dend)
text(x = 0, y = seq_along(lab), labels = lab, pos = 4, xpd = TRUE, cex = 1)
dev.off()

# ---- Figure 4.9: neighbour-joining tree as a cladogram ---------------
tr <- read.tree("R_outputs/08_BMP15_NJ_tree.nwk")

# The root has an extra node of zero length (Sus + ruminants); collapse it so
# that human, pig and the ruminants all join at the root, as in the NJ result.
extra <- which(tr$node.label == "")
if (length(extra) == 1) {
  node <- Ntip(tr) + extra
  tmp <- tr
  tmp$edge.length <- rep(1, nrow(tmp$edge))
  tmp$edge.length[tmp$edge[, 2] == node] <- 0
  tr <- di2multi(tmp, tol = 1e-8)
}
tr <- ladderize(tr, right = FALSE)

# Tip labels: italic species name + accession; Nigerian haplotypes in red
nice <- c(
  Homo_sapiens_NM_005448.2       = "Homo sapiens|NM_005448.2",
  Sus_scrofa_NM_001005155.2      = "Sus scrofa|NM_001005155.2",
  Ovis_aries_NM_001114767.2      = "Ovis aries|NM_001114767.2",
  Capra_hircus_NM_001285588.1    = "Capra hircus|NM_001285588.1",
  Bubalus_bubalis_XM_006059485.4 = "Bubalus bubalis|XM_006059485.4",
  Bos_mutus_XM_005898488.2       = "Bos mutus|XM_005898488.2",
  Bos_indicus_XM_070783943.1     = "Bos indicus|XM_070783943.1",
  Bison_bison_XM_010858078.1     = "Bison bison|XM_010858078.1",
  Bos_taurus_NM_001031752.1      = "Bos taurus|NM_001031752.1",
  Nigeria_H1_ref_C699            = "Nigeria H1 (reference, C at c.699)",
  Nigeria_H2_A699_del            = "Nigeria H2 (A at c.699 + deletion)",
  Nigeria_H3_A699_putative       = "Nigeria H3 (A at c.699, putative)")

png("R_outputs/Fig4_9_NJ_tree_cladogram.png", width = 2250, height = 1560, res = 300, family = "serif")
par(mar = c(1, 1, 1, 1))
plot(tr, use.edge.length = FALSE, show.tip.label = FALSE, edge.width = 2,
     edge.color = "grey15", x.lim = c(0, max(node.depth(tr)) * 2.1), no.margin = TRUE)
pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
for (i in seq_len(Ntip(tr))) {
  x <- pp$xx[i] + 0.15; y <- pp$yy[i]
  lb <- nice[tr$tip.label[i]]
  if (grepl("^Nigeria", lb)) {
    text(x, y, lb, adj = 0, cex = 0.85, font = 2, col = "#8B1A1A")
  } else {
    sp <- sub("\\|.*", "", lb); acc <- sub(".*\\|", "", lb)
    text(x, y, bquote(italic(.(sp)) ~ "(" * .(acc) * ")"), adj = 0, cex = 0.85)
  }
}
boot <- suppressWarnings(round(as.numeric(tr$node.label)))   # R rounds 46.5 to 46
boot_lab <- ifelse(is.na(boot), "", as.character(boot))
nodelabels(boot_lab, frame = "none", adj = c(1.2, -0.5), cex = 0.75, col = "#1f4e79")
dev.off()

cat("Saved: R_outputs/Fig4_7_dendrogram.png and R_outputs/Fig4_9_NJ_tree_cladogram.png\n")
cat("Figure 4.9 is a cladogram: branch lengths are not drawn to scale; bootstrap values (%) are at the nodes.\n")
})
while (!is.null(dev.list())) dev.off()

# ======================================================================
# STEP 14: 13_extra_figures
# ======================================================================
local({
# ======================================================================
# 13_extra_figures.R  -  figures drawn directly from the raw data:
#   1. chromatograms at c.699 (MT_58, MT_74, KR_132)
#   2. chromatograms around the c.788_791delAAGA deletion (MT_58, MT_65, KR_132)
#   3. schematic of the wild-type and truncated BMP15 proteins
#   4. haplotype network (H1, H3, H2) with the breeds of the animals
#      carrying each haplotype
# Uses base R only. Run after 00_trace_qc.R (needs the curated table).
# ======================================================================
# (settings already loaded above)
dir.create("R_outputs", showWarnings = FALSE)

# ---- ABIF reader (same as in 00_trace_qc.R) --------------------------
read_abif <- function(path) {
  b <- readBin(path, "raw", n = file.info(path)$size)
  if (rawToChar(b[1:4]) != "ABIF") stop(path, " is not an ABIF file")
  i32 <- function(off) readBin(b[(off + 1):(off + 4)], "integer", size = 4, endian = "big")
  i16 <- function(off) readBin(b[(off + 1):(off + 2)], "integer", size = 2, endian = "big", signed = TRUE)
  n_entries <- i32(6 + 12); dir_off <- i32(6 + 20); tags <- list()
  for (k in 0:(n_entries - 1)) {
    e <- dir_off + k * 28
    name <- rawToChar(b[(e + 1):(e + 4)]); num <- i32(e + 4); type <- i16(e + 8)
    cnt <- i32(e + 12); size <- i32(e + 16)
    if (size <= 0) next
    off <- if (size <= 4) e + 20 else i32(e + 20)
    tags[[paste0(name, ".", num)]] <- list(type = type, cnt = cnt, raw = b[(off + 1):(off + size)])
  }
  get_short <- function(key) readBin(tags[[key]]$raw, "integer", size = 2, n = tags[[key]]$cnt,
                                     endian = "big", signed = TRUE)
  order <- strsplit(rawToChar(tags[["FWO_.1"]]$raw[1:4]), "")[[1]]
  trace <- lapply(9:12, function(n) get_short(paste0("DATA.", n))); names(trace) <- order
  list(seq = rawToChar(tags[["PBAS.2"]]$raw), ploc = get_short("PLOC.2") + 1L, trace = trace)
}

ref_lines <- readLines(cfg$reference_cds)
ref <- strsplit(toupper(paste(ref_lines[!grepl("^>", ref_lines)], collapse = "")), "")[[1]]
COL <- c(A = "#1a9850", C = "#2166ac", G = "#222222", T = "#d73027")
MARK <- "#d6604d"

# one chromatogram panel between coding positions c0 and c1
chrom_panel <- function(sample, c0, c1, title, mark_pos = NULL, mark_span = NULL, mark_lab = "") {
  r <- read_abif(file.path(cfg$ab1_dir, paste0(sample, ".ab1")))
  i <- regexpr(cfg$anchor_c787, r$seq, fixed = TRUE)[1]
  if (i < 0) stop("Anchor sequence not found in ", sample)
  p787 <- i + nchar(cfg$anchor_c787) - 1           # base number of c.787 in this read
  idx <- function(cc) p787 + (cc - 787)
  a <- r$ploc[idx(c0)] - 6; b <- r$ploc[idx(c1)] + 6; xs <- a:b
  ymax <- max(sapply(r$trace, function(t) max(t[xs])))
  plot(NA, xlim = c(a, b), ylim = c(0, ymax * 1.3), axes = FALSE, xlab = "", ylab = "", xaxs = "i")
  if (!is.null(mark_span))
    rect(r$ploc[idx(mark_span[1])] - 5, 0, r$ploc[idx(mark_span[2])] + 5, ymax * 1.18,
         col = adjustcolor(MARK, 0.15), border = NA)
  for (bse in c("A", "C", "G", "T")) lines(xs, r$trace[[bse]][xs], col = COL[bse], lwd = 1.1)
  if (!is.null(mark_pos)) abline(v = r$ploc[idx(mark_pos)], col = MARK, lty = 2)
  lab_x <- r$ploc[idx(if (!is.null(mark_pos)) mark_pos else mark_span[1])]
  text(lab_x, ymax * 1.24, mark_lab, col = MARK, adj = 0, cex = 0.8)
  axis(1, at = r$ploc[idx(c0:c1)], labels = ref[c0:c1], tick = FALSE, line = -0.8,
       cex.axis = 0.75, family = "mono")
  mtext(title, side = 3, adj = 0, line = 0.2, cex = 0.75)
}

chrom_figure <- function(file, samples, titles, c0, c1, main, ...) {
  png(file, width = 2400, height = 1700, res = 250)
  par(mfrow = c(3, 1), mar = c(2.2, 1, 1.6, 1), oma = c(2, 0, 3, 0))
  for (k in seq_along(samples)) chrom_panel(samples[k], c0, c1, titles[k], ...)
  mtext(sprintf("Reference sequence (c.%d-%d)", c0, c1), side = 1, outer = TRUE, line = 0.5, cex = 0.8)
  mtext(main, side = 3, outer = TRUE, adj = 0, line = 1.2, cex = 0.9)
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
  legend("topright", legend = names(COL), col = COL, lwd = 2, horiz = TRUE, bty = "n", cex = 0.8)
  dev.off()
}

# ---- 1-2. Chromatograms ----------------------------------------------
chrom_figure("R_outputs/Fig_chromatogram_c699.png",
             c("MT_58_F_T", "MT_74_F_T", "KR_132_M_T"),
             c("MT_58 (female): C only", "MT_74 (female): C and A, double peak", "KR_132 (male): A only"),
             690, 708, "BMP15 c.699C>A: forward-strand chromatograms",
             mark_pos = 699, mark_lab = "c.699")
chrom_figure("R_outputs/Fig_chromatogram_deletion.png",
             c("MT_58_F_T", "MT_65_F_T", "KR_132_M_T"),
             c("MT_58 (non-carrier): single clean trace after c.787",
               "MT_65 (female carrier): trace splits into two after c.787",
               "KR_132 (male carrier): trace splits into two after c.787"),
             772, 812, "BMP15 c.788_791delAAGA: overlapping traces downstream of the deletion",
             mark_span = c(788, 791), mark_lab = "c.788_791 (AAGA)")

# ---- 3. Protein schematic -------------------------------------------
# Domain limits of the 403-aa precursor (from 03_protein_effect.R and the
# literature): signal peptide 1-18 (predicted), pro-region 19-274, furin
# site RRAR 275-278, mature domain 279-403. Deletion protein: 272 aa,
# residues 263-272 altered by the frameshift.
png("R_outputs/Fig_protein_schematic.png", width = 2250, height = 720, res = 300)
par(mar = c(2, 7, 1, 1))
plot(NA, xlim = c(1, 420), ylim = c(0, 3.1), axes = FALSE, xlab = "", ylab = "")
bar <- function(y, a, b, col) rect(a, y, b, y + 0.5, col = col, border = NA)
bar(2, 1, 18, "#9e9e9e"); bar(2, 18, 275, "#cfcfcf"); bar(2, 275, 279, "#d6604d"); bar(2, 279, 403, "#1b7837")
bar(0.6, 1, 18, "#9e9e9e"); bar(0.6, 18, 263, "#cfcfcf"); bar(0.6, 263, 272, "#e08214")
mtext("Wild type\n(403 aa)", side = 2, at = 2.25, las = 1, line = 0.5, cex = 0.8)
mtext("c.788_791del\n(272 aa)", side = 2, at = 0.85, las = 1, line = 0.5, cex = 0.8)
text(9, 2.75, "SP*", cex = 0.7, col = "#444444")
text(146, 2.75, "Pro-region", cex = 0.8, col = "#444444")
text(277, 1.7, "RRAR (aa 275-278)", cex = 0.7, col = "#b2182b")
text(341, 2.75, "Mature domain (aa 279-403)", cex = 0.8, col = "#1b7837")
text(276, 0.85, "altered aa 263-272, stop", adj = 0, cex = 0.75, col = "#444444")
axis(1, at = c(1, 100, 200, 300, 403), cex.axis = 0.7, col = "#555555", line = -0.5)
mtext("aa", side = 1, at = 418, line = 0.3, cex = 0.7)
dev.off()

# ---- 4. Haplotype network -------------------------------------------
# Animals carrying each haplotype (genotyped animals only):
#   H1 (c.699C, no deletion): c.699 call C or C/A
#   H3 (c.699A, no deletion): c.699 call A only
#   H2 (c.699A + deletion):   deletion carriers
rv <- read.csv(cfg$review_table, stringsAsFactors = FALSE)
rv <- rv[rv$Use_in_analysis == "Yes", ]
breeds <- c("Muturu", "Kuri", "White Fulani")
bcol <- c(Muturu = "#1b7837", Kuri = "#e08214", "White Fulani" = "#2166ac")
cnt <- function(keep) table(factor(rv$Breed[keep], levels = breeds))
hap <- list(
  H1 = list(x = 0, lab = "H1\nc.699C, no deletion", n = cnt(grepl("^C", rv$c699_call))),
  H3 = list(x = 4, lab = "H3 (putative)\nc.699A, no deletion", n = cnt(grepl("^A", rv$c699_call))),
  H2 = list(x = 8, lab = "H2\nc.699A + c.788_791del", n = cnt(rv$del4_status == "Carrier")))
cat("Animals carrying each haplotype:\n"); print(sapply(hap, function(h) h$n))

pie_at <- function(x, y, r, n, col) {
  start <- pi / 2; tot <- sum(n)
  for (k in seq_along(n)) {
    if (n[k] == 0) next
    th <- seq(start, start + 2 * pi * n[k] / tot, length.out = 100)
    polygon(c(x, x + r * cos(th)), c(y, y + r * sin(th)), col = col[k], border = "white", lwd = 1.5)
    start <- start + 2 * pi * n[k] / tot
  }
}
png("R_outputs/Fig_haplotype_network.png", width = 2000, height = 840, res = 200)
par(mar = c(0, 0, 2.5, 0))
plot(NA, xlim = c(-1.6, 9.6), ylim = c(-2.7, 2.5), asp = 1, axes = FALSE, xlab = "", ylab = "")
rad <- sapply(hap, function(h) 0.45 + 0.12 * sqrt(sum(h$n)))
links <- list(c("H1", "H3", "1 step\nc.699C>A"), c("H3", "H2", "1 step\n4-bp deletion"))
for (l in links) {
  xa <- hap[[l[1]]]$x; xb <- hap[[l[2]]]$x
  segments(xa + rad[l[1]], 0, xb - rad[l[2]], 0, col = "#555555", lwd = 2)
  segments((xa + xb) / 2, -0.15, (xa + xb) / 2, 0.15, col = "#555555", lwd = 2)
  text((xa + xb) / 2, 0.6, l[3], cex = 0.8, col = "#333333")
}
for (h in names(hap)) {
  pie_at(hap[[h]]$x, 0, rad[h], hap[[h]]$n, bcol)
  text(hap[[h]]$x, -rad[h] - 0.25, paste0(hap[[h]]$lab, "\n(n = ", sum(hap[[h]]$n), " animals)"),
       adj = c(0.5, 1), cex = 0.85)
}
legend("top", legend = breeds, pt.bg = bcol, pch = 21, col = bcol, pt.cex = 2, horiz = TRUE, bty = "n", cex = 0.85)
title("Haplotype network of BMP15 exon 2 in Nigerian cattle (circle size = number of animals carrying each haplotype; colours = breed)",
      adj = 0, cex.main = 0.8, font.main = 1)
dev.off()

cat("Saved: Fig_chromatogram_c699.png, Fig_chromatogram_deletion.png, Fig_protein_schematic.png and Fig_haplotype_network.png in R_outputs/\n")
})
while (!is.null(dev.list())) dev.off()

# ======================================================================
# STEP 15: 99_check_results
# ======================================================================
local({
# ======================================================================
# 99_check_results.R  –  do the outputs match the dissertation?
# Compares key numbers in R_outputs/ with the values reported in
# Chapter 4. Every line should say PASS.
# ======================================================================
# (settings already loaded above)
tol <- cfg$tolerance
res <- data.frame(Check = character(), Expected = numeric(), Obtained = numeric(), Result = character())
chk <- function(label, expected, obtained, tolerance = tol) {
  ok <- !is.na(obtained) && abs(as.numeric(obtained) - expected) <= tolerance
  res[nrow(res) + 1, ] <<- list(label, expected, as.numeric(obtained), ifelse(ok, "PASS", "CHECK"))
}
rd <- function(f) read.csv(file.path("R_outputs", f), stringsAsFactors = FALSE)
val <- function(df, keycol, key, col) df[df[[keycol]] == key, col][1]

t <- rd("Table4.4_deletion_carrier_frequency.csv")
chk("Deletion carriers, Muturu (of 19)", 8, val(t, "Breed", "Muturu", "Carriers"), 0)
chk("Deletion carriers, Kuri (of 8)", 6, val(t, "Breed", "Kuri", "Carriers"), 0)
chk("Deletion carriers, White Fulani (of 2)", 1, val(t, "Breed", "White Fulani", "Carriers"), 0)
chk("Carrier frequency, all breeds", 0.517, val(t, "Breed", "All breeds", "Carrier_freq"))

t <- rd("Table4.10_diversity_females.csv")
chk("Muturu females Ho", 0.538, val(t, "Breed", "Muturu", "Ho"))
chk("Muturu females He", 0.393, val(t, "Breed", "Muturu", "He"))
chk("Muturu PIC", 0.316, val(t, "Breed", "Muturu", "PIC"))
chk("Muturu HWE P", 0.499, val(t, "Breed", "Muturu", "HWE_p"))
chk("Kuri females Ho", 0.400, val(t, "Breed", "Kuri", "Ho"))
chk("Kuri females He", 0.480, val(t, "Breed", "Kuri", "He"))
chk("Kuri PIC", 0.365, val(t, "Breed", "Kuri", "PIC"))

t <- rd("Table4.11_allele_freq_sex_aware.csv")
chk("Freq. A (X chromosomes), Muturu", 0.226, val(t, "Breed", "Muturu", "Freq_A"))
chk("Freq. A (X chromosomes), Kuri", 0.667, val(t, "Breed", "Kuri", "Freq_A"))
chk("X chromosomes counted, all breeds", 45, val(t, "Breed", "All breeds", "X_chromosomes"), 0)

t <- rd("Table4.12_pairwise_differentiation.csv")
chk("Muturu vs Kuri Nei's Gst", 0.197, val(t, "Comparison", "Muturu vs Kuri", "Nei_Gst"))
chk("Muturu vs Kuri Hudson's FST", 0.285, val(t, "Comparison", "Muturu vs Kuri", "Hudson_FST"))
chk("Muturu vs Kuri Nei's D", 0.386, val(t, "Comparison", "Muturu vs Kuri", "Nei_distance"))
chk("Muturu vs Kuri Fisher P", 0.0117, val(t, "Comparison", "Muturu vs Kuri", "Fisher_p"), 0.001)

t <- rd("Table4.20_LD.csv")
chk("LD r-squared", 0.865, t$r_squared[1])
chk("LD D'", 1.000, t$D_prime[1])

amova <- list.files("R_outputs", pattern = "^Table4.21", full.names = TRUE)
if (length(amova)) {
  t <- read.csv(amova[1], stringsAsFactors = FALSE)
  chk("AMOVA PhiST, three breeds", 0.240, t$PhiST[t$Analysis == "Three breeds" & t$Source == "Among breeds"])
  chk("AMOVA PhiST, Muturu vs Kuri", 0.303, t$PhiST[t$Analysis == "Muturu vs Kuri" & t$Source == "Among breeds"])
}

t <- rd("Table4.19_diversity_TajimaD.csv")
chk("Tajima's D, Muturu", 0.684, val(t, "Population", "Muturu", "Tajima_D"))
chk("Tajima's D, Kuri", 1.066, val(t, "Population", "Kuri", "Tajima_D"))

cat("\n===== CHECK AGAINST THE DISSERTATION (Chapter 4) =====\n")
print(res, row.names = FALSE)
write.csv(res, "R_outputs/99_results_check.csv", row.names = FALSE)
cat("\n", sum(res$Result == "PASS"), "of", nrow(res), "checks passed.\n")
if (any(res$Result != "PASS")) cat("Lines marked CHECK differ from the dissertation - look at the matching table before changing anything.\n")
})
while (!is.null(dev.list())) dev.off()

# ======================================================================
# RECORD SOFTWARE VERSIONS AND INPUT CHECKSUMS
# ======================================================================
writeLines(capture.output(sessionInfo()), file.path("logs", "sessionInfo.txt"))
inp <- list.files("data", recursive = TRUE, full.names = TRUE)
write.csv(data.frame(file = inp, md5 = unname(tools::md5sum(inp))), file.path("logs", "input_checksums.csv"), row.names = FALSE)
cat("\nALL STEPS COMPLETED. Results are in R_outputs/.\n")
