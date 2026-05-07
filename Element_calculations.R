#!/usr/bin/env Rscript

Sys.setenv(TZ = 'America/New_York')

library(optparse)
library(readr)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(scales))
library(ggplot2)
library(gtools)
library(ggtext)
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(writexl))


option_list = list(
  make_option(c("-b", "--busco"), type="character", default=NULL,
              help="BUSCO full_table.tsv file", metavar="file.tsv"),
  make_option(c("-n", "--nigon"), type="character", default="gene2Nigon_busco20200927.tsv.gz",
              help="BUSCO id assignment to Nigons [default=%default]", metavar="file.tsv"),
  make_option(c("-a", "--Muller"), type="character", default=NULL,
              help="BUSCO id assignment to mullers", metavar="file.tsv"),
  make_option(c("-c", "--Stevens"), type="character", default=NULL,
              help="BUSCO id assignment to Stevens", metavar="file.tsv"),
  make_option(c("-f", "--fasta"), type="character", default=NULL,
              help="Genome FASTA file to calculate GC content [default=%default]", metavar="file.fa"),
  make_option(c("-w", "--windowSize"), type="integer", default=5e5,
              help="Window size (bp) to bin the BUSCO genes [default=%default]", metavar="integer"),
  make_option(c("-m", "--minimumGenesPerSequence"), type="integer", default=15,
              help="Sequences with fewer than this number of BUSCO genes will not be shown [default=%default]", metavar="integer"),
  make_option(c("-o", "--outPlot"), type="character", default="Nigons.jpeg",
              help="Output image file name for the facet barplot [default=%default]", metavar="file"),
  make_option(c("--stackedOutPlot"), type="character", default="stacked_barplot.png",
              help="Output image file name for the stacked barplot [default=%default]", metavar="file"),
  make_option(c("--tableOutput"), type="character", default="Stevens_breakdown.tsv",
              help="Output TSV file with Stevens element breakdown [default=%default]", metavar="file.tsv"),
  make_option(c("--height"), type="integer", default=6,
              help="Height of plot (in inches) [default=%default]", metavar="integer"),
  make_option(c("--width"), type="integer", default=5,
              help="Width of plot (in inches) [default=%default]", metavar="integer"),
  make_option(c("-s", "--species"), type="character", default="",
              help="Species name (for plot title, italicized) [default=%default]", metavar="Genus_species"),
  make_option(c("--threshold"), type="double", default=NULL,
              help="Threshold (between 0 and 1) for including multiple Stevens elements in label. If not set, default is dominant Stevens element only.", metavar="float"),
  make_option(c("--stacked"), action="store_true", default=FALSE,
              help="Generate a sorted stacked bar plot of physical lengths for Stevens elements [default=%default]")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

stevensDict <- read_tsv(opt$Stevens, col_types = c(col_character(), col_character()))
busco <- suppressWarnings(read_tsv(opt$busco,
                 col_names = c("Busco_id", "Status", "Sequence",
                               "start", "end", "strand", "Score", "Length",
                               "OrthoDB_url", "Description"),
                 col_types = c("ccciicdicc"),
                 comment = "#"))

windwSize <- opt$windowSize
minimumGenesPerSequence <- opt$minimumGenesPerSequence
spName <- opt$species
if(grepl("\\.", spName)) {
  spName <- paste0("*", sub("_", " ", spName), "*")
}

fbusco <- busco %>%
  filter(!Status %in% c("Missing")) %>%
  left_join(stevensDict, by = c("Busco_id" = "Orthogroups")) %>%
  mutate(Stevens = ifelse(is.na(Stevens), "-", Stevens),
         stPos = start) %>%
  filter(Stevens != "-")

stevensSummary <- fbusco %>%
  group_by(Sequence, Stevens) %>%
  filter(!is.na(stPos)) %>%
  mutate(ints = if (max(stPos, na.rm = TRUE) > windwSize) {
           as.numeric(as.character(cut(stPos,
                                       breaks = seq(0, max(stPos, na.rm = TRUE), windwSize),
                                       labels = seq(windwSize, max(stPos, na.rm = TRUE), windwSize))))
         } else {
           NA_real_
         }) %>%
  filter(!is.na(ints)) %>%
  count(ints, Stevens) %>%
  ungroup()

stevensAggregated <- stevensSummary %>%
  group_by(Sequence, Stevens) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  group_by(Sequence) %>%
  mutate(total_genes = sum(n),
         proportion = n / total_genes) %>%
  ungroup()

# Coordinate statistics
coordinateStats <- fbusco %>%
  filter(Stevens != "-") %>%
  group_by(Sequence, Stevens) %>%
  summarise(
    n_genes = n(),
    seq_length = max(end, na.rm = TRUE),
    start_mean = mean(start, na.rm = TRUE),
    end_mean = mean(end, na.rm = TRUE),
    start_sd = sd(start, na.rm = TRUE),
    end_sd = sd(end, na.rm = TRUE),
    start_min = min(start, na.rm = TRUE),
    start_max = max(start, na.rm = TRUE),
    end_min = min(end, na.rm = TRUE),
    end_max = max(end, na.rm = TRUE),
    .groups = "drop"
  )

stevensAggregatedStats <- stevensAggregated %>%
  left_join(coordinateStats, by = c("Sequence", "Stevens"))

# GC content (optional)
get_gc_content <- function(dna_string) {
  g <- letterFrequency(dna_string, "G", as.prob = FALSE)
  c <- letterFrequency(dna_string, "C", as.prob = FALSE)
  n <- letterFrequency(dna_string, "N", as.prob = FALSE)
  total <- length(dna_string) - n
  gc <- g + c
  return(ifelse(total > 0, gc / total, NA))
}

gc_data <- NULL
if (!is.null(opt$fasta)) {
  fasta <- readDNAStringSet(opt$fasta)
  gc_data <- data.frame(
    Sequence = names(fasta),
    GC_content = sapply(fasta, get_gc_content)
  )

  stevensAggregatedStats <- stevensAggregatedStats %>%
    left_join(gc_data, by = "Sequence")

  # Genome-wide summary stats
  seq_lengths <- width(fasta)
  gc_vals <- sapply(fasta, get_gc_content)
  valid_gc <- !is.na(gc_vals)

  genome_summary <- tibble(
    total_length = sum(seq_lengths),
    num_sequences = length(fasta),
    mean_seq_length = mean(seq_lengths),
    median_seq_length = median(seq_lengths),
    min_seq_length = min(seq_lengths),
    max_seq_length = max(seq_lengths),
    sd_seq_length = sd(seq_lengths),
    genome_gc_content = sum(letterFrequency(fasta[valid_gc], c("G", "C"))[,1] +
                             letterFrequency(fasta[valid_gc], c("G", "C"))[,2]) /
                        sum(seq_lengths[valid_gc] - letterFrequency(fasta[valid_gc], "N"))
  )
}

# === Export all summaries to a single Excel workbook ===
summary_list <- list(
  Genome_Stats = if (exists("genome_summary")) genome_summary else tibble(Note = "No genome FASTA provided"),
  Stevens_Aggregated = stevensAggregated,
  Stevens_Detailed = stevensAggregatedStats
)

write_xlsx(summary_list, path = "stevens_summary.xlsx")

