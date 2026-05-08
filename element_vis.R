#!/usr/bin/env Rscript

Sys.setenv(TZ = "America/New_York")

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gtools))
suppressPackageStartupMessages(library(ggtext))

option_list <- list(
  make_option(c("-b", "--busco"), type = "character", default = NULL,
              help = "BUSCO full_table.tsv file", metavar = "file.tsv"),
  make_option(c("-n", "--nigon"), type = "character", default = "gene2Nigon_busco20200927.tsv.gz",
              help = "BUSCO id assignment to Nigons [kept for compatibility; currently unused]", metavar = "file.tsv"),
  make_option(c("-a", "--Muller"), type = "character", default = NULL,
              help = "BUSCO id assignment to Muller elements", metavar = "file.tsv"),
  make_option(c("-c", "--Stevens"), type = "character", default = NULL,
              help = "BUSCO id assignment to Stevens elements", metavar = "file.tsv"),
  make_option(c("-f", "--fasta"), type = "character", default = NULL,
              help = "Optional genome FASTA file for sequence lengths and GC content", metavar = "file.fa"),
  make_option(c("-w", "--windowSize"), type = "integer", default = 5e5,
              help = "Window size (bp) to bin BUSCO genes [default=%default]", metavar = "integer"),
  make_option(c("-m", "--minimumGenesPerSequence"), type = "integer", default = 15,
              help = "Sequences with fewer than this many mapped BUSCO genes are not plotted [default=%default]", metavar = "integer"),
  make_option(c("-o", "--outPlot"), type = "character", default = "Elements.jpeg",
              help = "Output image file name for the windowed/facet barplot [default=%default]", metavar = "file"),
  make_option(c("--stackedOutPlot"), type = "character", default = "stacked_barplot.png",
              help = "Output image file name for the stacked chromosome-length barplot [default=%default]", metavar = "file"),
  make_option(c("--tableOutput"), type = "character", default = "element_breakdown.tsv",
              help = "Output TSV file with detailed element breakdown [default=%default]", metavar = "file.tsv"),
  make_option(c("--height"), type = "integer", default = 6,
              help = "Height of plot in inches [default=%default]", metavar = "integer"),
  make_option(c("--width"), type = "integer", default = 5,
              help = "Width of plot in inches [default=%default]", metavar = "integer"),
  make_option(c("-s", "--species"), type = "character", default = "",
              help = "Species name for plot title; use Genus_species to italicize [default=%default]", metavar = "Genus_species"),
  make_option(c("--threshold"), type = "double", default = NULL,
              help = "Threshold from 0-1 for including multiple elements in sequence labels; default = dominant element only", metavar = "float"),
  make_option(c("--stacked"), action = "store_true", default = FALSE,
              help = "Generate a sorted stacked bar plot of physical lengths for elements [default=%default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$busco)) {
  stop("You must provide --busco full_table.tsv", call. = FALSE)
}
if (is.null(opt$Stevens) && is.null(opt$Muller)) {
  stop("You must provide either --Stevens or --Muller mapping file.", call. = FALSE)
}
if (!is.null(opt$Stevens) && !is.null(opt$Muller)) {
  stop("Please provide only ONE of --Stevens or --Muller, not both.", call. = FALSE)
}
if (!is.null(opt$threshold) && (opt$threshold < 0 || opt$threshold > 1)) {
  stop("--threshold must be between 0 and 1.", call. = FALSE)
}

scheme <- if (!is.null(opt$Stevens)) "Stevens" else "Muller"
mapFile <- if (!is.null(opt$Stevens)) opt$Stevens else opt$Muller

message("Using element scheme: ", scheme)
message("Mapping file: ", mapFile)

# Read element mapping and standardize the first two columns.
elementDict <- read_tsv(mapFile, col_types = c(col_character(), col_character()), show_col_types = FALSE)
colnames(elementDict)[1:2] <- c("Orthogroups", "Element")
elementDict <- elementDict %>% select(Orthogroups, Element)

# Read BUSCO full table.
busco <- suppressWarnings(read_tsv(
  opt$busco,
  col_names = c("Busco_id", "Status", "Sequence",
                "start", "end", "strand", "Score", "Length",
                "OrthoDB_url", "Description"),
  col_types = c("ccciicdicc"),
  comment = "#",
  show_col_types = FALSE
))

windwSize <- opt$windowSize
minimumGenesPerSequence <- opt$minimumGenesPerSequence
spName <- opt$species
if (grepl("_", spName)) {
  spName <- paste0("*", sub("_", " ", spName), "*")
}

# Merge BUSCO hits with element mapping.
fbusco_all <- busco %>%
  filter(!Status %in% c("Missing")) %>%
  left_join(elementDict, by = c("Busco_id" = "Orthogroups")) %>%
  mutate(
    Element = ifelse(is.na(Element), "-", Element),
    stPos = start
  ) %>%
  filter(Element != "-")

if (nrow(fbusco_all) == 0) {
  stop("No BUSCO rows matched the supplied element mapping.", call. = FALSE)
}

# Keep only sequences with enough mapped BUSCO genes for plotting and output.
seq_keep <- fbusco_all %>%
  count(Sequence, name = "mapped_busco_genes") %>%
  filter(mapped_busco_genes >= minimumGenesPerSequence) %>%
  pull(Sequence)

fbusco <- fbusco_all %>% filter(Sequence %in% seq_keep)

if (nrow(fbusco) == 0) {
  stop("No sequences passed --minimumGenesPerSequence. Lower -m or check the BUSCO/mapping files.", call. = FALSE)
}

# Windowed counts per sequence/element.
elementSummary <- fbusco %>%
  group_by(Sequence, Element) %>%
  filter(!is.na(stPos)) %>%
  mutate(
    ints = if (max(stPos, na.rm = TRUE) > windwSize) {
      as.numeric(as.character(cut(
        stPos,
        breaks = seq(0, max(stPos, na.rm = TRUE) + windwSize, windwSize),
        labels = seq(windwSize, max(stPos, na.rm = TRUE) + windwSize, windwSize),
        include.lowest = TRUE
      )))
    } else {
      NA_real_
    }
  ) %>%
  filter(!is.na(ints)) %>%
  count(ints, Element) %>%
  ungroup()

# Aggregate counts and proportions per sequence/element.
elementAggregated <- fbusco %>%
  count(Sequence, Element, name = "n") %>%
  group_by(Sequence) %>%
  mutate(
    total_genes = sum(n),
    proportion = n / total_genes
  ) %>%
  ungroup()

# Coordinate statistics per Sequence/Element.
coordinateStats <- fbusco %>%
  group_by(Sequence, Element) %>%
  summarise(
    n_genes    = n(),
    seq_length = max(end, na.rm = TRUE),
    start_mean = mean(start, na.rm = TRUE),
    end_mean   = mean(end, na.rm = TRUE),
    start_sd   = sd(start, na.rm = TRUE),
    end_sd     = sd(end, na.rm = TRUE),
    start_min  = min(start, na.rm = TRUE),
    start_max  = max(start, na.rm = TRUE),
    end_min    = min(end, na.rm = TRUE),
    end_max    = max(end, na.rm = TRUE),
    .groups = "drop"
  )

elementAggregatedStats <- elementAggregated %>%
  left_join(coordinateStats, by = c("Sequence", "Element"))

# Optional GC content and real FASTA sequence lengths.
get_gc_content <- function(dna_string) {
  Biostrings::letterFrequency(dna_string, letters = c("G", "C"), as.prob = FALSE) %>%
    sum() -> gc
  n <- Biostrings::letterFrequency(dna_string, letters = "N", as.prob = FALSE)
  total <- length(dna_string) - n
  ifelse(total > 0, gc / total, NA_real_)
}

genome_summary <- tibble(Note = "No genome FASTA provided")
if (!is.null(opt$fasta)) {
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("The Biostrings package is required when using --fasta.", call. = FALSE)
  }

  fasta <- Biostrings::readDNAStringSet(opt$fasta)
  seq_lengths <- Biostrings::width(fasta)
  gc_vals <- sapply(fasta, get_gc_content)
  valid_gc <- !is.na(gc_vals)
  gc_counts <- Biostrings::letterFrequency(fasta[valid_gc], letters = c("G", "C"))
  n_counts <- Biostrings::letterFrequency(fasta[valid_gc], letters = "N")

  gc_data <- tibble(
    Sequence = names(fasta),
    fasta_seq_length = as.numeric(seq_lengths),
    GC_content = as.numeric(gc_vals)
  )

  elementAggregatedStats <- elementAggregatedStats %>%
    left_join(gc_data, by = "Sequence")

  genome_summary <- tibble(
    total_length = sum(seq_lengths),
    num_sequences = length(fasta),
    mean_seq_length = mean(seq_lengths),
    median_seq_length = median(seq_lengths),
    min_seq_length = min(seq_lengths),
    max_seq_length = max(seq_lengths),
    sd_seq_length = sd(seq_lengths),
    genome_gc_content = sum(gc_counts[, "G"] + gc_counts[, "C"]) /
      sum(seq_lengths[valid_gc] - n_counts[, "N"])
  )
}

# Dominant/mixed element labels.
if (!is.null(opt$threshold)) {
  dominantElement <- elementAggregated %>%
    group_by(Sequence) %>%
    filter(proportion >= opt$threshold) %>%
    summarise(
      DominantElement = paste0(sort(unique(Element)), collapse = "/"),
      .groups = "drop"
    )
} else {
  dominantElement <- elementAggregated %>%
    group_by(Sequence) %>%
    slice_max(n, with_ties = FALSE) %>%
    ungroup() %>%
    select(Sequence, Element) %>%
    rename(DominantElement = Element)
}

# Color palettes.
if (scheme == "Stevens") {
  cols <- c(
    "E"    = "#bcb8a7",
    "A"    = "#224160",
    "G"    = "#55a4ab",
    "D"    = "#604a6e",
    "H"    = "#a6913b",
    "B"    = "#68724f",
    "F"    = "#2b1b37",
    "C"    = "#fec735",
    "LG10" = "#fdd3cf",
    "X"    = "#cd4949"
  )
} else {
  cols <- c(
    "A" = "#FF0000",
    "B" = "#8A2BE2",
    "C" = "#0000FF",
    "D" = "#008000",
    "E" = "#FFFF00",
    "F" = "#FFA500"
  )
}

usedEls <- intersect(names(cols), unique(fbusco$Element))
unknownEls <- setdiff(unique(fbusco$Element), names(cols))
cols <- cols[usedEls]
if (length(unknownEls) > 0) {
  cols <- c(cols, setNames(hue_pal()(length(unknownEls)), unknownEls))
}

# Write summary tables.
out_prefix <- tools::file_path_sans_ext(opt$tableOutput)
write_tsv(elementAggregatedStats, opt$tableOutput)
write_tsv(elementAggregated, paste0(out_prefix, "_aggregated.tsv"))
write_tsv(genome_summary, paste0(out_prefix, "_genome_stats.tsv"))
message("Wrote detailed table: ", opt$tableOutput)
message("Wrote aggregated table: ", paste0(out_prefix, "_aggregated.tsv"))
message("Wrote genome stats table: ", paste0(out_prefix, "_genome_stats.tsv"))

# Windowed/facet plot.
if (nrow(elementSummary) > 0) {
  elementSummary <- elementSummary %>%
    left_join(dominantElement, by = "Sequence") %>%
    mutate(
      x_center = ints - windwSize,
      xmin = x_center - (windwSize / 2),
      xmax = x_center + (windwSize / 2)
    )

  plElements <- elementSummary %>%
    mutate(scaffold_f = factor(Sequence, levels = mixedsort(unique(Sequence)))) %>%
    ggplot(aes(fill = Element, y = n, x = ints - windwSize)) +
    facet_grid(scaffold_f ~ ., switch = "y") +
    geom_bar(position = "stack", stat = "identity") +
    ggtitle(spName) +
    theme_classic() +
    scale_y_continuous(breaks = scales::pretty_breaks(4), position = "right") +
    scale_x_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
    scale_fill_manual(values = cols) +
    guides(fill = guide_legend(ncol = 1, title = scheme)) +
    theme(
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      strip.text.y.left = element_text(angle = 0),
      text = element_text(size = 9),
      plot.title = ggtext::element_markdown(),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )

  ggsave(opt$outPlot, plElements, width = opt$width, height = opt$height, bg = "white")
  message("Wrote facet plot: ", opt$outPlot)
} else {
  warning("No sequence had positions exceeding --windowSize, so the windowed/facet plot was skipped.")
}

# Stacked physical barplot per chromosome.
if (opt$stacked) {
  sequence_coords <- busco %>%
    filter(Sequence %in% elementAggregated$Sequence) %>%
    group_by(Sequence) %>%
    summarise(max_gene_coordinate = max(end, na.rm = TRUE), .groups = "drop")

  # Use FASTA lengths for stacked bar heights when available; otherwise use max BUSCO coordinate.
  if (!is.null(opt$fasta) && exists("gc_data")) {
    sequence_coords <- sequence_coords %>%
      left_join(gc_data %>% select(Sequence, fasta_seq_length), by = "Sequence") %>%
      mutate(chrom_length = ifelse(!is.na(fasta_seq_length), fasta_seq_length, max_gene_coordinate))
  } else {
    sequence_coords <- sequence_coords %>% mutate(chrom_length = max_gene_coordinate)
  }

  df <- elementAggregated %>%
    left_join(sequence_coords %>% select(Sequence, chrom_length), by = "Sequence") %>%
    mutate(
      segment_height = proportion * chrom_length,
      physical_proportion = segment_height / chrom_length
    ) %>%
    group_by(Sequence) %>%
    arrange(Sequence, physical_proportion) %>%
    mutate(
      cum_bp = cumsum(segment_height),
      ymin_bp = lag(cum_bp, default = 0),
      ymax_bp = cum_bp
    ) %>%
    ungroup() %>%
    left_join(dominantElement, by = "Sequence") %>%
    mutate(SequenceRenamed = paste0(Sequence, " [", DominantElement, "]"))

  x_levels <- mixedsort(unique(df$SequenceRenamed))

  df <- df %>% mutate(x = as.numeric(factor(SequenceRenamed, levels = x_levels)))

  plElementsStacked <- ggplot(df, aes(
    xmin = x - 0.4,
    xmax = x + 0.4,
    ymin = ymin_bp,
    ymax = ymax_bp,
    fill = Element
  )) +
    geom_rect(color = NA) +
    scale_x_continuous(
      breaks = seq_along(x_levels),
      labels = x_levels,
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      labels = scales::label_number(scale = 1e-6, suffix = "M"),
      expand = c(0, 0),
      breaks = scales::pretty_breaks(n = 4)
    ) +
    scale_fill_manual(values = cols) +
    labs(
      title = paste(spName),
      x = "Sequence",
      y = "Chromosome Length (bp)",
      fill = scheme
    ) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
      text = element_text(size = 9),
      plot.title = ggtext::element_markdown(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )

  ggsave(
    opt$stackedOutPlot,
    plElementsStacked,
    width = opt$width,
    height = opt$height,
    dpi = 300,
    bg = "white"
  )
  message("Wrote stacked plot: ", opt$stackedOutPlot)
}
