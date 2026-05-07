#!/usr/bin/env Rscript

Sys.setenv(TZ = 'America/New_York')

library(optparse)
library(readr)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(scales))
library(ggplot2)
library(gtools)
library(ggtext)

option_list = list(
  make_option(c("-b", "--busco"), type = "character", default = NULL,
              help = "BUSCO full_table.tsv file", metavar = "file.tsv"),
  make_option(c("-n", "--nigon"), type = "character", default = "gene2Nigon_busco20200927.tsv.gz",
              help = "BUSCO id assignment to Nigons [default=%default]", metavar = "file.tsv"),
  make_option(c("-a", "--Muller"), type = "character", default = NULL,
              help = "BUSCO id assignment to Mullers", metavar = "file.tsv"),
  make_option(c("-c", "--Stevens"), type = "character", default = NULL,
              help = "BUSCO id assignment to Stevens", metavar = "file.tsv"),
  make_option(c("-w", "--windowSize"), type = "integer", default = 5e5,
              help = "Window size (bp) to bin the BUSCO genes [default=%default]", metavar = "integer"),
  make_option(c("-m", "--minimumGenesPerSequence"), type = "integer", default = 15,
              help = "Sequences with fewer than this number of BUSCO genes will not be shown [default=%default]", metavar = "integer"),
  make_option(c("-o", "--outPlot"), type = "character", default = "Elements.jpeg",
              help = "Output image file name for the facet barplot [default=%default]", metavar = "file"),
  make_option(c("--stackedOutPlot"), type = "character", default = "stacked_barplot.png",
              help = "Output image file name for the stacked barplot [default=%default]", metavar = "file"),
  make_option(c("--tableOutput"), type = "character", default = "element_breakdown.tsv",
              help = "Output TSV file with element breakdown [default=%default]", metavar = "file.tsv"),
  make_option(c("--height"), type = "integer", default = 6,
              help = "Height of plot (in inches) [default=%default]", metavar = "integer"),
  make_option(c("--width"), type = "integer", default = 5,
              help = "Width of plot (in inches) [default=%default]", metavar = "integer"),
  make_option(c("-s", "--species"), type = "character", default = "",
              help = "Species name (for plot title, italicized) [default=%default]", metavar = "Genus_species"),
  make_option(c("--threshold"), type = "double", default = NULL,
              help = "Threshold (0–1) for including multiple elements in label; default = dominant element only.", metavar = "float"),
  make_option(c("--stacked"), action = "store_true", default = FALSE,
              help = "Generate a sorted stacked bar plot of physical lengths for elements [default=%default]")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if (is.null(opt$Stevens) && is.null(opt$Muller)) {
  stop("You must provide either --Stevens or --Muller mapping file.")
}
if (!is.null(opt$Stevens) && !is.null(opt$Muller)) {
  stop("Please provide only ONE of --Stevens or --Muller, not both.")
}

if (!is.null(opt$Stevens)) {
  scheme <- "Stevens"
  mapFile <- opt$Stevens
} else {
  scheme <- "Muller"
  mapFile <- opt$Muller
}

message("Using element scheme: ", scheme)
message("Mapping file: ", mapFile)

# Read element mapping (generic)
elementDict <- read_tsv(mapFile, col_types = c(col_character(), col_character()))

# Make sure columns are named consistently
colnames(elementDict)[1:2] <- c("Orthogroups", "Element")

# Read BUSCO full table
busco <- suppressWarnings(read_tsv(
  opt$busco,
  col_names = c("Busco_id", "Status", "Sequence",
                "start", "end", "strand", "Score", "Length",
                "OrthoDB_url", "Description"),
  col_types = c("ccciicdicc"),
  comment = "#"
))

windwSize <- opt$windowSize
minimumGenesPerSequence <- opt$minimumGenesPerSequence
spName <- opt$species
if (grepl("\\.", spName)) {
  spName <- paste0("*", sub("_", " ", spName), "*")
}

# Merge BUSCO with mapping and filter to assigned elements
fbusco <- busco %>%
  filter(!Status %in% c("Missing")) %>%
  left_join(elementDict, by = c("Busco_id" = "Orthogroups")) %>%
  mutate(
    Element = ifelse(is.na(Element), "-", Element),
    stPos   = start
  ) %>%
  filter(Element != "-")

# Windowed counts per sequence/element
elementSummary <- fbusco %>%
  group_by(Sequence, Element) %>%
  filter(!is.na(stPos)) %>%
  mutate(
    ints = if (max(stPos, na.rm = TRUE) > windwSize) {
      as.numeric(as.character(
        cut(
          stPos,
          breaks = seq(0, max(stPos, na.rm = TRUE), windwSize),
          labels = seq(windwSize, max(stPos, na.rm = TRUE), windwSize)
        )
      ))
    } else {
      NA_real_
    }
  ) %>%
  filter(!is.na(ints)) %>%
  count(ints, Element) %>%
  ungroup()

# Aggregate counts and proportions per sequence/element
elementAggregated <- elementSummary %>%
  group_by(Sequence, Element) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  group_by(Sequence) %>%
  mutate(
    total_genes = sum(n),
    proportion  = n / total_genes
  ) %>%
  ungroup()

print(head(elementAggregated))

# Coordinate stats per Sequence–Element
coordinateStats <- fbusco %>%
  filter(Element != "-") %>%
  group_by(Sequence, Element) %>%
  summarise(
    n_genes    = n(),
    seq_length = max(end,    na.rm = TRUE),
    start_mean = mean(start, na.rm = TRUE),
    end_mean   = mean(end,   na.rm = TRUE),
    start_sd   = sd(start,   na.rm = TRUE),
    end_sd     = sd(end,     na.rm = TRUE),
    start_min  = min(start,  na.rm = TRUE),
    start_max  = max(start,  na.rm = TRUE),
    end_min    = min(end,    na.rm = TRUE),
    end_max    = max(end,    na.rm = TRUE),
    .groups    = "drop"
  )

elementAggregatedStats <- elementAggregated %>%
  left_join(coordinateStats, by = c("Sequence", "Element"))

if (!is.null(opt$threshold)) {
  dominantElement <- elementAggregated %>%
    group_by(Sequence) %>%
    filter(proportion >= opt$threshold) %>%
    summarise(
      DominantElement = paste0(sort(unique(Element)), collapse = "/"),
      .groups         = "drop"
    )
} else {
  dominantElement <- elementAggregated %>%
    group_by(Sequence) %>%
    slice_max(n, with_ties = FALSE) %>%
    ungroup() %>%
    select(Sequence, Element) %>%
    rename(DominantElement = Element)
}

#  Color palettes

if (scheme == "Stevens") {

  # Fixed palette for Stevens elements
  cols <- c(
    "E"    = "#bcb8a7",  # LG2
    "A"    = "#224160",  # LG3
    "G"    = "#55a4ab",  # LG4
    "D"    = "#604a6e",  # LG5
    "H"    = "#a6913b",  # LG6
    "B"    = "#68724f",  # LG7
    "F"    = "#2b1b37",  # LG8
    "C"    = "#fec735",  # LG9
    "LG10" = "#fdd3cf",  # LG10
    "X"    = "#cd4949"   # LGX
  )

  # Trim to only elements present
  usedEls <- intersect(names(cols), unique(fbusco$Element))
  cols <- cols[usedEls]

} else if (scheme == "Muller") {

  # Fixed palette for Muller elements
  cols <- c(
    "A" = "#FF0000",  # red
    "B" = "#8A2BE2",  # violet (BlueViolet)
    "C" = "#0000FF",  # blue
    "D" = "#008000",  # green
    "E" = "#FFFF00",  # yellow
    "F" = "#FFA500"   # orange
  )

  # Trim to only elements present
  usedEls <- intersect(names(cols), unique(fbusco$Element))
  cols <- cols[usedEls]

} else {

  # Fallback (debug) 
  elementLevels <- sort(unique(fbusco$Element))
  cols <- setNames(hue_pal()(length(elementLevels)), elementLevels)
}


# facet plot

elementSummary <- elementSummary %>%
  left_join(dominantElement, by = "Sequence")%>%
  mutate(
    x_center = ints - windwSize,
    xmin = x_center - (windwSize/2),
    xmax = x_center + (windwSize/2)
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
    axis.title.y      = element_blank(),
    axis.title.x      = element_blank(),
    strip.text.y.left = element_text(angle = 0),
    text              = element_text(size = 9),
    plot.title        = ggtext::element_markdown(),
    panel.border      = element_blank(),
    panel.background  = element_rect(fill = "white", color = NA),
    plot.background   = element_rect(fill = "white", color = NA)
  )

ggsave(opt$outPlot, plElements, width = opt$width, height = opt$height)

# Stacked physical barplot (per chromosome)

if (opt$stacked) {
  sequence_coords <- busco %>%
    filter(Sequence %in% elementAggregated$Sequence) %>%
    group_by(Sequence) %>%
    summarise(max_gene_coordinate = max(end), .groups = "drop")

  elementAggregated2 <- elementAggregated %>%
    left_join(sequence_coords, by = "Sequence") %>%
    mutate(segment_height = proportion * max_gene_coordinate)

  df <- elementAggregated2 %>%
    rename(chrom_length = max_gene_coordinate) %>%
    mutate(physical_proportion = segment_height / chrom_length) %>%
    group_by(Sequence) %>%
    # sort so biggest segment is at the top of the stack
    arrange(Sequence, desc(-physical_proportion)) %>%
    mutate(
      cum_bp  = cumsum(segment_height),
      ymin_bp = lag(cum_bp, default = 0),
      ymax_bp = cum_bp
    ) %>%
    ungroup() %>%
    left_join(dominantElement, by = "Sequence") %>%
    mutate(
      # label = sequence + dominant element on the SAME line
      SequenceRenamed = paste0(Sequence, " [", DominantElement, "]")
    )

  # define a *single* consistent ordering for the x axis
  x_levels <- mixedsort(unique(df$SequenceRenamed))

  df <- df %>%
    mutate(
      x = as.numeric(factor(SequenceRenamed, levels = x_levels))
    )

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
      x     = "Sequence",
      y     = "Chromosome Length (bp)",
      fill  = scheme
    ) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # <<< SLANTED LABELS HERE >>>
      axis.text.x      = element_text(angle = 45, vjust = 1, hjust = 1),
      text             = element_text(size = 9),
      plot.title       = ggtext::element_markdown(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA)
    )

  ggsave(
    opt$stackedOutPlot,
    plElementsStacked,
    width  = opt$width,
    height = opt$height,
    dpi    = 300,
    bg     = "white"
  )
}

