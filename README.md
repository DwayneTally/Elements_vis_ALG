# Elements Visualization Tool

This tool builds off of the original creator's tool (https://github.com/pgonzale60/vis_ALG), where it processes BUSCO results and originally assigned Nigon elements (Now includes Stevens or Muller elements) to generate chromosome visualizations and summary statistics for genome assemblies.

## Overview

`element_vis.R` uses BUSCO gene coordinates and an element mapping file to infer chromosome composition and generate visual summaries of conserved chromosome elements.

The script requires a BUSCO `full_table.tsv` generated using an appropriate lineage dataset (`endopterygota_odb10` for beetles).

Supported element schemes:
- **Stevens elements** (`-c` / `--Stevens`)
- **Muller elements** (`-a` / `--Muller`)
- **Nigon elements** (`-n` / `--Nigon`)
The script can generate:

- **barplots**
  - Original barplot distribution of the chromosomes.
  - Shows chromosome composition in genomic windows.
  - Chromosomes are labeled according to their dominant element by default.
  - If a threshold is supplied, chromosomes can be labeled with multiple elements.

- **Stacked chromosome-length barplots (`--stacked`)**
  - Shows total chromosome physical length "stacked" by inferred element composition.
  - Labels are formatted as:
    ```
    chromosome_name [dominant_element]
    ```
    or, when threshold labeling is enabled:
    ```
    chromosome_name [B/X]
    ```

- **TSV summary files**
  - Detailed per-sequence element breakdown
  - Aggregated chromosome composition
  - Genome-wide FASTA statistics (when FASTA provided)

## Installation

Create a conda environment:

```bash
ENV_NAME=element_vis

conda create -n $ENV_NAME \
    -c conda-forge \
    -c bioconda \
    r-base \
    r-optparse \
    r-readr \
    r-dplyr \
    r-scales \
    r-ggplot2 \
    r-gtools \
    r-ggtext \
    bioconductor-biostrings \
    -y

conda activate $ENV_NAME
```

Required R packages:
- optparse
- readr
- dplyr
- scales
- ggplot2
- gtools
- ggtext
- Biostrings (only required for `--fasta`)

## Threshold Behavior

The `--threshold` option accepts values between `0–1`.

By default:
- only the **dominant element** is used for chromosome labels.

With `--threshold`, any element meeting that proportion threshold is included.

Example:

If a chromosome contains:
- 90% Stevens B
- 10% Stevens X

Then:

```bash
--threshold 0.1
```

will label that chromosome:

```text
B/X
```



## Output Files

### 1. `element_breakdown.tsv`
Detailed per-sequence / per-element statistics.

Columns include:
- Sequence
- Element
- Number of mapped BUSCO genes (`n`)
- Total genes on that sequence
- Proportion of each element
- Number of genes assigned to that element
- BUSCO coordinate-derived sequence length
- Start/end coordinate summary statistics:
  - mean
  - standard deviation
  - min
  - max



### 2. `element_breakdown_aggregated.tsv`
Simplified per-sequence element composition table.

Contains:
- Sequence
- Element
- Number of mapped BUSCO genes
- Total mapped genes
- Proportion of each element



### 3. `element_breakdown_genome_stats.tsv`
Generated only when `--fasta` or `-f` is supplied.

Genome-wide assembly summary statistics:
- total assembly length
- number of sequences
- mean sequence length
- median sequence length
- min sequence length
- max sequence length
- standard deviation of sequence lengths
- genome GC content



## FASTA Option (`-f`, `--fasta`)

Providing a FASTA file enables:

- genome-wide assembly statistics
- GC content calculation
- use of true FASTA sequence lengths for stacked chromosome barplots

If no FASTA is supplied:
- stacked plots use BUSCO maximum gene coordinates as chromosome lengths
- genome summary statistics are not calculated



## Usage

### Stevens Elements
```bash
Rscript element_vis.R \
    -b full_table.tsv \
    -c stevens.tsv \
    -s Genus_species
```

### Muller Elements
```bash
Rscript element_vis.R \
    -b full_table.tsv \
    -a muller.tsv \
    -s Genus_species
```

### Generate stacked plot
```bash
Rscript element_vis.R \
    -b full_table.tsv \
    -c stevens.tsv \
    --stacked \
    --stackedOutPlot output_stacked.png
```

### Use FASTA for genome stats + chromosome lengths
```bash
Rscript element_vis.R \
    -b full_table.tsv \
    -c stevens.tsv \
    --fasta species.fna
```

### Threshold labeling
```bash
Rscript element_vis.R \
    -b full_table.tsv \
    -c stevens.tsv \
    --threshold 0.1
```



## Examples

### Standard Stevens visualization
```bash
Rscript Scripts/element_vis.R \
    -b ./Example/BUSCO/full_table.tsv \
    -c ./data/stevens.tsv \
    -s T.castaneum \
    -o T.castaneum.jpeg
```

### Stacked plot with threshold labeling
```bash
Rscript Scripts/element_vis.R \
    -b ./Example/BUSCO/full_table.tsv \
    -c ./data/stevens.tsv \
    -s T.castaneum \
    --stacked \
    --stackedOutPlot T.castaneum_stacked.png \
    --threshold 0.1
```

### FASTA-enabled genome summary
```bash
Rscript Scripts/element_vis.R \
    -b ./Example/BUSCO/full_table.tsv \
    -c ./data/stevens.tsv \
    -s T.castaneum \
    --fasta ./Example/fna/Tribolium_castaneum_GCF_000002335.3_genomic.fna
```

## Options
```bash
  "-b", "--busco", type = "character", default = NULL,
              help = "BUSCO full_table.tsv file", metavar = "file.tsv"),
  "-n", "--nigon", type = "character", default = "gene2Nigon_busco20200927.tsv.gz",
              help = "BUSCO id assignment to Nigons [kept for compatibility; currently unused]", metavar = "file.tsv"),
  "-a", "--Muller", type = "character", default = NULL,
              help = "BUSCO id assignment to Muller elements", metavar = "file.tsv"),
  "-c", "--Stevens", type = "character", default = NULL,
              help = "BUSCO id assignment to Stevens elements", metavar = "file.tsv"),
  "-f", "--fasta", type = "character", default = NULL,
              help = "Optional genome FASTA file for sequence lengths and GC content", metavar = "file.fa"),
  "-w", "--windowSize", type = "integer", default = 5e5,
              help = "Window size (bp) to bin BUSCO genes [default=%default]", metavar = "integer"),
  "-m", "--minimumGenesPerSequence", type = "integer", default = 15,
              help = "Sequences with fewer than this many mapped BUSCO genes are not plotted [default=%default]", metavar = "integer"),
  "-o", "--outPlot", type = "character", default = "Elements.jpeg",
              help = "Output image file name for the windowed/facet barplot [default=%default]", metavar = "file"),
  "--stackedOutPlot", type = "character", default = "stacked_barplot.png",
              help = "Output image file name for the stacked chromosome-length barplot [default=%default]", metavar = "file"),
  "--tableOutput", type = "character", default = "element_breakdown.tsv",
              help = "Output TSV file with detailed element breakdown [default=%default]", metavar = "file.tsv"),
  "--height", type = "integer", default = 6,
              help = "Height of plot in inches [default=%default]", metavar = "integer"),
  "--width", type = "integer", default = 5,
              help = "Width of plot in inches [default=%default]", metavar = "integer"),
  "-s", "--species", type = "character", default = "",
              help = "Species name for plot title; use Genus_species to italicize [default=%default]", metavar = "Genus_species"),
  "--threshold", type = "double", default = NULL,
              help = "Threshold from 0-1 for including multiple elements in sequence labels; default = dominant element only", metavar = "float"),
  "--stacked", action = "store_true", default = FALSE,
              help = "Generate a sorted stacked bar plot of physical lengths for elements [default=%default]")
```
## Notes
  - Genome-wide FASTA statistics tsv file will still be generated even if you don't give a fasta file input. the file will just contain 'Note No genome FASTA provided'
