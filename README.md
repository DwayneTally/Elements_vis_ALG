You need to run BUSCO on your genome assembly using the endopterygota_odb10 dataset. This script takes as input the resulting full_table.tsv.

Element_vis.R can make the regular barplots, stacked barplots, and also give the user a threshold option. The chromosome labels have been redone to reflect the dominant stevens element in the chormosome [DEFUALT]
 and assigns that stevens element to the chromosome label (in normal barplot png). In stacked barplots, it behaves the same execpt the chromosome label is <original chromosome sequence[stevens element]>. 
The color has also been redone to reflect the color scheme in Ryan et al (2024) paper. 
The threshold option takes numbers between 0-1 from the user and assigns that value as the threshold for the ratio of stevens elements per chromosome and renames the label based on that. 
I.E. if one chromosome is 90% Stevens B and 10% Stevens X and the user set the threshold to 0.1 then the new label stevens element label for that chromosme will be B/X.  

NOTE* to run Element_calculations.R shown in the example below, you'd need to unzip the fna file.

Element_calculations is a new helper script to do some quick statistical analysis on the genome. It takes in the same parameters as Stevens_vis but with the added option for fasta files. It outputs a multi sheet excel file.
The first sheet shows length, number of sequences, average, min, max, and sd length of the sequences (in bp) and GC%. The second page is breaks down the squences by giving you the total number of genes in the sequence and
how many genes that are associated with a specific stevens element are within that sequence, and the overall ratio of which stevens elements make up that chromosomes. The last page is a more detailed look of the second page where it gives number of genes, sequence length, start mean, end mean, start sd, end sd, start min, start max, end min, end max.

## Installation
Since this R tool builds off of vis_ALG (https://github.com/pgonzale60/vis_ALG), Please follow their installation guide. Other packages that might be required: optparse, readr, dplyr, scales, ggplot2, gtools,ggtext, Biostrings


## Usage
Rscript element_vis.R -b full_table.tsv -c stevens.tsv -s Genus_species -o output.pdf --stacked --stackedOutPlot output_stacked.pdf --threshold 0.5

Rscript Element_calculations.R -b full_table.tsv -c stevens.tsv -f species.fna


## Examples
Rscript element_vis.R -b ./Example/BUSCO/full_table.tsv -c ./data/stevens.tsv -s T.castaneum -o T.castaneum.pdf --stacked --stackedOutPlot T.castaneum_stacked.pdf --threshold 0.1


Rscript element_vis.R -b ./Example/BUSCO/full_table.tsv -c ./data/stevens.tsv -s T.castaneum -o T.castaneum.pdf --stacked --stackedOutPlot T.castaneum_stacked.pdf --threshold 0.1

Rscript element_vis.R -b ./Example/BUSCO/full_table.tsv -c ./data/stevens.tsv -s T.castaneum -o T.castaneum.pdf --threshold 0.1 -f Tribolium_castaneum_GCF_000002335.3_genomic.fna