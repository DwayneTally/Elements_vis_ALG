Rscript element_vis.R -b ./Example/BUSCO/full_table.tsv -c ./data/stevens.tsv -s T.castaneum -o T.castaneum.pdf --stacked --stackedOutPlot T.castaneum_stacked.pdf --threshold 0.1

Element_vis.R can make the regular barplots, stacked barplots, and also give the user a threshold option. The chromosome labels have been redone to reflect the dominant stevens element in the chormosome [DEFUALT]
 and assigns that stevens element to the chromosome label (in normal barplot png). In stacked barplots, it behaves the same execpt the chromosome label is <original chromosome sequence[stevens element]>. 
The color has also been redone to reflect the color scheme in Ryan et al (2024) paper. 
The threshold option takes numbers between 0-1 from the user and assigns that value as the threshold for the ratio of stevens elements per chromosome and renames the label based on that. 
I.E. if one chromosome is 90% Stevens B and 10% Stevens X and the user set the threshold to 0.1 then the new label stevens element label for that chromosme will be B/X.  

Rscript Element_calculations.R -b ./Example/BUSCO/full_table.tsv -c ./data/stevens.tsv -f ./Example/fna/Tribolium_castaneum_GCF_000002335.3_genomic.fna

Stevens_calculaitons is a new helper script to do some quick statistical analysis on the genome. It takes in the same parameters as Stevens_vis but with the added option for fasta files. It outputs a multi sheet excel file.
The first sheet shows length, number of sequences, average, min, max, and sd length of the sequences (in bp, so far) and GC%. The second page is breaks down the squences by giving you the total number of genes in the sequence and
how many genes that are associated with a specific stevens element are within that sequence, and the overall ratio of which stevens elements make up that chromosomes. The last page is a more detailed look of the second page where it gives number of genes, sequence length, start mean, end mean, start sd, end sd, start min, start max, end min, end max.
