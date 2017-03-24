# VarP
A tool (written in Python) to predict disease class based on gene panel sequencing data. This method has been developed as part of our lab's participation in a CAGI4 challenge. The challenge was to challenge was to determine which of 14 disease classes each of 106 patients has and the corresponding causal variants, given each patient’s gene panel sequencing data for 83 genes.

# Dependencies
* Varant (http://compbio.berkeley.edu/proj/varant/Home.html) - VarP is built on top of Varant tool and therefore users have to install Varant before running VarP.
* Python (>=v2.6)
* Python libraries - numpy, mayplotlib, tabix
* All the VCF files used should be Tabix indexed.

# Details on the code structure
QCModule folder contains all the scripts used to generate QC metrics and plot for the gene panel sequencing data.
* geneCoverageStat.py - Computes the average read depth for all the genes captured in the panel sequencing data
* makeControlVariants.py - Generates a tsv file that contains 1000 Genomes variants found in the given captured region bed file.
* plotGeneCoverage.py - Plot gene coverage for a given gene in the panel sequenencing datat
* variantCallStat.py - Generate variant call related metrices per sample such as - Ti/Tv, hom/het, # of SNVs/Indel (common, rare, novel), # of no call sites, # of low quality sites

VariantPrioritization folder contains the script that prioritizes variants based on minor allele frequency, mutation impact and inhertiance model, and a further embeds an VarP's algorithm to predict disease class for each patients.
* predictDisease.py - Identified potentially causative variants and disease class per sample.

SampleFiles folder contains a sample of all the data files used in VarP. The intention of providing this folder is to make the users aware of the file formates. They may reused some of these files and edit them with their own values.


