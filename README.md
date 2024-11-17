# Metagenomics-Data-Analysis
This page describes the steps to perform the NYUAD Core Bioinformatics hands-on workshop on Metagenomics data analysis and visulizing the results.

During this workshop, participants will learn about the concepts and steps involved in analyzing the metagenomic data obtained from short read high throughput sequencing.

Participants are expected to have some basic knowledge of metagenome/microbiome as well as some biological knowledge on the subject matter.

Although the material and the methods are designed to cater for the NYU Abu Dhabi High Performance Computing environment, it is possible to run this workshop on any other system provided that the neccessary software is installed, and the data is uploaded to that environment. It will also be neccessary to change the input and output directories/files to accomodate such an environment.

We will start off with fastq sequencing files (Illumina paired end short reads), and throughout the course of the workshop, we will learn how to process and analyze the data so that we end up with abundance tables that can be further used for making graphs and visulizing results. Moving forward, we will use the publically available redsea data downloaded from Androsiuk et al. 2022 and can be found here. (https://www.ncbi.nlm.nih.gov/bioproject/PRJNA289734/).

To summarize, we will:

Start off with QC/QT of the reads 

Reference based analysis
Taxonomic annotation with Kraken2 and MetaPhlAn 
Functional annotation with HUMAnN 
Assembly based analysis
Binning, bin refinement, bin taxonomly and functional profiling 

# Connecting to the HPC using a MAC machine #
* Open the "Terminal" app and type ssh NetID@jubail.abudhabi.nyu.edu
* Enter your NYU password when prompted and hit enter
* Once logged in, navigate to your personal "SCRATCH" directory or the any other subdirecotry where you want to run this workshop

# Connecting to the HPC using a Windows machine #
* Open the "Putty" app, and fill out the fields as follows Host name jubail.abudhabi.nyu.edu, Port=22, and then click on "Open".
* Enter your NetId, and your password when prompted.
* Once logged in, navigate to your personal "SCRATCH" directory or the any other subdirecotry where you want to run this workshop.

# Setting up the environment and copying the data #
We will be using the NYUAD High Performance Computing (HPC) cluster for this workshop, however, you can certainly run all of the analysis on any stand alone machine (server, personal laptop/Desktop etc.) provided that you have pre-installed the necessay software packages.

```
mkdir -p /scratch/$USER/TE_workshop 
cd /scratch/$USER/TE_workshop
rsync -avP /scratch/ma5877/TE_workshop/data_files/ .
```

To access jubail web:

```
https://ood.hpc.abudhabi.nyu.edu
```
# The data
After performing the above action, you will find the following files in your directory

The `redsea.sh` for hands-on training and `metagenomic_D1.yml` is a complete workflow for day 1 
squire_fetch: This contains the Repeatmasker, annotation , chromosome fasta files downloaded from UCSC for mouse genome (mm9)
squire_clean: Filtered repeatmasker file for repeats of interest after collapsing overlapping repeats. We obtained these files by preforming squire_clean step from the workflow. We have already performed this step and you don't have to repeat.
trimmed_reads: This directory contains the dataset that we will be using for this workshop, which are quality trimmed FASTQ sequencing files. They are from the Aging epigenomics of male mice and publicly available to download from the Short Read Archive (SRA [https://www.ncbi.nlm.nih.gov/sra]) using the following accessions SRS2248611, SRS2248612, SRS2248613, SRS2248599, SRS2248600, SRS2248601. These data have been downsized from their original dataset. The reason for this is due to time constraints. Running the complete dataset takes much longer than the time we have during this workshop.
yml file: This is the workflow file that we will make use of to run the analysis.
NOTE: There are a total 6 FASTQ files. Three files belong to the young mice (3-month) and other three came from the old mice (29-month). We will use young mice as control group to detect the differential TE expression in old mice.
