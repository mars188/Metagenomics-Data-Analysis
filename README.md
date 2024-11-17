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
