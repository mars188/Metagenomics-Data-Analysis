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
# Step 0: Copying the data and creating the required directories 
Run the following command on your terminal. 

```
mkdir -p clean_reads/redsea1 clean_reads/redsea2 clean_reads/redsea3 
cp /scratch/ma5877/metagenomics_workshop/redsea.sh .
cp /scratch/ma5877/metagenomics_workshop/metagenomic_D1.yml .

ln -s /scratch/ma5877/metagenomics_workshop/data/analysis/redsea1/metawrap_qc/final*fastq clean_reads/redsea1/
ln -s /scratch/ma5877/metagenomics_workshop/data/analysis/redsea2/metawrap_qc/final*fastq clean_reads/redsea2/
ln -s /scratch/ma5877/metagenomics_workshop/data/analysis/redsea3/metawrap_qc/final*fastq clean_reads/redsea3/

```

After performing the above action, you will find the following files in your directory

The `redsea.sh` script for hands-on training and 
the `metagenomic_D1.yml` is a complete workflow for day 1

# Step 1: Running the Kraken2 for taxonomic classification
Here is how to load the kraken2 module. 

 module purge
 module load all gencore/2
 module load kraken2/2.1.2

 Run the `kraken2 --help` to pull the detailed flags and then fill the emply spaces in the script below.
 ```
 mkdir kraken2
 kraken2 \
 -- /scratch/Reference_Genomes/In_house/Metagenomic/kraken2/krakendb_bacteria/ \
 -- \
 -- kraken2/redsea1_unclassified#.fastq \
 clean_reads/redsea1/final_pure_reads_1.fastq \
 clean_reads/redsea1/final_pure_reads_2.fastq \
 -- \
 -- kraken2/redsea1_k2report.txt \
 -- 24 \
 -- kraken2/redsea1.kraken

 kraken2 \
 -- /scratch/Reference_Genomes/In_house/Metagenomic/kraken2/krakendb_bacteria/ \
 -- \
 -- kraken2/redsea2_unclassified#.fastq \
 clean_reads/redsea2/final_pure_reads_1.fastq \
 clean_reads/redsea2/final_pure_reads_2.fastq \
 -- \
 -- kraken2/redsea2_k2report.txt \
 -- 24 \
 -- kraken2/redsea2.kraken


 kraken2 \
 -- /scratch/Reference_Genomes/In_house/Metagenomic/kraken2/krakendb_bacteria/ \
 -- \
 -- kraken2/redsea3_unclassified#.fastq \
 clean_reads/redsea3/final_pure_reads_1.fastq \
 clean_reads/redsea3/final_pure_reads_2.fastq \
 -- \
 -- kraken2/redsea3_k2report.txt \
 -- 24 \
 -- kraken2/redsea3.kraken
```
