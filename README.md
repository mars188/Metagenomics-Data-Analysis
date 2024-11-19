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
mkdir -p /scratch/$USER/metagenomics_workshop
mkdir -p /scratch/$USER/metagenomics_workshop/hands-on
mkdir -p /scratch/$USER/metagenomics_workshop/hands-on/clean_reads
cd /scratch/$USER/metagenomics_workshop/hands-on

```
# Step 0: Copying the data and creating the required directories 
Run the following command on your terminal. 

```
mkdir -p clean_reads/redsea1 clean_reads/redsea2 clean_reads/redsea3 
cp /scratch/ma5877/metagenomics_workshop/redsea.sh .

ln -s /scratch/ma5877/metagenomics_workshop/data/analysis/redsea1/metawrap_qc/final*fastq clean_reads/redsea1/
ln -s /scratch/ma5877/metagenomics_workshop/data/analysis/redsea2/metawrap_qc/final*fastq clean_reads/redsea2/
ln -s /scratch/ma5877/metagenomics_workshop/data/analysis/redsea3/metawrap_qc/final*fastq clean_reads/redsea3/

```

After performing the above action, you will find the following files in your directory

The `redsea.sh` script for hands-on training and 
the `metagenomic_D1.yml` is a complete workflow for day 1

To access jubail web:

```
https://ood.hpc.abudhabi.nyu.edu
```
# Step 1: Hands-on training with Kraken2 for taxonomic classification
Here is how to load the kraken2 module. 

```
module purge
module load all gencore/2
module load kraken2/2.1.2
```

Run the `kraken2 --help` to pull the detailed flags and then fill the empty spaces with asterik (*) in the script below.
 ```
kraken2 \
--* /scratch/Reference_Genomes/In_house/Metagenomic/kraken2/krakendb_bacteria/ \
--* \
--* kraken2/redsea1_unclassified#.fastq \
clean_reads/redsea1/final_pure_reads_1.fastq \
clean_reads/redsea1/final_pure_reads_2.fastq \
--* \
--* kraken2/redsea1_k2report.txt \
--* 24 \
--* kraken2/redsea1.kraken

kraken2 \
--* /scratch/Reference_Genomes/In_house/Metagenomic/kraken2/krakendb_bacteria/ \
--* \
--* kraken2/redsea2_unclassified#.fastq \
clean_reads/redsea2/final_pure_reads_1.fastq \
 clean_reads/redsea2/final_pure_reads_2.fastq \
--* \
--* kraken2/redsea2_k2report.txt \
--* 24 \
--* kraken2/redsea2.kraken


kraken2 \
--* /scratch/Reference_Genomes/In_house/Metagenomic/kraken2/krakendb_bacteria/ \
--* \
--* kraken2/redsea3_unclassified#.fastq \
clean_reads/redsea3/final_pure_reads_1.fastq \
clean_reads/redsea3/final_pure_reads_2.fastq \
--* \
--* kraken2/redsea3_k2report.txt \
--* 24 \
--* kraken2/redsea3.kraken
```
# Step 2: Running the Bracken for abundance estimation from kraken2 data at Phylum and Genus level

Here is how to load bracken and create the required directories
```
module purge
module load all gencore/2
module load bracken/2.8

mkdir -p bracken
mkdir -p bracken/phylum
mkdir -p bracken/genus
```

Run the `bracken --help` to pull the detailed flags and then fill the empty spaces with asterik (*) in the script below.
```
bracken \
-* /scratch/Reference_Genomes/In_house/Metagenomic/kraken2/krakendb_bacteria \
-* kraken2/redsea1_k2report.txt \
-* bracken/phylum/redsea1 \
-* 90 \
-* P \
-* 5 > bracken/phylum/P_log.txt

bracken \
-* /scratch/Reference_Genomes/In_house/Metagenomic/kraken2/krakendb_bacteria \
-* kraken2/redsea1_k2report.txt \
-* bracken/genus/redsea1 \
-* 90 \
-* G \
-* 5 > bracken/genus/G_log.txt

bracken \
-* /scratch/Reference_Genomes/In_house/Metagenomic/kraken2/krakendb_bacteria \
-* kraken2/redsea2_k2report.txt \
-* bracken/phylum/redsea2 \
-* 90 \
-* P \
-* 5 > bracken/phylum/P_log.txt

bracken \
-* /scratch/Reference_Genomes/In_house/Metagenomic/kraken2/krakendb_bacteria \
-* kraken2/redsea2_k2report.txt \
-* bracken/genus/redsea2 \
-* 90 \
-* G \
-* 5 > bracken/genus/G_log.txt

bracken \
-* /scratch/Reference_Genomes/In_house/Metagenomic/kraken2/krakendb_bacteria \
-* kraken2/redsea3_k2report.txt \
-* bracken/phylum/redsea3 \
-* 90 \
-* P \
-* 5 > bracken/phylum/P_log.txt

bracken \
-* /scratch/Reference_Genomes/In_house/Metagenomic/kraken2/krakendb_bacteria \
-* kraken2/redsea3_k2report.txt \
-* bracken/genus/redsea3 \
-* 90 \
-* G \
-* 5 > bracken/genus/G_log.txt
```
This removes the extra intermediate file not needed for further analysis. 
```
rm kraken2/*genuses.txt
rm kraken2/*phylums.txt
```
# Step 3: Merging the bracken files and visualising the abundance results 
Run the following script to create files and visualize the results. 

```
cd /scratch/$USER/metagenomics_workshop/hands-on
mkdir -p merge_brackenP
mkdir -p merge_brackenG
mkdir -p merge_bracken

cp bracken/phylum/redsea* merge_brackenP
cd merge_brackenP
python /scratch/gencore/ma5877/workflows/Tax_data_convert_kraken2.py \
-i redsea* \
-c 1 \
-o ../merge_bracken/All_data_bracken_merged_P.tsv

cd ../merge_bracken
sed s/Clade_name/Phylum/ All_data_bracken_merged_P.tsv > All_data_bracken_merged_P_new.tsv
cd ../
cp bracken/genus/redsea* merge_brackenG
cd merge_brackenG

python /scratch/gencore/ma5877/workflows/Tax_data_convert_kraken2.py \
-i redsea* \
-c 1 \
-o ../merge_bracken/All_data_bracken_merged_G.tsv

cd ../merge_bracken

sed s/Clade_name/Genus/ All_data_bracken_merged_G.tsv > All_data_bracken_merged_G_new.tsv

module purge
module load gencore Miniconda3/4.7.10
source activate /scratch/gencore/conda3/envs/rbase_4.2.1
cp /scratch/ma5877/metagenomics_workshop/scripts/plot_bracken_* .
R CMD BATCH --no-save plot_bracken_P.R
R CMD BATCH --no-save plot_bracken_G.R
conda deactivate
```
# Step 4: Metaphlan - For taxonomic classification of the reads, relative abundance and samples comparison with heatmap

Here is how to load metaphlan. 
```
mkdir metaphlan4
module purge
module load all gencore/2
module load metaphlan/4.0.1
cd /scratch/$USER/metagenomics_workshop/day1
```
Run the `metaphlan --help` to look at the detailed flags used for metaphlan. 

```
metaphlan \
clean_reads/redsea1/final_pure_reads_1.fastq,clean_reads/redsea1/final_pure_reads_2.fastq \
--input_type fastq \
--bowtie2out metaphlan4/redsea1.bowtie2out.txt \
--nproc 24 \
-o metaphlan4/redsea1_profile.txt

metaphlan \
clean_reads/redsea2/final_pure_reads_1.fastq,clean_reads/redsea2/final_pure_reads_2.fastq \
--input_type fastq \
--bowtie2out metaphlan4/redsea2.bowtie2out.txt \
--nproc 24 \
-o metaphlan4/redsea2_profile.txt

metaphlan \
clean_reads/redsea3/final_pure_reads_1.fastq,clean_reads/redsea3/final_pure_reads_2.fastq \
--input_type fastq \
--bowtie2out metaphlan4/redsea3.bowtie2out.txt \
--nproc 24 \
-o metaphlan4/redsea3_profile.txt
```
Now make the plots for visulising taxonomy at phylum and genus level. 

```
mkdir metaphlan4_plots
cp metaphlan4/*_profile.txt metaphlan4_plots
cd metaphlan4_plots
merge_metaphlan_tables.py *_profile.txt > merged_abundance_table.txt

grep -E "p__|clade" merged_abundance_table.txt | grep -v "c__" | sed "s/^.*|//g" | sed "s/clade_name[0-9]*-//g" > merged_abundance_table_phylum.txt
grep -E "g__|clade" merged_abundance_table.txt | grep -v "s__" | sed "s/^.*|//g" | sed "s/clade_name[0-9]*-//g" > merged_abundance_table_genus.txt

module purge
module load gencore Miniconda3/4.7.10
source activate /scratch/gencore/conda3/envs/hclust2_v99

hclust2.py \
-i merged_abundance_table_phylum.txt \
-o metaphlan4_abundance_heatmap_P.png \
--f_dist_f braycurtis \
--s_dist_f braycurtis \
--cell_aspect_ratio 1 \
-l \
--flabel_size 8 \
--slabel_size 8 \
--max_flabel_len 100 \
--max_slabel_len 100 \
--minv 0.1 \
--dpi 1200 \
--ftop 20

hclust2.py \
-i merged_abundance_table_genus.txt \
-o metaphlan4_abundance_heatmap_G.png \
--f_dist_f braycurtis \
--s_dist_f braycurtis \
--cell_aspect_ratio 0.3 \
-l \
--flabel_size 8 \
--slabel_size 8 \
--max_flabel_len 100 \
--max_slabel_len 100 \
--minv 0.1 \
--dpi 1200 \
--ftop 50
```
# Submitting the workflow with Biosails
We need to run two main commands of biosails as given below. 

```
cd /scratch/$USER/metagenomics_workshop

module purge
module load gencore
module load gencore_biosails

biox run -w metagenomic_D1.yml -o day1.sh

hpcrunner.pl submit_jobs -i day1.sh --project day1_run
```
