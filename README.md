# Mytilus SNParray
In this project we developed a 60K SNP array for four species of blue mussels  _Mytilus edulis, M. galloprovincialis, M. trossulus and M. chilensis _

## Scripts used for the development of the Blue Mussel 60K SNP array
All analysis were ran in ISCA - HPC from the UoE 


## 1. Quality check analysis

```
#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p pq # submit to the parallel queue
#SBATCH --time=100:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-189125 # research project to submit under
#SBATCH --nodes=8 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=jn378@exeter.ac.uk # email address
# Commands
module load fastp/0.19.7-foss-2018b
module load FastQC/0.10.1-Java-1.7.0_80

for f in `ls -1 *_r1.fq.gz | sed 's/_r1.fq.gz//' `
do
fastp -i ${f}_r1.fq.gz -I ${f}_r2.fq.gz -o ../3.fastp/${f}_r1_fastp.fq.gz -O ../3.fastp/${f}_r2_fastp.fq.gz -q 20 -l 100
--trim_poly_g
done

fastqc ../1.raw-data/*.gz -o .
#get stats of the files with seqkit

seqkit stats *.gz

#multiqc visualisation

module load MultiQC/1.2-intel-2017b-Python-2.7.14

multiqc .

#run fq_aggregate script

# Run this script in a directory containing zip files from fastqc. It aggregates images of each type in individual folder
s
# So looking across data is quick.
# http://www.danielecook.com/aggregate-fastqc-reports/
# usage:
# 1. Copy this script to directory containing fastqc reports
# 2. Run 'bash fastq_aggregate.sh'
# Run this script in a directory containing zip files from fastqc. It aggregates images of each type in individual folder
s
# So looking across data is quick.

zips=`ls *.zip`

for i in $zips; do
    unzip -o $i &>/dev/null;
done

fastq_folders=${zips/.zip/}

rm -rf fq_aggregated # Remove aggregate folder if present
mkdir fq_aggregated

# Rename Files within each using folder name.
for folder in $fastq_folders; do
    folder=${folder%.*}
    img_files=`ls ${folder}/Images/*png`;
    for img in $img_files; do
        img_name=$(basename "$img");
        img_name=${img_name%.*}
        new_name=${folder};
        mkdir -p fq_aggregated/${img_name};
        mv $img fq_aggregated/${img_name}/${folder/_fastqc/}.png;
    done;
done;


# Concatenate Summaries
for folder in $fastq_folders; do
    folder=${folder%.*}
    cat ${folder}/summary.txt >> fq_aggregated/summary.txt
done;

# Concatenate Statistics
for folder in $fastq_folders; do
    folder=${folder%.*}
    head -n 10 ${folder}/fastqc_data.txt | tail -n 7 | awk -v f=${folder/_fastqc/} '{ print $0 "\t" f }' >> fq_aggregated
/statistics.txt
    rm -rf ${folder}
done

```
##2. Aligned trimmed reads to genome
Align genome
'''
#!/bin/bash
#SBATCH -D .
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -p highmem
#SBATCH --time=100:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH -A Research_Project-189125
#SBATCH --job-name=bwa_repeat
#SBATCH --export=All
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=jn378@exeter.ac.uk # email address


### Script to run a bwa mem, bam and sort bam
# Usage: bwa_align_m.sh <ref file> <reads location> <alignment location>
# e.g. qsub 2_bwa_align_m.sh /research/lifesci_fraser/NERO/<GENOME> /research/lifesci_fraser/<POP>/Cleaned_read$



./bwa-mem2
module load SAMtools/1.9-foss-2018b
module load BCFtools/1.9-foss-2018b
./bwa-mem2 index ../galo_genome/mgal_01.fa.gz

done
'''


merge bams
'''
#!/bin/bash
#SBATCH -D .
#SBATCH -p pq
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH -A Research_Project-189125
#SBATCH --job-name=master_add_readgroups
#SBATCH --error=master_add_readgroups.err.txt
#SBATCH --output=master_add_readgroups.out.txt
#SBATCH --export=All
#SBATCH --array=0-20%20

### This script will add readgroups to the bam files based on info contained in the metadata file

module load SAMtools/1.9-foss-2018b

## Change the directories ##

metadata=/gpfs/ts0/projects/Research_Project-189125/snp/SNP-array-metadata-alexi-merge.txt
bam_path=/gpfs/ts0/projects/Research_Project-189125/snp/6.bam-mem2/



outbam=( `cat $metadata | cut -f 1` )
out=${outbam[(($SLURM_ARRAY_TASK_ID))]}



samtools merge ${out}_merged_sorted.raw.bam ${out}_sorted.raw.bam ${out}_a_sorted.raw.bam



done
'''
