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
```
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
```

Align samples with bwa

```
#!/bin/bash
#SBATCH -D .
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -p pq
#SBATCH --time=100:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH -A Research_Project-189125
#SBATCH --job-name=bwa_map_rerun
#SBATCH --export=All
#SBATCH --array=0-70%10
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=jn378@exeter.ac.uk # email address


### Script to run a bwa mem, bam and sort bam
# Usage: bwa_align_m.sh <ref file> <reads location> <alignment location>
# e.g. qsub 2_bwa_align_m.sh /research/lifesci_fraser/NERO/<GENOME> /research/lifesci_fraser/<POP>/Cleaned_reads ~/lustre
/wgs_data/raw_bams/



./bwa-mem2
module load SAMtools/1.9-foss-2018b
module load BCFtools/1.9-foss-2018b
##./bwa-mem2 index ../galo_genome/mgal_01.fa.gz

reference=/gpfs/ts0/home/jn378/mussels/snp/galo_genome/mgal_01.fa.gz
input_reads=/gpfs/ts0/home/jn378/mussels/snp/3.fastp
align=/gpfs/ts0/home/jn378/mussels/snp/7.bwa
bamfiles=/gpfs/ts0/projects/Research_Project-189125/snp/6.bam-mem2/

metadata=/gpfs/ts0/projects/Research_Project-189125/snp/snp-md-rerun.txt


read1=( `cat $metadata | cut -f 4` )
read1_array=$input_reads/${read1[(($SLURM_ARRAY_TASK_ID))]}

read2=( `cat $metadata | cut -f 5` )
read2_array=$input_reads/${read2[(($SLURM_ARRAY_TASK_ID))]}

outbam=( `cat $metadata | cut -f 1` )
out=${outbam[(($SLURM_ARRAY_TASK_ID))]}

echo "reference" $reference
echo "read1" $read1_array
echo "read2" $read2_array
echo "alignment" $align/${out}_unsorted.raw.sam

# ### Align with bwa mem ###

./bwa-mem2 mem -t 4 $reference $read1_array $read2_array > ${align}/${out}_unsorted.raw.sam

# ### Convert bam to sam, sort bam, index, flagstat ##

samtools view -bS ${align}/${out}_unsorted.raw.sam > ../6.bam-mem2/${out}_unsorted.raw.bam
samtools sort ../6.bam-mem2/${out}_unsorted.raw.bam -o ../6.bam-mem2/${out}_sorted.raw.bam
rm ../6.bam-mem2/${out}_unsorted.raw.bam
samtools index ../6.bam-mem2/${out}_sorted.raw.bam
samtools flagstat ../6.bam-mem2/${out}_sorted.raw.bam  > ../6.bam-mem2/${out}_mappingstats.txt

# ## Remove the sam files (I've commented this out for now just in case anything goes wrong with the sorting and indexing
)

# #rm ${align}/*.sam


## To check that the bams are not corrupted, run:

# samtools quickcheck -v *sorted.raw.bam > bad_bams.fofn   && echo 'all ok' || echo 'some files failed check, see bad_bam
s.fofn'
```

add RG tags


```
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
#SBATCH --array=0-200%32

### This script will add readgroups to the bam files based on info contained in the metadata file

module load picard/2.6.0-Java-1.8.0_131

## Change the directories ##

metadata=/gpfs/ts0/projects/Research_Project-189125/snp/rg_alexi.txt
bam_path=/gpfs/ts0/projects/Research_Project-205369/2.backup-bam-jenny/

# Metadata file for each population
# 1 sample_ID
# 2 read1
# 3 read2
# 4 clean1
# 5 clean2
# 6 instrument
# 7 seq_number
# 8 run_number
# 9 flow_cell
# 10 lane
# 11 barcode

simpleID_array=( `cat $metadata | cut -f 1` )
simpleID=${simpleID_array[(($SLURM_ARRAY_TASK_ID))]}

instrument_array=( `cat $metadata | cut -f 6` )
instrument=${instrument_array[(($SLURM_ARRAY_TASK_ID))]}

#seqnum_array=( `cat $metadata | cut -f 7` )
#seqnum=${seqnum_array[(($PBS_ARRAYID))]}

flow_cell_array=( `cat $metadata | cut -f 9` )
flow_cell=${flow_cell_array[(($PBS_ARRAYID))]}

lane_array=( `cat $metadata | cut -f 10` )
lane=${lane_array[(($PBS_ARRAYID))]}

barcode_array=( `cat $metadata | cut -f 11` )
barcode=${barcode_array[(($PBS_ARRAYID))]}

## In array ##
insampleID_array=( `cat $metadata | cut -f 1` )
insampleID=$bam_path/${insampleID_array[(($SLURM_ARRAY_TASK_ID))]}

# ## Out array
outsampleID_array=( `cat $metadata | cut -f 1` )
outsampleID=$bam_path/${outsampleID_array[(($SLURM_ARRAY_TASK_ID))]}

### Run picard tools AddreplaceRGs
#jp commands
#java -Xmx10g -Djava.io.tmpdir=/gpfs/ts0/scratch/jrp228 -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
#I=${insampleID}.sorted.bam O=${outsampleID}.sorted.rg.bam RGID=${simpleID}.${flow_cell}.${lane} RGLB=${simpleID}.${seqnu
m} \
#RGPL=${instrument} RGSM=${simpleID} RGPU=${flow_cell}.${lane}.${barcode}

#mussel commands
java -Xmx10g -Djava.io.tmpdir=/gpfs/ts0/scratch/jn378 -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I=${insampleID
}_sorted.raw.bam O=${outsampleID}.sorted.rg.bam RGID=${simpleID} RGLB=${simpleID}.1 RGPL=${instrument} RGSM=${simpleID} R
GPU=${flow_cell}.${lane}.${barcode}
```

merge bams
```
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
```

cut chunks of genome to call SNPs from

```
#!/bin/bash
#SBATCH -D .
#SBATCH -p pq
#SBATCH --time=00:5:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -A Research_Project-189125
#SBATCH --job-name=freebayes-jenny
#SBATCH --error=freebayes.err.txt
#SBATCH --output=freebayes.out.txt
#SBATCH --export=All
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=jn378@exeter.ac.uk # email address
#SBATCH --array=1-10

module purge
module load freebayes/1.3.1-GCCcore-8.2.0


## Fill in directories
bed=/gpfs/ts0/projects/Research_Project-189125/snp/galo_genome_2/genome.bed_chunk2.bed
data=/gpfs/ts0/projects/Research_Project-189125/snp/8.dedup/
REF=/gpfs/ts0/projects/Research_Project-189125/snp/galo_genome_2/mgal_unzip_genome.fa
scf_dir=/gpfs/ts0/projects/Research_Project-189125/snp/9.scf_dir
TMPDIR=/gpfs/ts0/scratch/jn378

# Use these to divide up the input file into the start and end chunks we want
START_LINE=$(( 100 * $SLURM_ARRAY_TASK_ID - 99 ))
END_LINE=$(( 100 * $SLURM_ARRAY_TASK_ID ))

# Make a chunk of our bedfile with the relevant 100 contigs in it
sed -n -e "${START_LINE},${END_LINE}p" $bed > ${bed}_chunk${SLURM_ARRAY_TASK_ID}.bed
new_bed=${bed}_chunk2${SLURM_ARRAY_TASK_ID}.bed
echo $new_bed
#genome3.bed has all contigs shorter than 10000 removed

#ls *.sorted.dups.bam > bam.fofn

#creates an array out of the scf names
#scf=( `cat $bed | cut -f 1` )
#scf_array=${scf[(($SLURM_ARRAY_TASK_ID))]}

#scf_array=( `cat $bed | cut -f 1` )
#start loop through all scfs
#for i in "${scf_array[@]}"
#do

#for scf in $(cut -f1 $new_bed)
#do

#test freebayes
# Run in parallel over the scf_array
#cut -f1 $new_bed | parallel -j8 --tmpdir $TMPDIR "freebayes --fasta-reference $REF --bam-list bam.fofn -g 1000 -C 4 -r {
} > $scf_dir/{}.vcf"
#echo done vcf

#Tidy up and remove the new bed
#rm -f $new_bed
#done
```

Call variants with freebayes
```
#!/bin/bash
#SBATCH -D .
#SBATCH -p pq
#SBATCH --time=40:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH -A Research_Project-189125
#SBATCH --job-name=notempdirb-jenny
#SBATCH --error=fbnotemp-%A_%a.err.txt
#SBATCH --output=fb-%A_%a.out.txt
#SBATCH --export=All
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=jn378@exeter.ac.uk # email address
#SBATCH --mail-type=fail         # send email if job fails
#SBATCH --array=254

module purge
module load freebayes/1.3.1-GCCcore-8.2.0


## Fill in directories
bed=/gpfs/ts0/projects/Research_Project-189125/snp/galo_genome_2/genome_2.bed
data=/gpfs/ts0/projects/Research_Project-189125/snp/8.dedup/
REF=/gpfs/ts0/projects/Research_Project-189125/snp/galo_genome_2/mgal_unzip_genome.fa
scf_dir=/gpfs/ts0/projects/Research_Project-189125/snp/9.scf_dir
TMPDIR=/gpfs/ts0/scratch/jn378

#we used genome_2 bed (removing first 1000 contigs because I was testing in these first 1000 with another script freebaye
s
#genome_2.bed and genome.bed bed files have filtered contings with coverage over 20% as well as those smaller than 10 000

#resulting in 23693 (or something like that) contigs to be analysed
#this little post is for me to remeber the parameters when writing the methods section :D


# Use these to divide up the input file into the start and end chunks we want
START_LINE=$(( 100  * $SLURM_ARRAY_TASK_ID - 99 ))
END_LINE=$(( 100 * $SLURM_ARRAY_TASK_ID ))

# Make a chunk of our bedfile with the relevant 100 contigs in it
sed -n -e "${START_LINE},${END_LINE}p" $bed > ${bed}_chunk${SLURM_ARRAY_TASK_ID}.bed
new_bed=${bed}_chunk${SLURM_ARRAY_TASK_ID}.bed

#genome3.bed has all contigs shorter than 10000 removed

ls *.sorted.dups.bam > bam.fofn

#creates an array out of the scf names
#scf=( `cat $new_bed | cut -f 1` )
#scf_array=${scf[(($SLURM_ARRAY_TASK_ID))]}

#scf_array=( `cat $bed | cut -f 1` )
#start loop through all scfs
#for i in "${scf_array[@]}"
#do

#for scf in $(cut -f1 $new_bed)
#do

#test freebayes
# Run in parallel over the scf_array
#scf=( `cat $bed | cut -f 1` )
#scf_array=${scf[(($SLURM_ARRAY_TASK_ID))]}

#start loop through all scfs
#for i in "${scf_array[@]}"
#do
#test freebayes

#freebayes --fasta-reference $REF --bam-list bam.fofn -g 1000 -C 4 -r ${scf_array} > $scf_dir/${scf_array}.vcf

cut -f1 $new_bed | parallel -j8  "freebayes --fasta-reference $REF --bam-list bam.fofn -g 1000 -C 4 -r {} > $scf_dir/{}.v
cf"
echo done vcf

#Tidy up and remove the new bed
rm -f $new_bed
```
