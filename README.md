# Mytilus SNParray
In this project we developed a 60K SNP array for four species of blue mussels  _Mytilus edulis, M. galloprovincialis, M. trossulus and M. chilensis 

publication: https://onlinelibrary.wiley.com/doi/full/10.1111/eva.13552

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


dedup and rename

```
#!/bin/bash
#SBATCH -D .
#SBATCH -p highmem
#SBATCH --time=05:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH -A Research_Project-189125
#SBATCH --job-name=master_dedup_merged
#SBATCH --error=master_dedup_merged.err.txt
#SBATCH --output=master_dedup_merged.out.txt
#SBATCH --export=All
#SBATCH --array=0-20%20
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=jn378@exeter.ac.uk # email address

## Change the directories and add the sampleID metadata file ##

module load picard/2.6.0-Java-1.8.0_131

samples=/gpfs/ts0/projects/Research_Project-189125/snp/metadata/alexi-merged.txt

bam_path=/gpfs/ts0/projects/Research_Project-189125/snp/7.merged/

## In array ##
insampleID_array=( `cat $samples |  cut -f 1` )
insampleID=$bam_path/${insampleID_array[(($SLURM_ARRAY_TASK_ID))]}

## Out array
outsampleID_array=( `cat $samples | cut -f 1` )
outsampleID=$bam_path/${outsampleID_array[(($SLURM_ARRAY_TASK_ID))]}

## Mark duplicates in bam files ##

java -Xmx10g -jar $EBROOTPICARD/picard.jar MarkDuplicates I=$insampleID.sorted.rg.bam REMOVE_DUPLICATES=true \
VALIDATION_STRINGENCY=LENIENT AS=true METRICS_FILE=$outsampleID.sorted.dups.metrics.txt OUTPUT=$outsampleID.sorted.dups.b
am TMP_DIR=/gpfs/ts0/scratch/jn378/tmp

### Index the Deduped merged bam files


java -Xmx10g -Djava.io.tmpdir=/gpfs/ts0/scratch/jn378/tmp -jar $EBROOTPICARD/picard.jar BuildBamIndex \
I=$outsampleID.sorted.dups.bam VALIDATION_STRINGENCY=LENIENT
(base) [jn378@login02 7.merged]$
(base) [jn378@login02 7.merged]$
(base) [jn378@login02 7.merged]$ ^C
(base) [jn378@login02 7.merged]$ ls
?  add-readgroups.sh  dedup-jenny.sh  dedup2.sh  merge_samples.sh  rename.sh
(base) [jn378@login02 7.merged]$ more rename.sh
for file in *_merged.sorted.rg.bam
do
    mv -i "${file}" "${file/_merged.sorted.rg.bam/.sorted.rg.bam}"
done
(base) [jn378@login02 7.merged]$ more rename.sh
for file in *_merged.sorted.rg.bam
do
    mv -i "${file}" "${file/_merged.sorted.rg.bam/.sorted.rg.bam}"
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

Qaulimap
```
#!/bin/bash
#SBATCH -D .
#SBATCH -p pq
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH -A Research_Project-189125
#SBATCH --job-name=qualimap
#SBATCH --error=qualimap.err.txt
#SBATCH --output=qualimap.out.txt
#SBATCH --export=All
#SBATCH --array=0-200%20

module purge
module load Java/1.8.0_92

simpleID_array=( `cat $metadata | cut -f 1` )
simpleID=${simpleID_array[(($SLURM_ARRAY_TASK_ID))]}

bam_path=/gpfs/ts0/projects/Research_Project-189125/snp/8.dedup/
samples=/gpfs/ts0/projects/Research_Project-189125/snp/metadata/SNP-array-metadata.txt


## In array ##
insampleID_array=( `cat $samples |  cut -f 1` )
insampleID=$bam_path/${insampleID_array[(($SLURM_ARRAY_TASK_ID))]}
## Out array
outsampleID_array=( `cat $samples | cut -f 1` )
outsampleID=$bam_path/${outsampleID_array[(($SLURM_ARRAY_TASK_ID))]}

qualimap=/gpfs/ts0/projects/Research_Project-189125/software/qualimap_v2.2.1
$qualimap/qualimap bamqc -bam ${insampleID}.sorted.dups.bam -outdir ${outsampleID}_qualimap --java-mem-size=32G -nt 2
```


merge scaffolds in vcf

```
#merge scaffol vcfs
#here we merged the vcf's generated for all invidual scaffolds and merge them in one main vcf file

while read -r name; do cat "./9.scf_dir/$name.vcf"; done < ./galo_genome_2/genome.bed | vcffirstheader | vcfstreamsort -w
 1000 | vcfuniq > merged_all.vcf

done
```
#filter for monomorphic flanking probes

```
awk '{print NR,$1,$2,two_before_current,one_before_current; two_before_current=one_before_current; one_before_current=$0}' <(grep -v '^#' snp_depth_filtered2.
recode.vcf | awk '{print $1,$2}') | awk '{if ((($2==$6)&&(($7+35)>$3)) || (($4==$6) && ($7-$5<35))) print $0"\tno\tcurrent="$7"\tprevious="$5"\tposterior="$3"
\tdiff_current_previos="$7-$5"\tdiff_posterior_current="$3-$7; else print $0"\tyes\tcurrent="$7"\tprevious="$5"\tposterior="$3"\tdiff_current_previos="$7-$5"\
tdiff_posterior_current="$3-$7}' | grep 'yes' | awk 'NF==13' | awk '{print $6, $7, $9}' | sed 's/current=//g' |awk '{if($2==$3) print $1, $2}' > mono_snps.txt
```


filter missing data pops

```
#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p pq # submit to the parallel queue
#SBATCH --time=5:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-189125 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --output=create_pops5.out.txt
#SBATCH --error=create_pops5.err.txt


module load VCFtools/0.1.16-foss-2018b-Perl-5.28.0
#module load Java/1.8.0_144
#module load GATK/3.8-0-Java-1.8.0_144



#update here to your directories
#vcf=freebayes_raw.vcf
#data=/gpfs/ts0/projects/Research_Project-205369/RobEllis/bams
#reference=/gpfs/ts0/projects/Undergraduate_Teaching-205369/genome/GCA_001676915.1_ASM167691v1_genomic.fna

main=/gpfs/ts0/projects/Research_Project-189125/snp/11.vcfs/snp_depth_filtered2.recode.vcf
out=/gpfs/ts0/projects/Research_Project-189125/snp/11.vcfs/pops/stats/
in=/gpfs/ts0/projects/Research_Project-189125/snp/11.vcfs/pops/


#divide into pop specific vcf
vcftools --vcf $main --keep $in/BBC.txt --max-missing 0.5 --out $out/BBC_50.filtered --recode --remove-filtered-all
#vcftools --vcf $main --keep $in/BBC.txt --max-missing 0.8 --out $out/BBC_80.filtered --recode --remove-filtered-all

#vcftools --vcf $main --keep $in/BBH.txt --max-missing 0.8 --out $out/BBH_80.filtered --recode --remove-filtered-all
vcftools --vcf $main --keep $in/BBH.txt --max-missing 0.5 --out $out/BBH_50.filtered --recode --remove-filtered-all

#vcftools --vcf $main --keep $in/CAF.txt --max-missing 0.8 --out $out/CAF_80.filtered --recode --remove-filtered-all
vcftools --vcf $main --keep $in/CAF.txt --max-missing 0.5 --out $out/CAF_50.filtered --recode --remove-filtered-all

#vcftools --vcf $main --keep $in/CC.txt --max-missing 0.8 --out $out/CC_80.filtered --recode --remove-filtered-all
vcftools --vcf $main --keep $in/CC.txt --max-missing 0.5 --out $out/CC_50.filtered --recode --remove-filtered-all

#vcftools --vcf $main --keep $in/CH.txt --max-missing 0.8 --out $out/CH_80.filtered --recode --remove-filtered-all
vcftools --vcf $main --keep $in/CH.txt --max-missing 0.5 --out $out/CH_50.filtered --recode --remove-filtered-all

#vcftools --vcf $main --keep $in/CRS.txt --max-missing 0.8 --out $out/CRS_80.filtered --recode --remove-filtered-all
vcftools --vcf $main --keep $in/CRS.txt --max-missing 0.5 --out $out/CRS_50.filtered --recode --remove-filtered-all

#vcftools --vcf $main --keep $in/CSI.txt --max-missing 0.8 --out $out/CSI_80.filtered --recode --remove-filtered-all
vcftools --vcf $main --keep $in/CSI.txt --max-missing 0.5 --out $out/CSI_50.filtered --recode --remove-filtered-all

#vcftools --vcf $main --keep $in/CV.txt --max-missing 0.8 --out $out/CV_80.filtered --recode --remove-filtered-all
vcftools --vcf $main --keep $in/CV.txt --max-missing 0.5 --out $out/CV_50.filtered --recode --remove-filtered-all

#vcftools --vcf $main --keep $in/DDE.txt --max-missing 0.8 --out $out/DDE_80.filtered --recode --remove-filtered-all
vcftools --vcf $main --keep $in/DDE.txt --max-missing 0.5 --out $out/DDE_50.filtered --recode --remove-filtered-all

#vcftools --vcf $main --keep $in/EX.txt --max-missing 0.8 --out $out/EX_80.filtered --recode --remove-filtered-all
vcftools --vcf $main --keep $in/EX.txt --max-missing 0.5 --out $out/EX_50.filtered --recode --remove-filtered-all

#vcftools --vcf $main --keep $in/BO.txt --max-missing 0.8 --out $out/BO_80.filtered --recode --remove-filtered-all
vcftools --vcf $main --keep $in/BO.txt --max-missing 0.5 --out $out/BO_50.filtered --recode --remove-filtered-all

#vcftools --vcf $main --keep $in/FAL.txt --max-missing 0.8 --out $out/FAL_80.filtered --recode --remove-filtered-all
vcftools --vcf $main --keep $in/FAL.txt --max-missing 0.5 --out $out/FAL_50.filtered --recode --remove-filtered-all

#vcftools --vcf $main --keep $in/Fin.txt --max-missing 0.8 --out $out/Fin_80.filtered --recode --remove-filtered-all
vcftools --vcf $main --keep $in/Fin.txt --max-missing 0.5 --out $out/Fin_50.filtered --recode --remove-filtered-all

#vcftools --vcf $main --keep $in/GA.txt --max-missing 0.8 --out $out/GA_80.filtered --recode --remove-filtered-all
vcftools --vcf $main --keep $in/GA.txt --max-missing 0.5 --out $out/GA_50.filtered --recode --remove-filtered-all

#vcftools --vcf $main --keep $in/GK.txt --max-missing 0.8 --out $out/GK_80.filtered --recode --remove-filtered-all
vcftools --vcf $main --keep $in/GK.txt --max-missing 0.5 --out $out/GK_50.filtered --recode --remove-filtered-all

#vcftools --vcf $main --keep $in/GR.txt --max-missing 0.8 --out $out/GR_80.filtered --recode --remove-filtered-all
vcftools --vcf $main --keep $in/GR.txt --max-missing 0.5 --out $out/GR_50.filtered --recode --remove-filtered-all

#vcftools --vcf $main --keep $in/IC.txt --max-missing 0.8 --out $out/IC_80.filtered --recode --remove-filtered-all
vcftools --vcf $main --keep $in/IC.txt --max-missing 0.5 --out $out/IC_50.filtered --recode --remove-filtered-all

#vcftools --vcf $main --keep $in/LF.txt --max-missing 0.8 --out $out/LF_80.filtered --recode --remove-filtered-all
vcftools --vcf $main --keep $in/LF.txt --max-missing 0.5 --out $out/LF_50.filtered --recode --remove-filtered-all

#vcftools --vcf $main --keep $in/OA.txt --max-missing 0.8 --out $out/OA_80.filtered --recode --remove-filtered-all
vcftools --vcf $main --keep $in/OA.txt --max-missing 0.5 --out $out/OA_50.filtered --recode --remove-filtered-all

#vcftools --vcf $main --keep $in/SBH.txt --max-missing 0.8 --out $out/SBH_80.filtered --recode --remove-filtered-all
vcftools --vcf $main --keep $in/SBH.txt --max-missing 0.5 --out $out/SBH_50.filtered --recode --remove-filtered-all

#vcftools --vcf $main --keep $in/PM.txt --max-missing 0.8 --out $out/PM_80.filtered --recode --remove-filtered-all
vcftools --vcf $main --keep $in/PM.txt --max-missing 0.5 --out $out/PM_50.filtered --recode --remove-filtered-all

#vcftools --vcf $main --keep $in/SW.txt --max-missing 0.8 --out $out/SW_80.filtered --recode --remove-filtered-all
vcftools --vcf $main --keep $in/SW.txt --max-missing 0.5 --out $out/SW_50.filtered --recode --remove-filtered-all

#vcftools --vcf $main --keep $in/VG.txt --max-missing 0.8 --out $out/VG_80.filtered --recode --remove-filtered-all
vcftools --vcf $main --keep $in/VG.txt --max-missing 0.5 --out $out/VG_50.filtered --recode --remove-filtered-all
```


collaps pops

```
#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p pq
#SBATCH --time=6:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-189125 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=2 # specify number of processors per node
#SBATCH --output=collapse50.out.txt
#SBATCH --error=collapse50.err.txt


module load VCFtools/0.1.16-foss-2018b-Perl-5.28.0
module load BCFtools/1.9-foss-2018b


#ls *.gz > merge-list.txt

#bcftools merge -l merge-list.txt -Oz -o merged_80.vcf.gz



#vcftools --gzvcf merged_80.vcf.gz --maf 0.05 --recode --remove-filtered-all --out merged_80_maf_0.05

#vcftools --gzvcf merged_80.vcf.gz --maf 0.01 --recode --remove-filtered-all --out merged_80_maf_0.01

#in this step we are are creating vcf files for each population containing the INTERCEPTING SNPs among populations
bcftools isec -p snps BBC_50.gz BBH_50.gz BO_50.gz CAF_50.gz CC_50.gz CRS_50.gz CSI_50.gz CV_50.gz DDE_50.gz EX_50.gz FAL
_50.gz Fin_50.gz GA_50.gz GK_50.gz GR_50.gz IC_50.gz LF_50.gz OA_50.gz PM_50.gz SBH_50.gz SW_50.gz VG_50.gz CH_50.gz n =
23





done
```

```
#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p pq # submit to the parallel queue
#SBATCH --time=30:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-189125 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --output=maf.out.txt
#SBATCH --error=maf.err.txt


module load VCFtools/0.1.16-foss-2018b-Perl-5.28.0
module load BCFtools/1.9-foss-2018b


#ls *.gz > merge-list.txt

bcftools merge -l merge-list.txt -Oz -o merged_50.vcf.gz



vcftools --gzvcf merged_50.vcf.gz --maf 0.05 --recode --remove-filtered-all --out merged_50_maf_0.05

vcftools --gzvcf merged_50.vcf.gz --maf 0.01 --recode --remove-filtered-all --out merged_50_maf_0.01

done
```

filter unique snps
```
#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p pq
#SBATCH --time=48:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-189125 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=2 # specify number of processors per node
#SBATCH --output=collapse50.out.txt
#SBATCH --error=collapse50.err.txt


module load VCFtools/0.1.16-foss-2018b-Perl-5.28.0
module load BCFtools/1.9-foss-2018b


#ls *.gz > merge-list.txt

#bcftools merge -l merge-list.txt -Oz -o merged_80.vcf.gz



#vcftools --gzvcf merged_80.vcf.gz --maf 0.05 --recode --remove-filtered-all --out merged_80_maf_0.05

#vcftools --gzvcf merged_80.vcf.gz --maf 0.01 --recode --remove-filtered-all --out merged_80_maf_0.01


#in this script we are creating vcf files for each pop, containing the UNIQUE snps in each populations - SNPs that are pr
esent in these pops alone
bcftools isec -p unique_snps/bbc -C BBC_50.gz BBH_50.gz BO_50.gz CAF_50.gz CC_50.gz CRS_50.gz CSI_50.gz CV_50.gz DDE_50.g
z EX_50.gz FAL_50.gz Fin_50.gz GA_50.gz GK_50.gz GR_50.gz IC_50.gz LF_50.gz OA_50.gz PM_50.gz SBH_50.gz SW_50.gz VG_50.gz
 CH_50.gz

```

pop stats

```

#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p pq # submit to the parallel queue
#SBATCH --time=2:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-189125 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --output=pop_stats.out.txt
#SBATCH --error=pop_stats.err.txt
#SBATCH --array=0-25%25


module load VCFtools/0.1.16-foss-2018b-Perl-5.28.0
module load Java/1.8.0_144
module load GATK/3.8-0-Java-1.8.0_144

out=/gpfs/ts0/projects/Research_Project-189125/snp/11.vcfs/pops/stats/pop50/stats-pop50/

ls *.vcf > pops.fofn

simpleID_array=( `cat pops.fofn` )
simpleID=${simpleID_array[(($SLURM_ARRAY_TASK_ID))]}

#divide into pop specific vcf
#vcftools --vcf $main --keep ${simpleID} --out $out/${simpleID}.filtered --recode --remove-filtered-all
#vcftools --gzvcf allele_filtered.recode.vcf.gz --freq2 --out a_freq --max-alleles 2

#calculate mean depth per individual
vcftools --vcf ${simpleID} --site-mean-depth --out $out/${simpleID}
#calculate depth per individual
vcftools --vcf ${simpleID} --depth --out $out/${simpleID}
#calculate quality per SNP call
vcftools --vcf ${simpleID} --site-quality --out $out/${simpleID}
#calculate heterozigosity
vcftools --vcf ${simpleID} --het --out $out/${simpleID}
#calculate missing data per individual
vcftools --vcf ${simpleID} --missing-individual --out $out/${simpleID}
#calculate missing data per site
vcftools --vcf ${simpleID} --missing-site --out $out/${simpleID}
#calulate missing dat per individual
vcftools --vcf ${simpleID} --missing-indv --out $out/${simpleID}
```




###to call final set snps
```
module load BEDTools/2.27.1-foss-2018b

#genome=/gpfs/ts0/projects/Research_Project-189125/snp/galo_genome_2/mgal_unzip_genome.bed
#in=/gpfs/ts0/projects/Research_Project-189125/snp/11.vcfs/snps-array/maxmiss50/merged_mono_flank_maf_0.01_miss0.05.filtered.vcf
#in2=/gpfs/ts0/projects/Research_Project-189125/snp/11.vcfs/snps-array/maxmiss50/LF_50.filtered.recode.vcf

grep -v "^#" trimed_maf_0.01_flank_miss_0.5.vcf | awk -v OFS='\t' '{ print $1"\t"$2"\t"$2}' > snps.bed

#bedtools flank -i snps.bed -g mgal_unzip_genome.bed -b 35 > flank_seq.bed

#at this stage mg_00001 was giving issues (bedtools didnt like that it started at the begining of the genome)
#so I've created a bed file containing the faulty reads of mg_00001 (grep -v mg_00001) 
#created a bed file without the mg_00001 sequences
#ran this file on bedtools
#summed the proper ends of the flanking regions in excel
#combined both files and went to the next step
bedtools getfasta -fi mgal_unzip_genome.fa -bed flank_seq.bed > flank-seq.fasta
#bedtools getfasta -fi /gpfs/ts0/projects/Research_Project-189125/snp/galo_genome_2/mgal_unzip_genome.fa -bed new.bed > flank-seq.fasta 

#cat flank-seq.fasta | sed -n '1~4p' > left.seq.bed
#cat flank-seq.fasta | sed -n '2~4p' > left2.seq.bed
#cat flank-seq.fasta | sed -n '3~4p' > right.seq.bed
#cat flank-seq.fasta | sed -n '4~4p' > right2.seq.bed
#paste left.seq.bed left2.seq.bed right.seq.bed right2.seq.bed | awk '{ print $1"\t"$2"\t"$3"\t"$4 }' > final-flank-seq.bed 


#grep -v "^#" snps-filtered-from-ends.vcf | awk '{ print $1"\t["$4"/"$5"]\t"$1"\t"$2"\t"$4"\t"$5 }' > trial2.txt 
#paste trial2.txt final-flank-seq.bed | awk '{ print $1"\t"$8""$2""$10 }' > list-snps-flanks.txt
#paste list-snps-flanks.txt trial2.txt | awk '{ print $1"\t"$2"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9 }' > excel.txt
#cat excel.txt | awk '{gsub("[0-9]*","",$1)}1'|  awk '{gsub(/mg_/,"myt_"NR,$1); print}' > final-snp-set.txt
#awk -F mytilus '{$1= FS"\t" $1;}1' < final-snp-set.txt > snp-set.txt
#awk '{print $1"\t"$2"\t"$3"\t""NA""\t""Standard""\t" $4 "\t" $5 "\t"$6"\t"$7}' snp-set.txt > snp-set-2.txt
#awk 'BEGIN{print "Organism" "\t" "snpid" "\t" "Seventyonemer" "\t" "Tiling_order" "\t" "Importance" "\t" "Chromosome" "\t" "Position" "\t" "Ref_allele" "\t" "Alt_allele"}1' snp-set-2.txt > snp-set-final-all.txt
```
