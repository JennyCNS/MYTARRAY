# Mytilus SNParray
In this project we developed a 60K SNP array for four species of blue mussels <\i> Mytilus edulis, M. galloprovincialis, M. trossulus and M. chilensis <\i>
# Scripts used for the development of the Blue Mussel 60K SNP array
# All analysis were ran in ISCA - HPC from the UoE 


# 1. Quality check analysis


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

for f in `ls -1 *_r1.fq.gz | sed 's/_r1.fq.gz//' `
do
fastp -i ${f}_r1.fq.gz -I ${f}_r2.fq.gz -o ../3.fastp/${f}_r1_fastp.fq.gz -O ../3.fastp/${f}_r2_fastp.fq.gz -q 20 -l 100
--trim_poly_g
done

