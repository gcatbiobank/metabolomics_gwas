#! /bin/bash
#$ -N mGWAS
#$ -V
#$ -t 1-777
#$ -tc 20
#$ -o job_logs/sge-output.$TASK_ID.log
#$ -e job_logs/sge-error.$TASK_ID.log
#$ -l h_vmem=16G

module load apps/R-3.5.0
Rscript mGWAS.R
module purge
