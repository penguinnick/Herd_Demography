#!/bin/bash
#------------------------
#SBATCH --account=gratis
#------------------------
#SBATCH --partition=cpu-invest
#SBATCH --qos=job_cpu_preemptable
#SBATCH --job-name="preempt Payne sheep SensAnalysis"
#SBATCH --time=6:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-user=nicholas.triozzi@unibe.ch
#SBATCH --mail-type=all

# Your code below this line
# cd Herd_Demography_UBELIX/
module load Workspace_Home
module load R
# export R_LIBS=$HOME/R/x86_64-pc-linux-gnu-library/4.4/
Rscript --vanilla R/XX_SENSITIVITYpayne_sheep.R