#!/bin/bash
#SBATCH --partition=short
#SBATCH --job-name=deduper3
#SBATCH --output=deduper3.out
#SBATCH --error=deduper3.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --mem=20G

ml purge
ml samtools/1.5 
module load easybuild intel/2017a python3/3.6.1

samtools sort -O sam /projects/bgmp/shared/deduper/Dataset3.sam > /projects/bgmp/adoe/deduper/Dataset3.sam

python3 deduper.py -s /projects/bgmp/adoe/deduper/Dataset3.sam -u /projects/bgmp/adoe/deduper/STL96.txt