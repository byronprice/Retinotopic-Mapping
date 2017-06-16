#!/bin/bash -l

# 24 hour time limit
#$ -l h_rt=12:00:00
#$ -N RetinoMaps
#$ -j y
#$ -o RetinoMaps.txt
#$ -l mem_total=16G
#$ -m ea
#$ -M byron.h.price@gmail.com


module load matlab/2017a

matlab -nodisplay -singleCompThread -r "addpath /usr3/graduate/bhprice;MapRetinotopy(34092,20170319);exit"
