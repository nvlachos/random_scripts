#!/bin/sh -l

#$ -o mashdist.out
#$ -e mashdist.err
#$ -N mashdist
#$ -cwd
#$ -q short.q

ml Mash/2.0

query=$(basename $1 | cut -d'.' -f1)
for j in ${2}/*.fasta; do
	ref=$(basename $j | cut -d'.' -f1)
	mash dist ${1} ${j} > ${2}/${query}_${ref}.dist
done
cat ${2}/*.dist > ${2}/${query}.dists
sort -k3 -n -o "${2}/${query}.dists" "${2}/${query}.dists"
