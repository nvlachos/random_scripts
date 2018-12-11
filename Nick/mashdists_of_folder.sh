#!/bin/sh -l

#$ -o mashdist.out
#$ -e mashdist.err
#$ -N mashdist
#$ -cwd
#$ -q all.q

for i in ${1}/*.fasta; do
	query=$(basename $i | cut -d'.' -f1)
	for j in ${2}/*.fasta; do
		ref=$(basename $j | cut -d'.' -f1)
		mash dist ${i} ${j} > ${1}/${query}_${ref}.dist
	done
	cat ${1}/*.dist > ${1}/${query}.dists
done
