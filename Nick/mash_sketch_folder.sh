#!/bin/sh -l

#$ -o mashdist.out
#$ -e mashdist.err
#$ -N mashdist
#$ -cwd
#$ -q short.q

ml Mash/2.0

echo ":${1}:"
if [[ ! -d ${1} ]]; then
	echo "${1} does not exist"
fi
mkdir ${1}/sketches
for j in ${1}/*.fasta; do
	mash sketch ${j}
done
mv ${1}/*.msh ${1}/sketches
