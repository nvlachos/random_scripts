#!/bin/sh -l

#$ -o mmshfolder.out
#$ -e mmshfolder.err
#$ -N mmshfolder
#$ -cwd
#$ -q all.q


module load Mash/2.0

echo ":${1}:"
if [[ ! -d ${1} ]]; then
	echo "${1} does not exist"
fi
command="mash sketch -o ${2}"
#mkdir ${1}/sketches
for j in ${1}/*.fna; do
	 #filename=$(basename ${j})
	 command="${command} ${j}"
done
#mv ${1}/*.msh ${1}/sketches
echo ${command}
${command}
