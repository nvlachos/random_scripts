#!/bin/sh -l

#$ -o mmshfolder.out
#$ -e mmshfolder.err
#$ -N mmshfolder
#$ -cwd
#$ -q all.q


module load Mash/2.0

#echo ":${1}:"
today=$(date '+%Y%m%d')
/scicomp/reference/public-references/refseq-bacteria
#if [[ ! -d ${1} ]]; then
#	echo "${1} does not exist"
#fi
command="mash sketch -o ${2}"
counter=0
#mkdir ${1}/sketches
for i in /scicomp/reference/public-references/refseq-bacteria/*; do
	echo "----- ${i} ------"
	for j in ${i}/*; do
		#filename=$(basename ${j})
		#command="${command} ${j}"
		if [[ -d "${j}" ]]; then 
			assembly_name=$(basename ${j})
			#echo "${j}/${assembly_name}_genomic.fna.gz"
			if [[ -f "${j}/${assembly_name}_genomic.fna.gz" ]]; then
				#echo "${counter}-${j}/${assembly_name}_genomic.fna.gz"
				command="${command} ${j}/${assembly_name}_genomic.fna.gz"
			fi
			counter=$(( counter + 1 ))
		fi
	done
done
#mv ${1}/*.msh ${share}/DBs/aniDB/
echo ${command}
echo ${command} > "${share}/make_refeq_command.txt"
${command}
