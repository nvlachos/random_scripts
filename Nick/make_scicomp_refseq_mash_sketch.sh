#!/bin/sh -l

#$ -o mmshref.out
#$ -e mmshref.err
#$ -N mmshref
#$ -cwd
#$ -q short.q


module load Mash/2.0

#
# This will create  a sketch of everything available in the refseq database found on scicomp
#
# Usage is ./make_scicomp_refseq_mash_sketch.sh Output_directory
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
elif [[ -z "${1}" ]]; then
	echo "Empty Genus_species supplied to $0, exiting"
	exit 1
# Gives the user a brief usage and help section if requested with the -h option argument
elif [[ "${1}" = "-h" ]]; then
	echo "Usage is ./make_scicomp_refseq_mash_sketch.sh	output_directory"
	echo "Output is saved to ${1}"
	exit 0
fi

command="mash sketch -o ${1}"
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
#mv ${1}/*.msh ${local_DBs}/aniDB/
echo ${command}
echo ${command} > "${1}/make_refeq_command.txt"
${command}
