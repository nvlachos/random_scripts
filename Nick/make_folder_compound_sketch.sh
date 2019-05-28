#!/bin/sh -l

#$ -o mmshfolder.out
#$ -e mmshfolder.err
#$ -N mmshfolder
#$ -cwd
#$ -q short.q


module load Mash/2.0

#
# Script to make a list of all matching isolates in our database with match genus species
#
# Usage ./make_folder_compound_sketch.sh output_folder folder_of_fasta_to_make_sketch_from
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
elif [[ -z "${1}" ]]; then
	echo "Empty output directory supplied to $0, exiting"
	exit 1
# Gives the user a brief usage and help section if requested with the -h option argument
elif [[ "${1}" = "-h" ]]; then
	echo "Usage is ./make_folder_compound_sketch.sh	output_directory folder_of_fastas_to_sketch"
	echo "Output is saved to in ${processed}/sample_name/ANI"
	exit 0
elif [ -z "$2" ]; then
	echo "Empty fasta folder location supplied to $0"
	exit 1
fi

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
