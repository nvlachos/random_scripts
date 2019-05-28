#!/bin/sh -l

#$ -o fix_species.out
#$ -e fix_species.err
#$ -N fix_species
#$ -cwd
#$ -q short.q


#
# This will check all files in a directory to see if the species was extracted properly (e.g, no sp.)
#
# Usage is ./fix_species.sh directory_of_fastas extension_to_look_for
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
elif [[ -z "${1}" ]]; then
	echo "Empty directory supplied to $0, exiting"
	exit 1
# Gives the user a brief usage and help section if requested with the -h option argument
elif [[ "${1}" = "-h" ]]; then
	echo "Usage is ./make_scicomp_refseq_mash_sketch.sh	output_directory"
	echo "Output is saved to ${1} & ${2}"
	exit 0
fi


#echo ":${1}:"
if [[ ! -d ${1} ]]; then
	echo "${1} does not exist"
fi

for i in ${1}/*; do
	for j in ${i}/*; do
		genus=$(basename ${j} | cut -d'_' -f1)
		species=$(basename ${j} | cut -d'_' -f2)
		accession=$(basename ${j} | cut -d'_' -f3,4)

			if [[ "${species}" = "sp." ]]; then
				header=$(head -n1 ${j})
				species=$(echo "${header}" | cut -d' ' -f4 | sed 's/,//g')
				echo "${genus}-${species}-${accession}"
				echo "saving ${j} to ${i}/${genus}_${species}_${accession}"
				mv "${j}" "${i}/${genus}_${species}_${accession}"
			fi
		done
done
