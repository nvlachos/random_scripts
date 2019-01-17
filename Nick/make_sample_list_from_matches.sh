#!/bin/bash -l

#$ -o matcher.out
#$ -e matcher.err
#$ -N matcher
#$ -cwd
#$ -q all.q

#Import the config file with shortcuts and settings
. ./config.sh

#
# Script to make a list of all matching isolates in our database which matches genus species
#
# Usage ./make_sample_list_from_matches.sh Genus_species path_and_name_of_output_file
#
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to make_sample_list_from_matches.sh, exiting"
	exit 1
elif [[ -z "${1}" ]]; then
	echo "Empty Genus_species supplied to make_sample_list_from_matches.sh.sh, exiting"
	exit 1
# Gives the user a brief usage and help section if requested with the -h option argument
elif [[ "${1}" = "-h" ]]; then
	echo "Usage is ./make_sample_list_from_matches.sh	Genus_species_to_match	path_and_name_of_output_file"
	echo "Output is saved to in ${processed}/sample_name/ANI"
	exit 0
elif [ -z "$2" ]; then
	echo "Empty output location supplied to make_sample_list_from_matches.sh "
	exit 1
fi


NOW=$(date "+%m-%d-%Y")
genus=$(echo ${1} | cut -d'_' -f1)
species=$(echo ${1} | cut -d'_' -f2)
sample_list="${2}"
> "${sample_list}"

for path in ${processed}/*; do
    if [ ! -d "${path}" ]; then
			echo "$path is not counted as path"
			continue
		fi
    run_ID="$(basename "${path}")"
    for isolate in $path/*; do
			[ -d "${isolate}" ] || continue # if not a directory, skip
			isolate_name="$(basename "${isolate}")"
			if [[ ! -f "${processed}/${run_ID}/${isolate_name}/${isolate_name}.tax" ]]; then
				"${shareScript}/determine_taxID.sh" "${isolate_name}" "${run_ID}"
			fi
			echo "Trying1 ${processed}/${run_ID}/${isolate_name}/${isolate_name}.tax"
			sample_genus=$(tail -n3 "${processed}/${run_ID}/${isolate_name}/${isolate_name}.tax")
			echo "Trying2 ${sample_genus}"
			sample_genus=$(echo "${smaple_genus}" | head -n1 | cut -d'	' -f2)
			echo "Trying3 ${sample_genus}"
			sample_species=$(tail -n2 "${processed}/${run_ID}/${isolate_name}/${isolate_name}.tax" | head -n1 | cut -d'	' -f2)
			echo "Trying4 ${sample_species}"
			if [[ -z "${sample_species}" ]] || [[ "${sample_genus}" = "Not_assigned" ]] || [[ "${sample_genus}" = "N/A" ]]; then
				rm ${processed}/${run_ID}/${isolate_name}/${isolate_name}.tax
				"${shareScript}/determine_taxID.sh" "${isolate_name}" "${run_ID}"
				echo "Trying5 ${processed}/${run_ID}/${isolate_name}/${isolate_name}.tax"
				sample_genus=$(tail -n2 "${processed}/${run_ID}/${isolate_name}/${isolate_name}.tax")
				echo "Trying6 ${sample_genus}"
				sample_genus=$(echo "${smaple_genus}" | head -n1 | cut -d'	' -f2)
				echo "Trying7 ${sample_genus}"
				sample_species=$(tail -n1 "${processed}/${run_ID}/${isolate_name}/${isolate_name}.tax" | cut -d'	' -f2)
				echo "Trying8 ${sample_species}"
			fi
			if [[ "${1,,}" == "${sample_genus,,}_${sample_species,,}" ]]; then
				echo "Match-${1} at ${run_ID}/${isolate_name}"
				echo "${run_ID}/${isolate_name}" >> "${sample_list}"
			else
				echo "No match of ${1} to ${sample_genus,,} ${sample_species} from ${isolate_name} in ${run_ID}/${isolate_name}"
			fi
		done
done

exit
