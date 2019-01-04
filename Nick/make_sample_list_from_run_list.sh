#!/bin/bash -l

#$ -o do_all.out
#$ -e do_all.err
#$ -N rekrak
#$ -cwd
#$ -q all.q

#Import the config file with shortcuts and settings
. ./config.sh

#
# Script to make a list of all isolates in our database
#
# Usage ./make_sample_list_from_run_list.sh path_and_name_of_sample_output_file path_and_name_of_directory_output_file
#
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to make_sample_list_from_run_list.sh, exiting"
	exit 1
elif [[ -z "${1}" ]]; then
	echo "Empty Genus_species supplied to make_sample_list_from_run_list.sh, exiting"
	exit 1
# Gives the user a brief usage and help section if requested with the -h option argument
elif [[ "${1}" = "-h" ]]; then
	echo "Usage is ./make_sample_list_from_run_list.sh	path_and_name_of_SAMPLE_list_output_file path_and_name_of_directory_list_output_file"
	echo "Output is saved to ${1} & ${2}"
	exit 0
fi

dir_list=""

if [[ ! -z "${2}" ]]; then
	dir_list="${2}"
	> "${dir_list}"
fi

sample_list="${1}"
> "${sample_list}"

for path in ${processed}/*; do
    [ -d "${path}" ] || continue # if not a directory, skip
    run_ID="$(basename "${path}")"
	if [[ "${dir_list}" != "" ]]; then
		echo "${run_ID}" >> "${dir_list}"
	fi
	if [[ -f ${path}/${run_id}_list.txt ]]; then
		cat "${path}/${run_id}_list.txt" >> "${sample_list}"
	else
		echo "${run_id} does not have a list file"
	fi
	#for isolate in $path/*; do
	#	[ -d "${isolate}" ] || continue # if not a directory, skip
	#	isolate_name="$(basename "${isolate}")"
	#	echo "${run_ID}/${isolate_name}" >> "${sample_list}"
	done
done

exit
