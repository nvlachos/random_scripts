#!/bin/bash -l

#$ -o samp_list.out
#$ -e sampl_list.err
#$ -N sample_list
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
. ./config.sh

#
# Script to make a list of all isolates in our database
#
# Usage ./make_sample_list_from_run_list.sh path_and_name_of_sample_output_file path_and_name_of_directory_output_file -r (optional to remove FAILED isolates)
#
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
	#if [[ "${dir_list}" != "" ]]; then
		echo "${run_ID}" >> "${dir_list}"
	#fi
	if [[ -f "${path}/${run_ID}_list.txt" ]]; then
		cat "${path}/${run_ID}_list.txt" >> "${sample_list}"
	elif [[ -f "${path}/${run_ID}_list_ordered.txt" ]]; then
		cat "${path}/${run_ID}_list_ordered.txt" >> "${sample_list}"
	else
		for isolate in ${path}/*; do
			[ -d "${isolate}" ] || continue # if not a directory, skip
	 		isolate_name="$(basename "${isolate}")"
	 		echo "${run_ID}/${isolate_name}" >> "${path}/${run_ID}_list.txt"
		done
		if [[ -f "${path}/${run_ID}_list.txt" ]]; then
			cat "${path}/${run_ID}_list.txt" >> "${sample_list}"
		fi
fi
done

if [[ "${3}" = "-r" ]]; then
	mv "${sample_list}" "${sample_list}_temp"
	> "${sample_list}"
	while IFS=read -r var  || [ -n "$var" ]; do
		if [[ "${var}" = *"_FAILED" ]]; then
			echo "Removing ${var} from the list"
		else
			echo "${var}" >> "${sample_list}"
		fi
	done < "${sample_list}_temp"
	rm "${sample_list}_temp"
fi



exit
