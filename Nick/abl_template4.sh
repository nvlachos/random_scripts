#!/bin/sh -l

#$ -o abl-template3.out
#$ -e abl-template3.err
#$ -N abl-template3
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
. ./config.sh
#Import the module file that loads all necessary mods
#. "${mod_changers}/pipeline_mods"

#
# Usage ./act_by_list_template.sh path_to_list
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./act_by_list_template.sh path_for_list_file"
	exit 0
fi

# Loop through and act on each sample name in the passed/provided list
while IFS= read -r var; do
	sample_name=$(echo "${var}" | cut -d'/' -f2 | tr -d '[:space:]')
	project=$(echo "${var}" | cut -d'/' -f1 | tr -d '[:space:]')

	SIZES1=$(md5sum "${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R1_001.fastq.gz" | cut -d' ' -f1)
	SIZEA1=$(md5sum "${processed}/AdrianMissingQuaisars/${sample_name}/FASTQs/${sample_name}_R1_001.fastq.gz" | cut -d' ' -f1)
	SIZES2=$(md5sum "${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R2_001.fastq.gz" | cut -d' ' -f1)
	SIZEA2=$(md5sum "${processed}/AdrianMissingQuaisars/${sample_name}/FASTQs/${sample_name}_R2_001.fastq.gz" | cut -d' ' -f1)

	if [[ "${SIZES1}" == "${SIZEA1}" ]]; then
		R1="SAME"
	else
		R1="${sample_name}:NOT_SAME (${SIZES1}:${SIZEA1})"
	fi

	if [[ "${SIZES2}" == "${SIZEA2}" ]]; then
		R2="SAME"
	else
		R2="${sample_name}:NOT_SAME (${SIZES2}:${SIZEA2})"
	fi

	echo "${R1}:${R2}"

done < "${1}"

echo "All isolates completed"
global_end_time=$(date "+%m-%d-%Y @ %Hh_%Mm_%Ss")
#Script exited gracefully (unless something else inside failed)
printf "%s %s" "Act_by_list_template.sh has completed check of snd MLSTs" "${global_end_time}" | mail -s "act_by_list complete" nvx4@cdc.gov
exit 0
