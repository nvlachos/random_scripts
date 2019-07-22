#!/bin/sh -l

#$ -o abl-dmlst.out
#$ -e abl-dmlst.err
#$ -N abl-dmlst
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
. ./config.sh
#Import the module file that loads all necessary mods
#. "${mod_changers}/pipeline_mods"

#
# Usage ./abl_dmlst.sh path_to_list path_for_output_file [Abaum|Ecoli]
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./abl_dmlst.sh path_for_list_file [Abaum|Ecoli]"
	exit 0
fi

# Loop through and act on each sample name in the passed/provided list
while IFS= read -r var  || [ -n "$var" ]; do
	sample_name=$(echo "${var}" | cut -d'/' -f2 | tr -d '[:space:]')
	project=$(echo "${var}" | cut -d'/' -f1 | tr -d '[:space:]')
	if [[ "${3}" == "Abaum" ]]; then
		if [[ -s "${processed}/${project}/${sample_name}/MLST/${sample_name}.mlst" ]] && [[ -s "${processed}/${project}/${sample_name}/MLST/${sample_name}_abaumannii.mlst" ]]; then
			echo "${processed}/${project} has both torsten abaums"
		elif [[ -s "${processed}/${project}/${sample_name}/MLST/${sample_name}.mlst" ]]; then
			echo "${processed}/${project} has only abaumannii_2"
		elif [[ -s "${processed}/${project}/${sample_name}/MLST/${sample_name}_abaumannii.mlst" ]]; then
			echo "${processed}/${project} has only abaumannii"
		else
			echo "${processed}/${project} has no mlst files"
		fi
	elif [[ "${3}" == "Ecoli" ]]; then
		if [[ -s "${processed}/${project}/${sample_name}/MLST/${sample_name}.mlst" ]] && [[ -s "${processed}/${project}/${sample_name}/MLST/${sample_name}_ecoli_2.mlst" ]]; then
			echo "${processed}/${project} has both torsten ecolis"
		elif [[ -s "${processed}/${project}/${sample_name}/MLST/${sample_name}.mlst" ]]; then
			echo "${processed}/${project} has only ecoli"
		elif [[ -s "${processed}/${project}/${sample_name}/MLST/${sample_name}_ecoli_2.mlst" ]]; then
			echo "${processed}/${project} has only ecoli_2"
		else
			echo "${processed}/${project} has no mlst files"
		fi
	fi
done < "${1}"

echo "All isolates completed"
global_end_time=$(date "+%m-%d-%Y @ %Hh_%Mm_%Ss")
#Script exited gracefully (unless something else inside failed)
printf "%s %s" "Act_by_list_dmlst.sh has completed check of snd MLSTs" "${global_end_time}" | mail -s "act_by_list complete" nvx4@cdc.gov
exit 0
