#!/bin/sh -l

#$ -o abl-ARC.out
#$ -e abl-ARC.err
#$ -N abl-AR_checker
#$ -cwd
#$ -q all.q

#Import the config file with shortcuts and settings
. ./config.sh
#Import the module file that loads all necessary mods
. "${mod_changers}/pipeline_mods"

#
# Usage ./act_by_list_AR_completion_check.sh path_to_list path_for_output_file
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to act_by_list.sh, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./act_by_list_AR_completion_check.sh path_to_list_file path_for_output_file"
	exit 0
elif [[ -z "${2}" ]]; then
	echo  "No output file input, exiting..."
	exit 1
fi

# Loop through and act on each sample name in the passed/provided list
counter=0
> ${2}
while IFS= read -r var; do
	sample_name=$(echo "${var}" | cut -d'/' -f2 | tr -d '[:space:]')
	project=$(echo "${var}" | cut -d'/' -f1 | tr -d '[:space:]')
	#echo "${counter} - ${processed}/${project}/${sample_name}/MLST/"
	info1=$(tail -n1 "${processed}/${project}/${sample_name}/MLST/${sample_name}.mlst")
	info2=""
	if [[ -f "${processed}/${project}/${sample_name}/MLST/${sample_name}_abaumannii.mlst" ]]; then
		info2=$(tail -n1 "${processed}/${project}/${sample_name}/MLST/${sample_name}_abaumannii.mlst")
	elif [[ -f "${processed}/${project}/${sample_name}/MLST/${sample_name}_ecoli_2.mlst" ]]; then
		info2=$(tail -n1 "${processed}/${project}/${sample_name}/MLST/${sample_name}_ecoli_2.mlst")
	fi
	#echo "${info}"
	info_out=$(echo "${info}" | cut -d'	' -f3-)
	#echo "${info_out}"
	if [[ "${info_out}" == *","* ]] || [[ "${info_out}" == *"/"* ]]; then
	#		echo "${project}/${sample}	dual!!!" >> "${2}"
		echo "${project}/${sample_name}	${info1}"
	else
		#echo "${project}/${sample}	unknown?!" >> "${2}"
		#echo "${project}/${sample}	unknown ?!"
		:
	fi
	if [[ ! -z "${info2}" ]]; then
		info_out2=$(echo "${info2}" | cut -d'	' -f3-)
		#echo "${info_out}"
		if [[ "${info_out2}" == *","* ]] || [[ "${info_out2}" == *"/"* ]]; then
		#		echo "${project}/${sample}	dual!!!" >> "${2}"
			echo "${project}/${sample_name}	${info2}"
		else
			#echo "${project}/${sample}	unknown?!" >> "${2}"
			#echo "${project}/${sample}	unknown ?!"
			:
		fi
	fi
done < "${1}"

echo "All isolates completed"
global_end_time=$(date "+%m-%d-%Y @ %Hh_%Mm_%Ss")
#Script exited gracefully (unless something else inside failed)
printf "%s %s" "Act_by_list.sh has completed check of snd MLSTs" "${global_end_time}" | mail -s "act_by_list 2MLST complete" nvx4@cdc.gov
exit 0
