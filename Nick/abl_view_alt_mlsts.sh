#!/bin/sh -l

#$ -o abl-vaMLST.out
#$ -e abl-vaMLST.err
#$ -N abl-vaMLST
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
. ./config.sh
#Import the module file that loads all necessary mods
. "${mod_changers}/pipeline_mods"

#
# Usage ./abl_view_alt_mlsts.sh path_to_list path_for_output_file
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./abl_view_alt_mlsts.sh path_to_list_file path_for_output_file"
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
	info1=$(head -n1 "${processed}/${project}/${sample_name}/MLST/${sample_name}.mlst")
	if [[ "${info1}" == *"}" ]]; then
		echo "Cutting trailing }"
		echo "Before - ${info1}"
		info1=$(echo "${info1}" | rev | cut -d'}' -f2 | rev)
		echo "After - ${info1}"
	fi
	info2=""
	if [[ -f "${processed}/${project}/${sample_name}/MLST/${sample_name}_abaumannii.mlst" ]]; then
		info2=$(head -n1 "${processed}/${project}/${sample_name}/MLST/${sample_name}_abaumannii.mlst")
	elif [[ -f "${processed}/${project}/${sample_name}/MLST/${sample_name}_ecoli_2.mlst" ]]; then
		info2=$(head -n1 "${processed}/${project}/${sample_name}/MLST/${sample_name}_ecoli_2.mlst")
	fi
	echo "${info1}"
	info_out=$(echo "${info1}" | cut -d'	' -f3-)
	type=$(echo ${info1} | cut -d' ' -f3)
	#echo "${info_out}"
	if [[ "${info_out}" == *","* ]] || [[ "${info_out}" == *"/"* ]] || [[ "${info_out}" == *"|"* ]]; then
		echo "1-DUAL - ${type} - ${info1}" >> "${2}"
		echo "1-${project}/${sample_name} ${type}	${info1}"
	elif [[ "${type}" == "-" ]]; then
		echo "1-UNKN - ${type} - ${info1}" >> "${2}"
		echo "1-${project}/${sample_name}	${type} ${info1}"
	else
		echo "1-Normal? - ${info1}" >> "${2}"
		echo "1-${project}/${sample} ${tpe} ${info1}"
		:
	fi
	if [[ ! -z "${info2}" ]]; then
		type2=$(echo ${info2} | cut -d'	' -f3)
		info_out2=$(echo "${info2}" | cut -d' ' -f3-)
		#echo "${info_out}"
		if [[ "${info_out2}" == *","* ]] || [[ "${info_out2}" == *"/"* ]] || [[ "${info_out2}" == *"|"* ]]; then
			echo 2-"DUAL - ${type2} - ${info2}" >> "${2}"
			echo "2-${project}/${sample_name}	${type2} ${info2}"
		elif [[ "${type2}" == "-" ]]; then
			echo "2-UNKN - ${type2} - ${info2}" >> "${2}"
			echo "2-${project}/${sample_name}	${type2} ${info2}"
		else
			echo "2-Normal? - ${type2} - ${info2}" >> "${2}"
			echo "2-${project}/${sample} ${type2} ${info2}"
			:
		fi
	fi
done < "${1}"

echo "All isolates completed"
global_end_time=$(date "+%m-%d-%Y @ %Hh_%Mm_%Ss")
#Script exited gracefully (unless something else inside failed)
printf "%s %s" "abl_view_alt_mlsts.sh has completed check of snd MLSTs" "${global_end_time}" | mail -s "act_by_list complete" nvx4@cdc.gov
exit 0
