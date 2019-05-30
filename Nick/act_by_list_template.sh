#!/bin/sh -l

#$ -o abl-template.out
#$ -e abl-template.err
#$ -N abl-template
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
. ./config.sh
#Import the module file that loads all necessary mods
. "${mod_changers}/pipeline_mods"

#
# Usage ./act_by_list_template.sh path_to_list path_for_output_file
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./act_by_list_template.sh path_for_output_file"
	exit 0
elif [[ -z "${2}" ]]; then
	echo  "No output file input, exiting..."
	exit 1
fi

#sed 's/,/\//g' "${mlst_file_array[2]"

# Loop through and act on each sample name in the passed/provided list
counter=0
> ${2}
while IFS= read -r var; do
	sample_name=$(echo "${var}" | cut -d'/' -f2 | tr -d '[:space:]')
	project=$(echo "${var}" | cut -d'/' -f1 | tr -d '[:space:]')
	#sed -i 's/281,1839/281\/1839/g' "${processed}/${project}/${sample_name}/MLST/${sample_name}_abaumannii.mlst"
	if [[ -f ${processed}/${project}/${sample_name}/MLST/${sample_name}.mlst ]]; then
		mlst_line=$(head -n1 ${processed}/${project}/${sample_name}/MLST/${sample_name}.mlst)
		IFS='	' read -r -a mlst_file_array <<< "$mlst_line"
		#echo "Test original: ${mlst_file_array[2]}"
		if [[ "${mlst_file_array[2]}" == *","* ]]; then
			mlst_file_array[2]=${mlst_file_array[2]/,/\/}
		fi
		if [[ "${mlst_file_array[2]}" == *"|"* ]]; then
			mlst_file_array[2]=${mlst_file_array[2]/\|/\/}
		fi
		echo -e "${mlst_file_array[@]}\n" > ${processed}/${project}/${sample_name}/MLST/${sample_name}.mlst
	fi
	if [[ -f ${processed}/${project}/${sample_name}/MLST/${sample_name}_abaumannii.mlst ]]; then
		mlst_line=$(head -n1 ${processed}/${project}/${sample_name}/MLST/${sample_name}_abaumannii.mlst)
		IFS='	' read -r -a mlst_file_array <<< "$mlst_line"
		echo "Test abaumannii: ${mlst_file_array[2]}"
		if [[ "${mlst_file_array[2]}" == *","* ]]; then
			mlst_file_array[2]=${mlst_file_array[2]/,/\/}
		fi
		if [[ "${mlst_file_array[2]}" == *"|"* ]]; then
			mlst_file_array[2]=${mlst_file_array[2]/\|/\/}
		fi
		echo -e "${mlst_file_array[@]}\n" > ${processed}/${project}/${sample_name}/MLST/${sample_name}_abaumannii.mlst
	fi
	if [[ -f ${processed}/${project}/${sample_name}/MLST/${sample_name}_ecoli_2.mlst ]]; then
		mlst_line=$(head -n1 ${processed}/${project}/${sample_name}/MLST/${sample_name}_ecoli_2.mlst)
		IFS='	' read -r -a mlst_file_array <<< "$mlst_line"
		echo "Test ecoli: ${mlst_file_array[2]}"
		if [[ "${mlst_file_array[2]}" == *","* ]]; then
			mlst_file_array[2]=${mlst_file_array[2]/,\/}
		fi
		if [[ "${mlst_file_array[2]}" == *"|"* ]]; then
			mlst_file_array[2]=${mlst_file_array[2]/\|/\/}
		fi
		echo -e "${mlst_file_array[@]}\n" > ${processed}/${project}/${sample_name}/MLST/${sample_name}_ecoli_2.mlst
	fi

	#echo "${counter}"

	#counter=$(( counter + 1 ))
done < "${1}"

echo "All isolates completed"
global_end_time=$(date "+%m-%d-%Y @ %Hh_%Mm_%Ss")
#Script exited gracefully (unless something else inside failed)
printf "%s %s" "Act_by_list_template.sh has completed check of snd MLSTs" "${global_end_time}" | mail -s "act_by_list complete" nvx4@cdc.gov
exit 0
