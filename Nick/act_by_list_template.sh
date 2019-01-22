#!/bin/sh -l

#$ -o act_by_list_barebones_4.out
#$ -e act_by_list_barebones_4.err
#$ -N abl-prokka
#$ -cwd
#$ -q all.q

#Import the config file with shortcuts and settings
. ./config.sh
#Import the module file that loads all necessary mods
. "${mod_changers}/pipeline_mods"

#
# Usage ./act_by_list_AR_completion_check.sh path_to_list ResGANNOT_identifier(YYYYMMDD)
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to act_by_list.sh, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./act_by_list_AR_completion_check.sh path_to_list_file ResGANNOT_Identifier(YYYYMMDD)"
	exit 0
fi

# Loop through and act on each sample name in the passed/provided list
counter=0
while IFS= read -r var; do
	sample_name=$(echo "${var}" | cut -d'/' -f2 | tr -d '[:space:]')
	project=$(echo "${var}" | cut -d'/' -f1 | tr -d '[:space:]')
	if [[ -f "${processed}/${project}/${sample_name}/c-sstar/${sample_name}.ResGANNOT_20181204.gapped_80_sstar_summary.txt" ]]; then
		while IFS= read -r var; do
				gene=$(echo "${var}" | cut -d'	' -f3)
				id=$(echo "${var}" | cut -d'	' -f7)
				length=$(echo "${var}" | cut -d'	' -f10)
				break
		done < "${processed}/${project}/${sample_name}/c-sstar/${sample_name}.ResGANNOT_20181204.gapped_80_sstar_summary.txt"
		if [[ "${gene}" = "blasst" ]]; then
			echo "${project}/${sample_name}	[${id}/${length}]"
		fi
	fi













	# if [[ -f "${processed}/${project}/${sample_name}/MLST/${sample_name}_abaumannii.mlst" ]]; then
	# 	ST_type=$(tail -n1 "${processed}/${project}/${sample_name}/MLST/${sample_name}_abaumannii.mlst" | cut -d'	' -f3)
	# 	echo "${ST_type}"
	# 	if [[ "${ST_type}" = "208"* ]]; then
	# 		echo "${project}/${sample_name}"
	# 		echo "${ST_type}	${project}/${sample_name}" >> "${2}"
	# 	fi
	# fi
done < "${1}"
echo "All isolates completed"
global_end_time=$(date "+%m-%d-%Y @ %Hh_%Mm_%Ss")
#Script exited gracefully (unless something else inside failed)
printf "%s %s" "Act_by_list.sh has completed check for AR completion " "${global_end_time}" | mail -s "act_by_list complete" nvx4@cdc.gov
exit 0
