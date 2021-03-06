#!/bin/sh -l

#$ -o abl-ARC.out
#$ -e abl-ARC.err
#$ -N abl-plasCATTER
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
. ./config.sh
#Import the module file that loads all necessary mods
. "${mod_changers}/pipeline_mods"

#
# Usage ./act_by_list_concat_csstar_plasFlow_files.sh path_to_list date_of_ResGANNCBI_DB path_for_output_file
#

number='^[0-9]+$'

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./act_by_list_concat_csstar_plasFlow_files.sh path_to_list_file ResGANNCBI_Identifier(YYYYMMDD) path_for_output_file"
	exit 0
elif ! [[ ${2} =~ $number ]] || [[ -z "${2}" ]]; then
	echo "${2} is not a number or is empty. Please input correct date for ResGANNCBI DB...exiting"
	exit 2
elif [[ ! -f "${3}" ]]; then
	echo "${3} does not exist"
	outdir=$(dirname "${3}")
	if [[ ! -d "${outdir}" ]]; then
		echo "${outdir} does not exist"
		mkdir -p "${outdir}"
	fi
	>${3}
fi

# Loop through and act on each sample name in the passed/provided list
counter=0
while IFS= read -r var || [ -n "$var" ]; do
	sample_name=$(echo "${var}" | cut -d'/' -f2 | tr -d '[:space:]')
	project=$(echo "${var}" | cut -d'/' -f1 | tr -d '[:space:]')
	#echo "${counter} - ${processed}/${project}/${sample_name}/MLST/"
	if [[ -f "${processed}/${project}/${sample_name}/c-sstar_plasFlow/${sample_name}.ResGANNCBI_${2}.gapped_40_sstar_summary.txt" ]]; then
		cat "${processed}/${project}/${sample_name}/c-sstar_plasFlow/${sample_name}.ResGANNCBI_${2}.gapped_40_sstar_summary.txt" >> ${3}
	else
		echo "${project}/${sample_name} does not have c-sstar_plasFlow 40 file"
	fi
done < "${1}"

echo "All isolates completed"
global_end_time=$(date "+%m-%d-%Y @ %Hh_%Mm_%Ss")
#Script exited gracefully (unless something else inside failed)
printf "%s %s" "act_by_list_concat_csstar_plasFlow_files.sh has completed check of snd MLSTs" "${global_end_time}" | mail -s "act_by_list complete" nvx4@cdc.gov
exit 0
