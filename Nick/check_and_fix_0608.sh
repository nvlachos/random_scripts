#!/bin/sh -l

#$ -o sp0608.out
#$ -e sp0608.err
#$ -N sp0608
#$ -cwd
#$ -q all.q

#Import the config file with shortcuts and settings
. ./config.sh
#Import the module file that loads all necessary mods
. "${mod_changers}/pipeline_mods"

#
# Usage ./check_and_fix_0608.sh path_to_list
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
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
	if [[ -s "${processed}/${project}/${sample_name}/c-sstar/${sample_name}.ResGANNOT_20180608.gapped_98_sstar_summary.txt" ]]; then
		header=$(head -n1 "${processed}/${project}/${sample_name}/c-sstar/${sample_name}.ResGANNOT_20180608.gapped_98_sstar_summary.txt")
		if [[ "${header}" = *"anti-microbial"* ]]; then
			:
		else
			echo "Fixing - ${counter} - ${project}/${sample_name}"
			mv "${processed}/${project}/${sample_name}/c-sstar/${sample_name}.ResGANNOT_20180608.gapped_98_sstar_summary.txt" "${processed}/${project}/${sample_name}/c-sstar/${sample_name}.ResGANNOT_20180608.gapped_98_sstar_summary.txt_old"
			cat "${processed}/${project}/${sample_name}/c-sstar/${sample_name}.ResGANNOT_20180608.gapped_98_sstar_summary.txt_old" | tr -s " " "\t" > "${processed}/${project}/${sample_name}/c-sstar/${sample_name}.ResGANNOT_20180608.gapped_98_sstar_summary.txt"
			rm "${processed}/${project}/${sample_name}/c-sstar/${sample_name}.ResGANNOT_20180608.gapped_98_sstar_summary.txt_old"
		fi
	fi
	counter=$(( counter + 1 ))
done < "${1}"
echo "All isolates completed"
global_end_time=$(date "+%m-%d-%Y @ %Hh_%Mm_%Ss")
#Script exited gracefully (unless something else inside failed)
printf "%s %s" "Act_by_list.sh has completed check for AR completion " "${global_end_time}" | mail -s "act_by_list complete" nvx4@cdc.gov
exit 0
