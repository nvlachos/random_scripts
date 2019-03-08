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
	if [[ -f ${processed}/${project}/${sample_name}/srst2/${sample_name}_ResGANNOT_20180608__genes__ResGANNOT_20180608_srst2__results.txt ]]; then
		mv ${processed}/${project}/${sample_name}/srst2/${sample_name}_ResGANNOT_20180608__genes__ResGANNOT_20180608_srst2__results.txt ${processed}/${project}/${sample_name}/srst2/${sample_name}__genes__ResGANNOT_20180608_srst2__results.txt
	fi
	if [[ -f ${processed}/${project}/${sample_name}/srst2/${sample_name}_ResGANNOT_20180608__fullgenes__ResGANNOT_20180608_srst2__results.txt ]]; then
		mv ${processed}/${project}/${sample_name}/srst2/${sample_name}_ResGANNOT_20180608__fullgenes__ResGANNOT_20180608_srst2__results.txt ${processed}/${project}/${sample_name}/srst2/${sample_name}__genes__ResGANNOT_20180608_srst2__results.txt
	fi
	if [[ -f ${processed}/${project}/${sample_name}/srst2/${sample_name}_S1_L001_R1_001.fastq.gz ]]; then
		rm ${processed}/${project}/${sample_name}/srst2/${sample_name}_S1_L001_R1_001.fastq.gz
	fi
	if [[ -f ${processed}/${project}/${sample_name}/srst2/${sample_name}_S1_L001_R2_001.fastq.gz ]]; then
		rm ${processed}/${project}/${sample_name}/srst2/${sample_name}_S1_L001_R2_001.fastq.gz
	fi
done < "${1}"
echo "All isolates completed"
global_end_time=$(date "+%m-%d-%Y @ %Hh_%Mm_%Ss")
#Script exited gracefully (unless something else inside failed)
printf "%s %s" "Act_by_list.sh has completed check for AR completion " "${global_end_time}" | mail -s "act_by_list complete" nvx4@cdc.gov
exit 0
