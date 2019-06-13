#!/bin/sh -l

#$ -o check_nn.out
#$ -e check_nn.err
#$ -N check_nn
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
. ./config.sh
#Import the module file that loads all necessary mods
. "${mod_changers}/pipeline_mods"

#
# Usage ./check_node_naming_by_list.sh path_to_list ResGANNOT_identifier(YYYYMMDD)
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./check_node_naming_by_list.sh path_to_list_file ResGANNOT_Identifier(YYYYMMDD)"
	exit 0
fi

# Loop through and act on each sample name in the passed/provided list
counter=0
while IFS=read -r var  || [ -n "$var" ]; do
	sample_name=$(echo "${var}" | cut -d'/' -f2 | tr -d '[:space:]')
	project=$(echo "${var}" | cut -d'/' -f1 | tr -d '[:space:]')
	if  [[ -s ${processed}/${project}/${sample_name}/Assembly/${sample_name}_scaffolds_trimmed.fasta ]]; then
		identifier=$(head -n1 ${processed}/${project}/${sample_name}/Assembly/${sample_name}_scaffolds_trimmed.fasta | rev | cut -d'_' -f4- | rev)
		if [[ "${identifier}" = ">NODE" ]]; then
			echo "${counter} - ${project}/${sample_name}"
		elif [[ "${identifier}" = ">${sample_name}" ]]; then
			#echo "${counter}"
			:
		else
			echo "${counter} - ${project}/${sample_name} - I DONT KNOW - ${identifier}"
		fi
	else
		#echo "XXX-${counter} - ${project}/${sample_name} - NO_ASSEMBLY"
		:
	fi
	counter=$(( counter + 1 ))
done < "${1}"
echo "All isolates completed"
global_end_time=$(date "+%m-%d-%Y @ %Hh_%Mm_%Ss")
#Script exited gracefully (unless something else inside failed)
printf "%s %s" "check_node_naming_by_list.sh has completed check for AR completion " "${global_end_time}" | mail -s "node checker complete" nvx4@cdc.gov
exit 0
