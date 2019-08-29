#!/bin/sh -l

#$ -o abl-template1.out
#$ -e abl-template1.err
#$ -N abl-template1
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
. ./config.sh
#Import the module file that loads all necessary mods
#. "${mod_changers}/pipeline_mods"

#
# Usage ./act_by_list_template.sh path_to_list path_for_output_file
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
	echo ""${processed}/${project}/${sample_name}"
	if [[ -f "${processed}/${project}/${sample_name}/GAMA/${sample_name}_ResGANNCBI_20190826.GAMA" ]]; then
		echo "1"
		mv "${processed}/${project}/${sample_name}/GAMA/${sample_name}_ResGANNCBI_20190826.GAMA" "${processed}/${project}/${sample_name}/GAMA/${sample_name}.ResGANNCBI_20190826.GAMA"
	fi
	if [[ -f "${processed}/${project}/${sample_name}/GAMA/${sample_name}_ResGANNOT_20190812.GAMA" ]]; then
		echo "2"
		mv "${processed}/${project}/${sample_name}/GAMA/${sample_name}_ResGANNOT_20190812.GAMA" "${processed}/${project}/${sample_name}/GAMA/${sample_name}.ResGANNOT_20190812.GAMA"
	fi
	if [[ -f "${processed}/${project}/${sample_name}/GAMA_plasFlow/${sample_name}_ResGANNCBI_20190826.GAMA" ]]; then
		echo "3"
		mv "${processed}/${project}/${sample_name}/GAMA_plasFlow/${sample_name}_ResGANNCBI_20190826.GAMA" "${processed}/${project}/${sample_name}/GAMA_plasFlow/${sample_name}.ResGANNCBI_20190826.GAMA"
	fi
	if [[ -f "${processed}/${project}/${sample_name}/GAMA_plasFlow/${sample_name}_ResGANNOT_20190812.GAMA" ]]; then
		echo "4"
		mv "${processed}/${project}/${sample_name}/GAMA_plasFlow/${sample_name}_ResGANNOT_20190812.GAMA" "${processed}/${project}/${sample_name}/GAMA_plasFlow/${sample_name}.ResGANNOT_20190812.GAMA"
	fi
done < "${1}"

echo "All isolates completed"
global_end_time=$(date "+%m-%d-%Y @ %Hh_%Mm_%Ss")
#Script exited gracefully (unless something else inside failed)
printf "%s %s" "Act_by_list_template.sh has completed check of snd MLSTs" "${global_end_time}" | mail -s "act_by_list complete" nvx4@cdc.gov
exit 0
