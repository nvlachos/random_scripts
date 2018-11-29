#!/bin/sh -l

#$ -o act_by_list_barebones.out
#$ -e act_by_list_barebones.err
#$ -N ablb
#$ -cwd
#$ -q all.q

#Import the config file with shortcuts and settings
. ./config.sh
#Import the module file that loads all necessary mods
. "${mod_changers}/pipeline_mods"

#
# Usage ./act_by_list.sh list_name(currently has to be placed in /scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR folder)
#
# script changes depending on what needs to be run through the list
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to act_by_list.sh, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./rename_by_list.sh path_to_list_file(single sample ID per line, e.g. B8VHYproject_ID/old_name:new_name)"
		exit 0
fi

# Loop through and act on each sample name in the passed/provided list
while IFS= read -r var; do
	sample_name=$(echo "${var}" | cut -d'/' -f2 | cut -d':' -f1 | tr -d '[:space:]')
	project=$(echo "${var}" | awk -F/ '{ print $1}' | tr -d '[:space:]')
	new_sample_name=$(echo ${var} | cut -d':' -f2)
	echo "old-${sample_name}:new-${new_sample_name}"
	OUTDATADIR="${processed}/${project}/${sample_name}"
	echo "doing normal for ${sample_name}"
	"${shareScript}/rename_sample.sh" "${sample_name}" "${new_sample_name}" "${project}"
done < "${share}/${1}"
echo "All isolates completed"
global_end_time=$(date "+%m-%d-%Y @ %Hh_%Mm_%Ss")
#Script exited gracefully (unless something else inside failed)
printf "%s %s" "Act_by_list.sh has completed whatever it was doing at" "${global_end_time}" | mail -s "act_by_list complete" nvx4@cdc.gov
exit 0
