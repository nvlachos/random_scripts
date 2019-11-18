#!/bin/sh -l

#$ -o run_sum.out
#$ -e run_sum.err
#$ -N run_sum
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi
. ./config.sh

#
# Description: Creates a summary file for the run and prints out a one word status and a short description of each step being reported
#
# Usage ./run_sum.sh input_list output_filename
#
# Output loction: parameter
#
# Modules required: None
#
# v1.0 (11/18/2019)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
elif [[ -z "${1}" ]] || [[ ! -f "${1}" ]]; then
	echo "Empty list supplied to run_sum.sh, exiting"
	exit 1
elif [[ "${1}" = "-h" ]]; then
	echo "Usage is ./run_sum.sh miseq_run_ID -vo(optional)"
	echo "Output is saved to ${processed}/miseq_run_ID"
	exit 0
fi

echo "Creating run summary at ${runsumdate}"

# Run validate_piperun.sh on every sample in the list and cat output into one summary file
while IFS= read -r samples || [ -n "$samples" ]; do
	file=$(echo "${samples}" | awk -F/ '{ print $2}' | tr -d '[:space:]')
	proj=$(echo "${samples}" | awk -F/ '{ print $1}' | tr -d '[:space:]')
	"${shareScript}/validate_piperun_aniDB2.sh" "${file}" "${proj}" > "${processed}/${proj}/${file}/${file}_pipeline_stats_aniDB2.txt"
	cat "${processed}/${proj}/${file}/${file}_pipeline_stats_aniDB2.txt" >> "${2}"
done < ${1}
