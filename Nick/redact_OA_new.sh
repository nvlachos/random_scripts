#!/bin/sh -l

#$ -o redact_OA.out
#$ -e redact_OA.err
#$ -N redact_OA
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
. ./config.sh
#Import the module file that loads all necessary mods
# . "${mod_changers}/pipeline_mods"

#
# Usage ./redact_OA.sh project_ID analysis_ID path_to_crosswalk_file(separated by colons)
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./redact_OA.sh project_ID analysis_ID path_to_crosswalk_file"
	exit 0
elif [[ -z "${2}" ]]; then
	echo  "No analysis ID input, exiting..."
	exit 1
elif [[ -z "${3}" ]]; then
	echo  "No crosswalk file input, exiting..."
	exit 1
elif [[ ! -f ${3} ]]; then
	echo "${3} does not exist"
	exit 1
fi

#sed 's/,/\//g' "${mlst_file_array[2]"

# Loop through and act on each sample name in the passed/provided list

if [[ -d ${Phyl_OA}/${1} ]]; then
	if  [[ -d ${Phyl_OA}/${1}/${2} ]]; then
		echo "Redacting Phylogeny files: ${Phyl_OA}/${1}/${2}"
		if [[ -f ${Phyl_OA}/${1}/${2}/${2}_snvMatrix_redacted.tsv ]]; then
			rm ${Phyl_OA}/${1}/${2}/${2}_snvMatrix_redacted.tsv
		fi
		cp ${Phyl_OA}/${1}/${2}/${2}_snvMatrix.tsv ${Phyl_OA}/${1}/${2}/${2}_snvMatrix_redacted.tsv
		if [[ -f ${Phyl_OA}/${1}/${2}/${2}_SNVPhyl_redacted.newick ]]; then
			rm ${Phyl_OA}/${1}/${2}/${2}_SNVPhyl.newick ${Phyl_OA}/${1}/${2}/${2}_SNVPhyl_redacted.newick
		fi
		cp ${Phyl_OA}/${1}/${2}/${2}_SNVPhyl.newick ${Phyl_OA}/${1}/${2}/${2}_SNVPhyl_redacted.newick
		while IFS= read -r var  || [ -n "$var" ]; do
			original_name=$(echo "${var}" | cut -d',' -f1 | cut -d'/' -f2 | tr -d '[:space:]')
			redacted_name=$(echo "${var}" | cut -d',' -f2 | tr -d '[:space:]')
			sed -i "s/${original_name}/${redacted_name}/g" ${Phyl_OA}/${1}/${2}/${2}_SNVPhyl_redacted.newick
			sed -i "s/${original_name}/${redacted_name}/g" ${Phyl_OA}/${1}/${2}/${2}_snvMatrix_redacted.tsv
		done < ${3}
	else
		echo "Phylo: ${1} exists, but ${2} is missing"
	fi
	echo "Redacting OA files: ${Phyl_OA}/${1}"
	if [[ -f "${Phyl_OA}/${1}/${1}_AR_plasmid_report.csv" ]]; then
		if [[ -f "${Phyl_OA}/${1}/${1}_AR_plasmid_report_redacted.csv" ]]; then
			rm "${Phyl_OA}/${1}/${1}_AR_plasmid_report_redacted.csv"
		fi
		cp "${Phyl_OA}/${1}/${1}_AR_plasmid_report.csv" "${Phyl_OA}/${1}/${1}_AR_plasmid_report_redacted.csv"
		if [[ -f "${Phyl_OA}/${1}/${1}_MASH_redacted.newick" ]]; then
			rm "${Phyl_OA}/${1}/${1}_MASH_redacted.newick"
		fi
		cp "${Phyl_OA}/${1}/${1}_MASH.newick" "${Phyl_OA}/${1}/${1}_MASH_redacted.newick"
		while IFS= read -r var  || [ -n "$var" ]; do
			original_name=$(echo "${var}" | cut -d',' -f1 | cut -d'/' -f2 | tr -d '[:space:]')
			original_project=$(echo "${var}" | cut -d'/' -f1 | tr -d '[:space:]')
			redacted_name=$(echo "${var}" | cut -d',' -f2 | tr -d '[:space:]')
			sed -i "s/${original_name}/${redacted_name}/g" ${Phyl_OA}/${1}/${1}_redacted.newick
			sed -i "s/${original_name}/${redacted_name}/g" ${Phyl_OA}/${1}/${1}_AR_plasmid_report_redacted.csv
			sed -i "s/${original_project}/NA/g" ${Phyl_OA}/${1}/${1}_AR_plasmid_report_redacted.csv
		done < ${3}
	else
		"OA: ${Phyl_OA}/${1}/${1}_AR_plasmid_report.csv does not exist, cannot redact OA"
	fi
else
	echo "${Phyl_OA}/${1} does not exist"
fi

echo "Redaction complete"
global_end_time=$(date "+%m-%d-%Y @ %Hh_%Mm_%Ss")
#Script exited gracefully (unless something else inside failed)
printf "%s %s" "redact_OA.sh has completed " "${global_end_time}" | mail -s "redact_OA complete" nvx4@cdc.gov
exit 0
