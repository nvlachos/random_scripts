#!/bin/sh -l

#$ -o redact_OA.out
#$ -e redact_OA.err
#$ -N redact_OA
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
. ./config.sh
#Import the module file that loads all necessary mods
. "${mod_changers}/pipeline_mods"

#
# Usage ./redact_OA.sh project_ID analysis_ID path_to_crosswalk_file(separated by colons)
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./redact_OA.sh path_for_output_file"
	exit 0
elif [[ -z "${2}" ]]; then
	echo  "No analysis ID input, exiting..."
	exit 1
elif [[ -z "${3}" ]] || [[ ! -f ${3} ]]; then
	echo  "No crosswalk file input, exiting..."
	exit 1
fi

#sed 's/,/\//g' "${mlst_file_array[2]"

# Loop through and act on each sample name in the passed/provided list

if [[ -d /scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/Phylogeny_analyses/${1} ]]; then
	if  [[ -d /scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/Phylogeny_analyses/${1}/${2} ]]; then
		echo "Redacting Phylogeny folder: /scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/Phylogeny_analyses/${1}"
		cp /scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/Phylogeny_analyses/${1}/output/snvMatrix.tsv snvMatrix_redacted.tsv
		cp /scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/Phylogeny_analyses/${1}/output/phylogeneticTree.newick phylogeneticTree_redacted.newick
		while IFS= read -r var; do
			original_name=$(echo "${var}" | cut -d':' -f1) | tr -d '[:space:]')
			redacted_name=$(echo "${var}" | cut -d':' -f2) | tr -d '[:space:]')
			sed -i "s/${original_name}/${redacted_name}/g" /scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/Phylogeny_analyses/${1}/output/phylogeneticTree_redacted.newick
			sed -i "s/${original_name}/${redacted_name}/g" /scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/Phylogeny_analyses/${1}/output/snvMatrix_redacted.tsv
		done < ${3}
	else
		echo "Phylo: ${1} exists, but ${2} is missing"
	fi
else
	echo "Phylo: /scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/Phylogeny_analyses/${1} does not exist"
fi
	if [[ -d /scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/Projects/${1} ]]; then
		if [[ -d /scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/Projects/${1}/${2} ]]; then
			echo "Redacting Project folder: /scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/Phylogeny/${1}"
			cp /scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/Projects/${1}/${2}/${2}.nwk ${2}_redacted.nwk
			cp /scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/Phylogeny_analyses/${1}/${2}/${2}_AR_plasmid_report.csv ${2}_AR_plasmid_report_redacted.csv
			while IFS= read -r var; do
				original_name=$(echo "${var}" | cut -d':' -f1) | tr -d '[:space:]')
				redacted_name=$(echo "${var}" | cut -d':' -f2) | tr -d '[:space:]')
				sed -i "s/${original_name}/${redacted_name}/g" /scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/Projects/${1}/${2}/${2}_redacted.nwk
				sed -i "s/${original_name}/${redacted_name}/g" cp /scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/Phylogeny_analyses/${1}/${2}/${2}_AR_plasmid_report_redacted.csv
			done < ${3}
		else
			echo "Projects: ${1} exists, but ${2} is missing"
		fi
	else
		echo "Projects: /scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/Projects/${1} does not exist"
	fi
done < "${1}"

echo "Redaction complete"
global_end_time=$(date "+%m-%d-%Y @ %Hh_%Mm_%Ss")
#Script exited gracefully (unless something else inside failed)
printf "%s %s" "redact_OA.sh has completed " "${global_end_time}" | mail -s "redact_OA complete" nvx4@cdc.gov
exit 0
