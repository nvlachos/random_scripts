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
	echo "No argument supplied to blast_list.sh, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./blast_list_to_reference.sh path_to_list_file(single sample ID per line, e.g. B8VHY/1700128 (it must include project id also)) absolute_output_folder absolute_path_to_reference"
	echo "Output location varies depending on which tasks are performed but will be found somewhere under ${processed}"
	exit 0
fi

# Loop through and act on each sample name in the passed/provided list

makeblastdb -in ${3} -out ${3} -dbtype nucl

if [[ ! -d ${2} ]]; then
	mkdir ${2}
fi

while IFS= read -r var; do
	sample_name=$(echo "${var}" | awk -F/ '{ print $2}' | tr -d '[:space:]')
	project=$(echo "${var}" | awk -F/ '{ print $1}' | tr -d '[:space:]')
	SOURCEDATADIR="${processed}/${project}/${sample_name}"

	blastOut=${2}/${sample_name}.blast
	blastn -query ${SOURCEDATADIR}/Assembly/${sample_name}_scaffolds_trimmed.fasta -db ${3} -out $blastOut -word_size 10 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen"

	echo "completed ${var} in loop at ${current_time}"
done < "${1}"
echo "All runs completed"
global_end_time=$(date "+%m-%d-%Y @ %Hh_%Mm_%Ss")
#Script exited gracefully (unless something else inside failed)
printf "%s %s" "Act_by_list.sh has completed whatever it was doing at" "${global_end_time}" | mail -s "act_by_list complete" nvx4@cdc.gov
exit 0
