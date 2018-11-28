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
	echo "Usage is ./blast_fasta_to_directory_of_assemblies.sh input_folder absolute_path_to_fasta"
	echo "Output location varies depending on which tasks are performed but will be found somewhere under ${processed}"
	exit 0
fi

# Loop through and act on each sample name in the passed/provided list

makeblastdb -in ${2} -out ${2} -dbtype nucl

if [[ ! -d ${1} ]]; then
	echo "Directory of fastas (${1}) don't exist... exiting"
fi

bestlist="${1}/best_hits.tsv"
> ${bestlist}

for assembly in ${1}/*;
do
	if [[ "${assembly}" == *".fasta" ]]; then
			echo "Attempting to blast ${assembly}"
			blastOut=${assembly}.blast
			blastn -query ${assembly} -db ${2} -out $blastOut -word_size 10 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen"
			sort -k4,4r -k3,3nr -o ${blastOut} ${blastOut}
			$(head -n1 ${blastOut}) >> $bestlist
			echo "completed ${assembly} in loop at ${current_time}"
	else
		echo "${assembly} is not a fasta file"
	fi
done

echo "All runs completed"
global_end_time=$(date "+%m-%d-%Y @ %Hh_%Mm_%Ss")
#Script exited gracefully (unless something else inside failed)
printf "%s %s" "Act_by_list.sh has completed whatever it was doing at" "${global_end_time}" | mail -s "act_by_list complete" nvx4@cdc.gov
exit 0
