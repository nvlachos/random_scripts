#!/bin/sh -l

#$ -o blast_fasta.out
#$ -e blast_fasta.err
#$ -N blast_fasta
#$ -cwd
#$ -q all.q

#Import the config file with shortcuts and settings
. ./config.sh
#Import the module file that loads all necessary mods
. "${mod_changers}/pipeline_mods"

#
# Usage ./blast_fasta_to_directory_of_assemblies.sh list_name(currently has to be placed in /scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR folder) fasta_to_blast directory_of_fastas_to_blast_against
#
# script changes depending on what needs to be run through the list
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./blast_fasta_to_directory_of_assemblies.sh path_to_list path_to_fasta_file path_to_directory_of_fastas"
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
	if [[ "${assembly}" == *".fasta" ]] || [[ "${assembly}" == *".fna" ]]; then
			echo "Attempting to blast ${assembly}"
			blastOut="${assembly}.blast"
			blastn -query ${assembly} -db ${2} -out $blastOut -word_size 10 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen"
			sort -k4,4r -k3,3nr -o ${blastOut} ${blastOut}
			#echo "Truly?"
			info=$(head -n1 ${blastOut})
			echo "${info}"
			echo "${info}" >> ${bestlist}
			#echo "After..."
			echo "completed ${assembly} in loop at ${current_time}"
	else
		echo "${assembly} is not a fasta file"
	fi
done

echo "All runs completed"
global_end_time=$(date "+%m-%d-%Y @ %Hh_%Mm_%Ss")
#Script exited gracefully (unless something else inside failed)
printf "%s %s" "blast_fasta_to_directory_of_assemblies.sh has completed whatever it was doing at" "${global_end_time}" | mail -s "blast_fasta_to_directory_of_assemblies complete" nvx4@cdc.gov
exit 0
