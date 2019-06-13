#!/bin/sh -l

#$ -o abl-ggbm.out
#$ -e abl-ggbm.err
#$ -N abl-ggbm
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
. ./config.sh
#Import the module file that loads all necessary mods
. "${mod_changers}/pipeline_mods"

#
# Takes a list of isolates and pulls out the sequences of the requested gene from the csstar blast output
#
# Usage ./act_by_list_get_gene_sequence_from_blast_matches.sh path_to_list gene_to_find Output_directory
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./act_by_list_get_gene_sequence_from_blast_matches.sh path_to_list_file(single sample ID per line, e.g. B8VHY/1700128 (it must include project id also)) gene_to_find output_folder"
	echo "Output location varies depending on which tasks are performed but will be found somewhere under ${processed}"
	exit 0
elif [[ -z ${1} ]]; then
	echo "Empty list given to act_by_list_get_gene_sequence_from_blast_matches...exiting"
	exit
elif [[ -z ${2} ]]; then
	echo "Empty gene given to act_by_list_get_gene_sequence_from_blast_matches...exiting"
	exit
elif [[ -z ${3} ]]; then
	echo "Empty output folder given to act_by_list_get_gene_sequence_from_blast_matches...exiting"
	exit
fi

# Loop through and act on each sample name in the passed/provided list
counter=1
mkdir "${3}"

while IFS= read -r var || [ -n "$var" ]; do
	sample_name=$(echo "${var}" | cut -d'/' -f2 | tr -d '[:space:]')
	project=$(echo "${var}" | cut -d'/' -f1 | tr -d '[:space:]')
	if [[ -f "${processed}/${project}/${sample_name}/c-sstar/${sample_name}.${resGANNOT_srst2_filename}.gapped_98_sstar_summary.txt" ]]; then
		ar_file="${processed}/${project}/${sample_name}/c-sstar/${sample_name}.${resGANNOT_srst2_filename}.gapped_98_sstar_summary.txt"
		blast_file="${processed}/${project}/${sample_name}/c-sstar/${resGANNOT_srst2_filename}_gapped/${sample_name}_scaffolds_trimmed.blastn.tsv"
	else
		ls "${processed}/${project}/${sample_name}/c-sstar/"
		break
	fi
	#echo "Test: ${ar_file}"
	while IFS='' read -r line || [[ -n "$line" ]]; do
		#echo "${line}"
		gene=$(echo "${line}" | cut -d'	' -f5)
		#echo "${gene}"
		if [[ "${gene}" = "${2}" ]]; then
			allele=$(echo "${line}" | cut -d'	' -f2)
			seqmatch="NA"
			contig="NA"
			while IFS='' read -r line || [[ -n "$line" ]]; do
				blast_allele=$(echo "${line}" | cut -d'	' -f1)
				blast_allele=${blast_allele,,}
				#echo ${allele},${blast_allele}
				if [[ "${blast_allele}" = *"${allele}"* ]]; then
					seqmatch=$(echo "${line}" | cut -d'	' -f15)
					contig=$(echo "${line}" | cut -d'	' -f2)
					#echo "Match-${seqmatch}"
					break
				fi
			done < ${blast_file}
			echo -e ">${counter}_${project}_${sample_name}_${contig}_${gene}_${allele}\n${seqmatch}" >> ${3}/${2}.fasta
			#counter=$(( counter + 1 ))
		fi
		#counter=$(( counter + 1 ))
	done < "${ar_file}"
	counter=$(( counter + 1 ))
done < "${1}"
echo "All isolates completed"
global_end_time=$(date "+%m-%d-%Y @ %Hh_%Mm_%Ss")
#Script exited gracefully (unless something else inside failed)
printf "%s %s" "act_by_list_get_gene_sequence_from_blast_matches.sh has completed ${2} " "${global_end_time}" | mail -s "act_by_list complete" nvx4@cdc.gov
exit 0
