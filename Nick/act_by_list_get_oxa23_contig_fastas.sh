#!/bin/sh -l

#$ -o act_by_list_barebones_4.out
#$ -e act_by_list_barebones_4.err
#$ -N ablb4
#$ -cwd
#$ -q all.q

#Import the config file with shortcuts and settings
. /scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/scripts/config.sh
#Import the module file that loads all necessary mods
. "${mod_changers}/pipeline_mods"

#
# Usage ./act_by_list.sh list_name(currently has to be placed in /scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR folder) Description of list function
#
# script changes depending on what needs to be run through the list
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to act_by_list.sh, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./act_by_list.sh path_to_list_file(single sample ID per line, e.g. B8VHY/1700128 (it must include project id also)) gene_to_find output_folder"
	echo "Output location varies depending on which tasks are performed but will be found somewhere under ${processed}"
	exit 0
elif [[ -z ${1} ]]; then
	echo "Empty list given to get_contig_with_geneX.sh...exiting"
	exit
elif [[ -z ${2} ]]; then
	echo "Empty gene given to get_contig_with_geneX.sh...exiting"
	exit
elif [[ -z ${3} ]]; then
	echo "Empty output folder given to get_contig_with_geneX.sh...exiting"
	exit
fi

# Loop through and act on each sample name in the passed/provided list
counter=1
mkdir "${share}/${3}"

while IFS= read -r var; do
	sample_name=$(echo "${var}" | cut -d'/' -f2 | tr -d '[:space:]')
	project=$(echo "${var}" | cut -d'/' -f1 | tr -d '[:space:]')
	if [[ -f "${processed}/${project}/${sample_name}/c-sstar/${sample_name}.cdiff_gyrs_srst2.gapped_98_sstar_summary.txt" ]]; then
		ar_file="${processed}/${project}/${sample_name}/c-sstar/${sample_name}.cdiff_gyrs_srst2.gapped_98_sstar_summary.txt"
		blast_file="${processed}/${project}/${sample_name}/c-sstar/cdiff_gyrs_srst2_gapped/${sample_name}_scaffolds_trimmed.blastn.tsv"
	elif [[ -f "${processed}/${project}/${sample_name}/c-sstar/${sample_name}.${resGANNOT_srst2_filename}.gapped_98_sstar_summary.txt" ]]; then
		ar_file="${processed}/${project}/${sample_name}/c-sstar/${sample_name}.${resGANNOT_srst2_filename}.gapped_98_sstar_summary.txt"
		blast_file="${processed}/${project}/${sample_name}/c-sstar/${resGANNOT_srst2_filename}_gapped/${sample_name}_scaffolds_trimmed.blastn.tsv"
	elif [[ -f "${processed}/${project}/${sample_name}/c-sstar/${sample_name}.ResGANNOT_20180724.gapped_98_sstar_summary.txt" ]]; then
		ar_file="${processed}/${project}/${sample_name}/c-sstar/${sample_name}.ResGANNOT_20180724.gapped_98_sstar_summary.txt"
		blast_file="${processed}/${project}/${sample_name}/c-sstar/ResGANNOT_20180724_gapped/${sample_name}_scaffolds_trimmed.blastn.tsv"
	elif [[ -f "${processed}/${project}/${sample_name}/c-sstar/${sample_name}.ResGANNOT_20180608.gapped_98_sstar_summary.txt" ]]; then
		ar_file="${processed}/${project}/${sample_name}/c-sstar/${sample_name}.ResGANNOT_20180608.gapped_98_sstar_summary.txt"
		blast_file="${processed}/${project}/${sample_name}/c-sstar/ResGANNOT_20180608_gapped/${sample_name}_scaffolds_trimmed.blastn.tsv"
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
			echo -e ">${counter}_${project}_${sample_name}_${contig}_${gene}_${allele}\n${seqmatch}" >> ${share}/${2}.fasta
			#python "${shareScript}/Contig_Writer_Exe.py" "${processed}/${project}/${sample_name}/Assembly/${sample_name}_scaffolds_trimmed.fasta" "${contig}" "${share}/oxa_23_fastas/${sample_name}_oxa_${oxa}.fasta"
			#counter=$(( counter + 1 ))
		fi
		#counter=$(( counter + 1 ))
	done < "${ar_file}"
	counter=$(( counter + 1 ))
done < "${share}/${1}"
echo "All isolates completed"
global_end_time=$(date "+%m-%d-%Y @ %Hh_%Mm_%Ss")
#Script exited gracefully (unless something else inside failed)
printf "%s %s" "Act_by_list.sh has completed ${2} " "${global_end_time}" | mail -s "act_by_list complete" nvx4@cdc.gov
exit 0
