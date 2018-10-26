#!/bin/sh -l

#$ -o act_by_list_barebones_4.out
#$ -e act_by_list_barebones_4.err
#$ -N ablb4
#$ -cwd
#$ -q all.q

#Import the config file with shortcuts and settings
. ./config.sh
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
if [[ ! -d "${share}/${3}" ]]; then
	mkdir "${share}/${3}"
fi

while IFS= read -r var; do
	sample_name=$(echo "${var}" | cut -d'/' -f2 | tr -d '[:space:]')
	project=$(echo "${var}" | cut -d'/' -f1 | tr -d '[:space:]')
	if [[ -f "${processed}/${project}/${sample_name}/c-sstar/${sample_name}.${resGANNOT_srst2_filename}.gapped_98_sstar_summary.txt" ]]; then
		ar_file="${processed}/${project}/${sample_name}/c-sstar/${sample_name}.${resGANNOT_srst2_filename}.gapped_98_sstar_summary.txt"
	elif [[ -f "${processed}/${project}/${sample_name}/c-sstar/${sample_name}.ResGANNOT_20180724.gapped_98_sstar_summary.txt" ]]; then
		ar_file="${processed}/${project}/${sample_name}/c-sstar/${sample_name}.ResGANNOT_20180724.gapped_98_sstar_summary.txt"
	elif [[ -f "${processed}/${project}/${sample_name}/c-sstar/${sample_name}.ResGANNOT_20180608.gapped_98_sstar_summary.txt" ]]; then
		ar_file="${processed}/${project}/${sample_name}/c-sstar/${sample_name}.ResGANNOT_20180608.gapped_98_sstar_summary.txt"
	else
		ls "${processed}/${project}/${sample_name}/c-sstar/"
	fi
	#echo "Test: ${ar_file}"
	while IFS='' read -r line || [[ -n "$line" ]]; do
		#echo "${line}}"
		gene=$(echo "${line}" | tr -s "[:blank:]" | cut -d' ' -f5)
		if [[ "${gene}" = "${2}" ]]; then
			contig=$(echo "${line}" | tr -s "[:blank:]" | cut -d' ' -f6)
			contig=$(echo "${contig/node/NODE}")
			#echo "${counter}:${project}:${sample_name}:${gene}:${contig}" >> ${share}/crab_oxa23_contigs.txt
			python "${shareScript}/Contig_Writer_Exe.py" "${processed}/${project}/${sample_name}/Assembly/${sample_name}_scaffolds_trimmed.fasta" "${contig}" "${share}/${3}/${sample_name}_${gene}.fasta"
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
