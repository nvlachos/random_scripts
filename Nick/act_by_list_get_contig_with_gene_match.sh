#!/bin/sh -l

#$ -o abl_scon.out
#$ -e abl_scon.err
#$ -N scon
#$ -cwd
#$ -q all.q

#Import the config file with shortcuts and settings
. ./config.sh
#Import the module file that loads all necessary mods
. "${mod_changers}/pipeline_mods"

#
# Usage ./act_by_list_get_contig_with_gene_match.sh path_tolist gene_name_to_find Output_directory
#
# Finds any isolate in the list that has the matching gene name and writes the matching contig to a file
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to./act_by_list_get_contig_with_gene_match.sh, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./act_by_list.sh path_to_list_file(single sample ID per line, e.g. B8VHY/1700128 (it must include project id also)) gene_to_find output_folder"
	exit 0
elif [[ -z ${1} ]]; then
	echo "Empty list given to ./act_by_list_get_contig_with_gene_match.sh...exiting"
	exit
elif [[ -z ${2} ]]; then
	echo "Empty gene given to ./act_by_list_get_contig_with_gene_match.sh...exiting"
	exit
elif [[ -z ${3} ]]; then
	echo "Empty output folder given to ./act_by_list_get_contig_with_gene_match.sh...exiting"
	exit
fi

# create output directory
counter=1
if [[ ! -d "${3}" ]]; then
	mkdir -p "${3}"
fi

need_qsub="false"
> "${3}/csstar_to_do.txt"

# Loop through and act on each sample name in the passed/provided list
while IFS= read -r var; do
	sample_name=$(echo "${var}" | cut -d'/' -f2 | tr -d '[:space:]')
	project=$(echo "${var}" | cut -d'/' -f1 | tr -d '[:space:]')
	if [[ -f "${processed}/${project}/${sample_name}/c-sstar/${sample_name}.${resGANNOT_srst2_filename}.gapped_98_sstar_summary.txt" ]]; then
		ar_file="${processed}/${project}/${sample_name}/c-sstar/${sample_name}.${resGANNOT_srst2_filename}.gapped_98_sstar_summary.txt"
	# elif [[ -f "${processed}/${project}/${sample_name}/c-sstar/${sample_name}.ResGANNOT_20181003.gapped_98_sstar_summary.txt" ]]; then
	# 	ar_file="${processed}/${project}/${sample_name}/c-sstar/${sample_name}.ResGANNOT_20181003.gapped_98_sstar_summary.txt"
	# elif [[ -f "${processed}/${project}/${sample_name}/c-sstar/${sample_name}.ResGANNOT_20180817.gapped_98_sstar_summary.txt" ]]; then
	# 	ar_file="${processed}/${project}/${sample_name}/c-sstar/${sample_name}.ResGANNOT_20180817.gapped_98_sstar_summary.txt"
	# elif [[ -f "${processed}/${project}/${sample_name}/c-sstar/${sample_name}.ResGANNOT_20180724.gapped_98_sstar_summary.txt" ]]; then
	# 	ar_file="${processed}/${project}/${sample_name}/c-sstar/${sample_name}.ResGANNOT_20180724.gapped_98_sstar_summary.txt"
	# elif [[ -f "${processed}/${project}/${sample_name}/c-sstar/${sample_name}.ResGANNOT_20180608.gapped_98_sstar_summary.txt" ]]; then
	# 	ar_file="${processed}/${project}/${sample_name}/c-sstar/${sample_name}.ResGANNOT_20180608.gapped_98_sstar_summary.txt"
	else
		echo "${project}/${sample_name}" >> "${3}/csstar_to_do.txt"
		need_qsub="true"
	fi
	#echo "Test: ${ar_file}"
	if [[ "${need_qsub}" = "true" ]]; then
		while IFS='' read -r line || [[ -n "$line" ]]; do
			#echo "${line}}"
			gene=$(echo "${line}" | tr -s "[:blank:]" | cut -d' ' -f5)
			if [[ "${gene}" = "${2}" ]]; then
				contig=$(echo "${line}" | tr -s "[:blank:]" | cut -d' ' -f6)
				contig=$(echo "${contig/node/NODE}")
				python "${shareScript}/Contig_Writer_Exe.py" "${processed}/${project}/${sample_name}/Assembly/${sample_name}_scaffolds_trimmed.fasta" "${contig}" "${3}/${sample_name}_${gene}.fasta"
				#counter=$(( counter + 1 ))
			fi
			#counter=$(( counter + 1 ))
		done < "${ar_file}"
	fi
	counter=$(( counter + 1 ))
done < "${1}"
if [[ "${need_qsub}" = "true" ]]; then
	echo "Some samples needed to be run with the newest ResGANNOT DB, submitting to qsub and then rerunning"
	${shareScript}/abl_mass_qsub_csstar.sh "${3}/csstar_to_do.txt" 25
	${shareScript}/"act_by_list_get_contig_with_gene_match.sh" ${1} ${2} ${3}
fi
echo "All isolates completed"
global_end_time=$(date "+%m-%d-%Y @ %Hh_%Mm_%Ss")
#Script exited gracefully (unless something else inside failed)
printf "%s %s" "_get_contig_with_gene_match.sh has completed ${2} " "${global_end_time}" | mail -s "act_by_list complete" nvx4@cdc.gov
exit 0
