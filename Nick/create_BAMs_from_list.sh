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
	echo "No argument supplied to act_by_list.sh, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./act_by_list.sh path_to_list_file(single sample ID per line, e.g. B8VHY/1700128 (it must include project id also)) output_directory path_to_reference_sequence"
	echo "Output location varies depending on which tasks are performed but will be found somewhere under ${processed}"
	exit 0
fi

# Loop through and act on each sample name in the passed/provided list

ref_base=$(basename ${3} | rev | cut -d'.' -f2 | rev)
ref_dir=$(dirname ${3})


bowtie2-build ${3} ${ref_dir}/${ref_base}.bt2 #/scicomp/groups/OID/NCEZID/DHQP/CEMB/MiSeqAnalysisFiles/hypVir/pVir-CR-HvKP4.bt2

if [[ ! -d ${processed}/ ]]; then
	mkdir -p ${3}/BAMs
fi

while IFS= read -r var; do
	sample_name=$(echo "${var}" | awk -F/ '{ print $2}' | tr -d '[:space:]')
	project=$(echo "${var}" | awk -F/ '{ print $1}' | tr -d '[:space:]')
	OUTDATADIR="${processed}/${project}/${sample_name}"

	bowtie2 -f -p 10 -x ${ref_dir}/${ref_base} -U ${OUTDATADIR}/Assembly/${sample_name}_scaffolds_trimmed.fasta > ${2}/BAMs/{sample_name}.sam

	samtools view -bS ${2}/BAMs/${sample_name}.sam | samtools sort - ${2}/BAMs/${sample_name}.sam

	samtools index ${2}/BAMs/${sample_name}.bam ${2}/BAMs${sample_name}.bai

	perl "${shareScript}/estimate_core_genome_from_bam.pl" -bam "${2}/BAMs/" -genome "${3}" -depth 10

done < "${1}"
echo "All runs completed"
global_end_time=$(date "+%m-%d-%Y @ %Hh_%Mm_%Ss")
#Script exited gracefully (unless something else inside failed)
printf "%s %s" "Act_by_list.sh has completed whatever it was doing at" "${global_end_time}" | mail -s "act_by_list complete" nvx4@cdc.gov
exit 0
