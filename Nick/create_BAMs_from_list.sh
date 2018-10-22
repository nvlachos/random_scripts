#!/bin/sh -l

#$ -o act_by_list_barebones.out
#$ -e act_by_list_barebones.err
#$ -N ablb
#$ -cwd
#$ -q all.q

#Import the config file with shortcuts and settings
. /scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/scripts/config.sh
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
	echo "Usage is ./act_by_list.sh path_to_list_file(single sample ID per line, e.g. B8VHY/1700128 (it must include project id also))"
	echo "Output location varies depending on which tasks are performed but will be found somewhere under ${processed}"
	exit 0
fi

# Loop through and act on each sample name in the passed/provided list

bowtie2-build /scicomp/groups/OID/NCEZID/DHQP/CEMB/Tom_DIR/scripts/pVir-CR-HvKP4-FOIA2018/pVir-CR-HvKP4.fasta /scicomp/groups/OID/NCEZID/DHQP/CEMB/MiSeqAnalysisFiles/hypVir/pVir-CR-HvKP4.bt2

if [[ ! -d ${processed}/ ]]; then
	mkdir ${share}/projects/CR-Acinetobacter_baumanni/BAMs
fi

while IFS= read -r var; do
	sample_name=$(echo "${var}" | awk -F/ '{ print $2}' | tr -d '[:space:]')
	project=$(echo "${var}" | awk -F/ '{ print $1}' | tr -d '[:space:]')
	OUTDATADIR="${processed}/${project}/${sample_name}"
	
	bowtie2 -f -p 10 -x ${share}/projects/CR-Acinetobacter_baumanni/BAMs/.bt2 -U ${OUTDATADIR}/Assembly/${sample_name}_scaffolds_trimmed.fasta > ${share}/projects/CR-Acinetobacter_baumanni/BAMs/{sample_name}.sam

	samtools view -bS ${share}/projects/CR-Acinetobacter_baumanni/BAMs/${sample_name}.sam | samtools sort - ${share}/projects/CR-Acinetobacter_baumanni/BAMs/${sample_name}.sam
	
	samtools index ${share}/projects/CR-Acinetobacter_baumanni/BAMs/${sample_name}.bam ${share}/projects/CR-Acinetobacter_baumanni/BAMs${sample_name}.bai
	
	perl "${shareScript}/estimate_core_genome_from_bam.pl" -bam "${share}/projects/CR-Acinetobacter_baumanni/BAMs/" -genome "/scicomp/groups/OID/NCEZID/DHQP/CEMB/Tom_DIR/scripts/pVir-CR-HvKP4-FOIA2018/pVir-CR-HvKP4.fasta" -depth 10 
	#2> "${share}/Phylogeny_analyses/${1}/LyveSET/out/genome_core.txt"

	done < "${share}/${1}"
echo "All runs completed"
global_end_time=$(date "+%m-%d-%Y @ %Hh_%Mm_%Ss")
#Script exited gracefully (unless something else inside failed)
printf "%s %s" "Act_by_list.sh has completed whatever it was doing at" "${global_end_time}" | mail -s "act_by_list complete" nvx4@cdc.gov
exit 0
