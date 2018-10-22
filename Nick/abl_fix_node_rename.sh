#!/bin/sh -l

#$ -o ablmq-fnr.out
#$ -e ablmq-fnr.err
#$ -N ablmq-fnr
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
. /scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/scripts/config.sh
#Import the module file that loads all necessary mods
. "${mod_changers}/pipeline_mods"

#
# Usage ./act_by_list.sh list_name(currently has to be placed in /scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR folder) description of list function
#
# script changes depending on what needs to be run through the list
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to act_by_list.sh, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./abl_mass_qsub_template.sh path_to_list_file(single sample ID per line, e.g. B8VHY/1700128 (it must include project id also)) max_submissions"
	echo "Output location varies depending on which tasks are performed but will be found somewhere under ${processed}"
	exit 0
fi

arr=()
while IFS= read -r line || [[ "$line" ]];  do
  arr+=("$line")
done < ${1}

arr_size="${#arr[@]}"
echo "-${arr_size}:${arr[@]}-"

if [[ ! -z "${2}" ]]; then
	max_subs="${2}"
else
	max_subs=10
fi

# Loop through and act on each sample name in the passed/provided list
counter=0
max_subs=${2}
main_dir="${share}/mass_subs/csstar_subs"
if [[ ! -d "${share}/mass_subs/csstar_subs" ]]; then
	mkdir -p "${share}/mass_subs/csstar_subs/complete"
elif [[ ! -d  "${share}/mass_subs/csstar_subs/complete" ]]; then
	mkdir "${share}/mass_subs/csstar_subs/complete"
fi

start_time=$(date)

#while [ ${counter} -lt ${arr_size} ] ; do
#	sample=$(echo "${arr[${counter}]}" | cut -d'/' -f2)
#	project=$(echo "${arr[${counter}]}" | cut -d'/' -f1)
#	rm -r ${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed.fasta
#	rm -r ${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed_original.fasta
#	python ${shareScript}/removeShortContigs.py ${processed}/${project}/${sample}/Assembly/scaffolds.fasta 500
#	mv ${processed}/${project}/${sample}/Assembly/scaffolds.fasta.TRIMMED.fasta ${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed_original.fasta
#	python ${shareScript}/fasta_headers.py ${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed_original.fasta ${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed.fasta
#	counter=$(( counter + 1))
#done

while [ ${counter} -lt ${arr_size} ] ; do
	sample=$(echo "${arr[${counter}]}" | cut -d'/' -f2)
	project=$(echo "${arr[${counter}]}" | cut -d'/' -f1)
	rm -r ${processed}/${project}/${sample}/plasmidAssembly/${sample}_plasmid_scaffolds_trimmed.fasta
	rm -r ${processed}/${project}/${sample}/plasmidAssembly/${sample}_plasmid_scaffolds_trimmed_original.fasta
	python ${shareScript}/removeShortContigs.py ${processed}/${project}/${sample}/plasmidAssembly/scaffolds.fasta 500
	mv ${processed}/${project}/${sample}/plasmidAssembly/scaffolds.fasta.TRIMMED.fasta ${processed}/${project}/${sample}/plasmidAssembly/${sample}_plasmid_scaffolds_trimmed_original.fasta
	python "${shareScript}/fasta_headers.py" ${processed}/${project}/${sample}/plasmidAssembly/${sample}_plasmid_scaffolds_trimmed_original.fasta ${processed}/${project}/${sample}/plasmidAssembly/${sample}_plasmid_scaffolds_trimmed.fasta
	counter=$(( counter + 1))
done

