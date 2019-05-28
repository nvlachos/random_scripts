#!/bin/sh -l

#$ -o ablmq-fnr.out
#$ -e ablmq-fnr.err
#$ -N ablmq-fnr
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
. ./config.sh
#Import the module file that loads all necessary mods
. "${mod_changers}/pipeline_mods"

#
# Usage ./fix_node_rename.sh list_name(currently has to be placed in /scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR folder) description of list function
#
# script changes depending on what needs to be run through the list
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./abl_mass_qsub_template.sh path_to_list_file(single sample ID per line, e.g. B8VHY/1700128 (it must include project id also)) max_submissions"
	echo "Output location varies depending on which tasks are performed but will be found somewhere under ${processed}"
	exit 0
fi

# Create an array of all samples in the list
arr=()
while IFS= read -r line || [[ "$line" ]];  do
  arr+=("$line")
done < ${1}

arr_size="${#arr[@]}"
echo "-${arr_size}:${arr[@]}-"

# Fix all normal assemblies
while [ ${counter} -lt ${arr_size} ] ; do
	sample=$(echo "${arr[${counter}]}" | cut -d'/' -f2)
	project=$(echo "${arr[${counter}]}" | cut -d'/' -f1)
	if [[ -f ${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed.fasta ]]; then
		rm -r ${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed.fasta
	fi
	if [[ -f ${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed_original.fasta]]; then
		rm -r ${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed_original.fasta
	fi
	if [[ -f ${processed}/${project}/${sample}/Assembly/scaffolds.fasta ]]; then
		python ${shareScript}/removeShortContigs.py -i ${processed}/${project}/${sample}/Assembly/scaffolds.fasta -t 500 -s "normal_SPAdes"
		mv ${processed}/${project}/${sample}/Assembly/scaffolds.fasta.TRIMMED.fasta ${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed_original.fasta
		python ${shareScript}/fasta_headers.py -i ${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed_original.fasta -o ${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed.fasta
	else
		echo "No scaffolds file found for ${sample}"
	fi
	counter=$(( counter + 1))
done

# # Fix all plasmid assemblies
# while [ ${counter} -lt ${arr_size} ] ; do
# 	sample=$(echo "${arr[${counter}]}" | cut -d'/' -f2)
# 	project=$(echo "${arr[${counter}]}" | cut -d'/' -f1)
# 	if [[ -f ${processed}/${project}/${sample}/plasmidAssembly/${sample}_plasmid_scaffolds_trimmed.fasta ]]; then
# 		rm -r ${processed}/${project}/${sample}/plasmidAssembly/${sample}_plasmid_scaffolds_trimmed.fasta
# 	fi
# 	if [[ -f rm -r ${processed}/${project}/${sample}/plasmidAssembly/${sample}_plasmid_scaffolds_trimmed_original.fasta ]]; then
# 		rm -r ${processed}/${project}/${sample}/plasmidAssembly/${sample}_plasmid_scaffolds_trimmed_original.fasta
# 	fi
# 	if [[ -f ${processed}/${project}/${sample}/plasmidAssembly/scaffolds.fasta ]]; then
# 		python ${shareScript}/removeShortContigs.py -i ${processed}/${project}/${sample}/plasmidAssembly/scaffolds.fasta -t 500 -s "normal_SPAdes"
# 		mv ${processed}/${project}/${sample}/plasmidAssembly/scaffolds.fasta.TRIMMED.fasta ${processed}/${project}/${sample}/plasmidAssembly/${sample}_plasmid_scaffolds_trimmed_original.fasta
# 		python "${shareScript}/fasta_headers.py" -i ${processed}/${project}/${sample}/plasmidAssembly/${sample}_plasmid_scaffolds_trimmed_original.fasta -o ${processed}/${project}/${sample}/plasmidAssembly/${sample}_plasmid_scaffolds_trimmed.fasta
# 	else
# 		echo "No plasmid scaffolds found for ${sample}"
# 	counter=$(( counter + 1))
# done
