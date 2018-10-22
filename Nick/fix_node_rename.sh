#!/bin/sh -l

#$ -o fix_NODE.out
#$ -e fix_NODE.err
#$ -N fix_NODE
#$ -cwd
#$ -q all.q

#Import the config file with shortcuts and settings
. ./config.sh
#Import the module file that loads all necessary mods
. "${mod_changers}/pipeline_mods"

#
# Usage ./fix_node_rename.sh   sample_name	Run_ID
#
# script changes naming structure of SPAdes output to include isolate name for every contig and removes coverage info
#


# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to run_srst2.sh, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./node_rename.sh  sample_name Run_ID"
	echo "Output location is ${processed}/${2}/${1}/Assembly and ${processed}/${2}/plasmidAssembly"
	exit 0
fi

if [[ -s "${processed}/${2}/${1}/Assembly/scaffolds.fasta" ]];then
	if [[ -f "${processed}/${2}/${1}/Assembly/scaffolds.fasta.TRIMMED.fasta" ]]; then
		rm -r "${processed}/${2}/${1}/Assembly/scaffolds.fasta.TRIMMED.fasta"
	fi
	if [[ -f "${processed}/${2}/${1}/Assembly/${1}_scaffolds_trimmed.fasta" ]]; then
		rm -r "${processed}/${2}/${1}/Assembly/${1}_scaffolds_trimmed.fasta"
	fi
	if [[ -f "${processed}/${2}/${1}/Assembly/${1}_scaffolds_trimmed_original.fasta" ]]; then
		rm -r "${processed}/${2}/${1}/Assembly/${1}_scaffolds_trimmed_original.fasta"
	fi
	python "${shareScript}/removeShortContigs.py" "${processed}/${2}/${1}/Assembly/scaffolds.fasta" "500"
	mv "${processed}/${2}/${1}/Assembly/scaffolds.fasta.TRIMMED.fasta" "${processed}/${2}/${1}/Assembly/${1}_scaffolds_trimmed_original.fasta"
	#cp "${processed}/${2}/${1}/Assembly/${1}_scaffolds_trimmed_original.fasta" "${processed}/${2}/${1}/Assembly/${1}_scaffolds_trimmed.fasta"
	while IFS= read -r line;
	do
		first=${line::1}
		if [ "${first}" = ">" ]; then
			contig_ID=$(echo "${line}" | cut -d'_' -f2)
			length=$(echo "${line}" | cut -d'_' -f4)
			new_line=">${1}_${contig_ID}_length_${length}"
			echo "${new_line}" >> "${processed}/${2}/${1}/Assembly/${1}_scaffolds_trimmed.fasta"
		else
			echo "${line}" >> "${processed}/${2}/${1}/Assembly/${1}_scaffolds_trimmed.fasta"
		fi
	done < "${processed}/${2}/${1}/Assembly/${1}_scaffolds_trimmed_original.fasta"
else
	echo "${processed}/${2}/${1}/Assembly/${1}_scaffolds_trimmed.fasta doesn't have substance"
fi

if [[ -s "${processed}/${2}/${1}/plasmidAssembly/scaffolds.fasta" ]];then 
	if [[ -f "${processed}/${2}/${1}/plasmidAssembly/scaffolds.fasta.TRIMMED.fasta" ]]; then
		rm -r "${processed}/${2}/${1}/plasmidAssembly/scaffolds.fasta.TRIMMED.fasta"
	fi
	if [[ -f "${processed}/${2}/${1}/plasmidAssembly/${1}_plasmid_scaffolds_trimmed.fasta" ]]; then
		rm -r "${processed}/${2}/${1}/plasmidAssembly/${1}_plasmid_scaffolds_trimmed.fasta"
	fi
	if [[ -f "${processed}/${2}/${1}/plasmidAssembly/${1}_scaffolds_trimmed_original.fasta" ]]; then
		rm -r "${processed}/${2}/${1}/plasmidAssembly/${1}_scaffolds_trimmed_original.fasta"
	fi
	python "${shareScript}/removeShortContigs.py" "${processed}/${2}/${1}/plasmidAssembly/scaffolds.fasta" "500"
	mv "${processed}/${2}/${1}/plasmidAssembly/scaffolds.fasta.TRIMMED.fasta" "${processed}/${2}/${1}/plasmidAssembly/${1}_plasmid_scaffolds_trimmed_original.fasta"
	while IFS= read -r line;
	do
		# Grab first letter of line (indicating taxonomic level)
		first=${line::1}
		if [ "${first}" = ">" ]; then
			contig_ID=$(echo "${line}" | cut -d'_' -f2)
			length=$(echo "${line}" | cut -d'_' -f4)
			component=$(echo "${line}" | cut -d'_' -f8)
			new_line=">${1}_${contig_ID}_length_${length}_component_${component}"
			echo "${new_line}" >> "${processed}/${2}/${1}/plasmidAssembly/${1}_plasmid_scaffolds_trimmed.fasta"
		else
			echo "${line}" >> "${processed}/${2}/${1}/plasmidAssembly/${1}_plasmid_scaffolds_trimmed.fasta"
		fi
	done < "${processed}/${2}/${1}/plasmidAssembly/${1}_plasmid_scaffolds_trimmed_original.fasta"
else
	echo "${processed}/${2}/${1}/Assembly/${1}_plasmid_scaffolds_trimmed.fasta doesn't have substance"
fi

