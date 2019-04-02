#!/bin/sh -l

#$ -o fix_MLST.out
#$ -e fix_MLST.err
#$ -N fix_MLST
#$ -cwd
#$ -q all.q

#Import the config file with shortcuts and settings
. ./config.sh
#Import the module file that loads all necessary mods
. "${mod_changers}/pipeline_mods"

#
# Usage ./fix_MLST.sh list_name(currently has to be placed in /scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR folder) -type (1 for automatically detected MLST, 2 for forced MLST against a certain database)
#
# script changes depending on what needs to be run through the list
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to act_by_list.sh, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./fix_2ns_mlst.sh path_to_corrected_MLST_file(DB	MiSeq_ID	sample_ID	MLST	locus_1(type)	locus2(type) ... locus7(type)) -mlst_type (1 for normal, 2 for alternate(forced))"
	#echo "Output location varies depending on which tasks are performed but will be found somewhere under ${processed}"
	exit 0
elif [[ ! -f ${1} ]]; then
	echo "List file does not exist, exiting"
	exit 1
fi


# Loop through and act on each sample name in the passed/provided list
while IFS= read -r line_in; do
	line_in=$(echo "${line_in}" | tr -d '\n')
	sample=$(echo "${line_in}" | cut -d'/' -f2)
	project=$(echo "${line_in}" | cut -d'/' -f1)
	# Check main automated mlst file first
	if [[ -f "${processed}/${project}/${sample}/MLST/${sample}.mlst" ]]; then
		echo "Starting standard"
		STtype=$(tail -n1 "${processed}/${project}/${sample}/MLST/${sample}.mlst" | cut -d'	' -f3)
		mlst_DB=$(tail -n1 "${processed}/${project}/${sample}/MLST/${sample}.mlst" | cut -d'	' -f2)
		echo "ST-${STtype} from ${mlst_DB}"
		IFS='	' read -r -a source_allele_array <<< $(tail -n1 "${processed}/${project}/${sample}/MLST/${sample}.mlst" | cut -d'	' -f4-)
		allele_count=${#aource_allele_array}
		computed_allele_array=()
		for allele in ${source_allele_array[@]}; do
			allele_name=$(echo ${allele} | cut -d'(' -f1)
			allele_type=$(echo ${allele} | cut -d'(' -f2 | cut -d')' -f1)
			echo "${allele}:${allele_name}:${allele_type}"
			if [[ "${allele_type}" == *","* ]] || [[ "${allele_type}" == *"/"* ]]; then
				echo "DUAL ALLELE"
				IFS=',' read -r -a multilocus_array <<< "${allele_type}"
				locus_array=(${allele_name} ${multilocus_array[@]})
			else
				locus_array=(${allele_name} ${allele_type})
			fi
			echo "LA: ${locus_array[@]}"
		done
		computed_allele_array+=("${locus_array}")
	fi # Done with main mlst file
	counter=0
	echo "CAA=${#computed_allele_array}"
	for i in ${#computed_allele_array}; do
		echo "${counter}-${allele[$counter]}"
		counter=$((counter + 1 ))
	done
done < ${1}
