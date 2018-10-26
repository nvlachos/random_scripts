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
fi


# Loop through and act on each sample name in the passed/provided list
while IFS= read -r line_in; do
	line_in=$(echo "${line_in}" | tr -d '\n')
	if [[ "${2}" = "1" ]]; then
		IFS='	' read -r -a array_1st <<< "$line_in"
		assembly_source="${processed}/${array_1st[1]}/${array_1st[2]}/Assembly/${array_1st[2]}_scaffolds_trimmed.fasta"
		mlst_db="${array_alt[0]}"
		mlst_type="${array_1st[3]}"
		locus1="${array_1st[4]}"
		locus2="${array_1st[5]}"
		locus3="${array_1st[6]}"
		locus4="${array_1st[7]}"
		locus5="${array_1st[8]}"
		locus6="${array_1st[9]}"
		locus7="${array_1st[10]}"
		echo -e "${assembly_source}	${array_1st[0]}	${mlst_type}	${locus1}	${locus2}	${locus3}	${locus4}	${locus5}	${locus6}	${locus7}"
		echo -e "${assembly_source}	${array_1st[0]}	${mlts_type}	${locus1}	${locus2}	${locus3}	${locus4}	${locus5}	${locus6}	${locus7}" > "${processed}/${array_1st[1]}/${array_1st[2]}/MLST/${array_1st[2]}.mlst"
	elif [[ "${2}" = "2" ]]; then
		IFS='	' read -r -a array_alt <<< "$line_in"
		header=$(head -n1 "${processed}/${array_alt[1]}/${array_alt[2]}/MLST/${array_alt[2]}_${array_alt[0]}.mlst")
		assembly_source="${processed}/${array_alt[1]}/${array_alt[2]}/Assembly/${array_alt[2]}_scaffolds_trimmed.fasta"
		mlst_db="${array_alt[0]}"
		mlst_type="${array_alt[3]}"
		echo "${array_alt[@]}"
		locus1=$(echo "${array_alt[4]}" | cut -d'(' -f2 | cut -d')' -f1)
		locus2=$(echo "${array_alt[5]}" | cut -d'(' -f2 | cut -d')' -f1)
		locus3=$(echo "${array_alt[6]}" | cut -d'(' -f2 | cut -d')' -f1)
		locus4=$(echo "${array_alt[7]}" | cut -d'(' -f2 | cut -d')' -f1)
		locus5=$(echo "${array_alt[8]}" | cut -d'(' -f2 | cut -d')' -f1)
		locus6=$(echo "${array_alt[9]}" | cut -d'(' -f2 | cut -d')' -f1)
		locus7=$(echo "${array_alt[10]}" | cut -d'(' -f2 | cut -d')' -f1)
		echo -e "${header}\n${assembly_source}	${mlst_db}	${mlst_type}	${locus1}	${locus2}	${locus3}	${locus4}	${locus5}	${locus6}	${locus7}"
		echo -e "${header}\n${assembly_source}	${mlst_db}	${mlst_type}	${locus1}	${locus2}	${locus3}	${locus4}	${locus5}	${locus6}	${locus7}" > "${processed}/${array_alt[1]}/${array_alt[2]}/MLST/${array_alt[2]}_${array_alt[0]}.mlst"
	else
		echo "Unknown 2nd parameter: ${2}, should be 1 (for original MLST) or 2 (for alternate MLST)"
	fi
done < ${share}/${1}
