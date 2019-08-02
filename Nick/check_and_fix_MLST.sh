#!/bin/sh -l

#$ -o fix_MLST.out
#$ -e fix_MLST.err
#$ -N fix_MLST
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
. ./config.sh
#Import the module file that loads all necessary mods
. "${mod_changers}/pipeline_mods"

#
# Usage ./check_and_fix_MLST.sh list_name
#
# script changes depending on what needs to be run through the list
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./check_and_fix_MLST.sh path_to_list"
	#echo "Output location varies depending on which tasks are performed but will be found somewhere under ${processed}"
	exit 0
elif [[ ! -f ${1} ]]; then
	echo "List file does not exist, exiting"
	exit 1
fi


# Loop through and act on each sample name in the passed/provided list
while IFS= read -r line_in; do
	#line_in=$(echo "${line_in}" | tr -d '\n')
	sample=$(echo "${line_in}" | cut -d'/' -f2)
	project=$(echo "${line_in}" | cut -d'/' -f1)
	# Check main automated mlst file first
	for mlst_file in "${processed}/${project}/${sample}/MLST/"*; do
		echo "checking ${mlst_file}"
		if [[ "${mlst_file}" == *".mlst" ]]; then
			if [[ "${mlst_file}" == *"srst2"* ]]; then
				#echo "Dont have srst2 checker yet, need to find good srst2 files"
				#echo "${mlst_file}" >> /scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/srst2_mlsts.txt
				python3 "${shareScript}/check_and_fix_MLST.py" -i "${mlst_file}" -t "srst2"
			else
				echo "About to do ${mlst_file}"
				lines=$(wc -l < ${mlst_file})
				if [[ "${lines}" -eq 2 ]]; then
					echo "Need to redo mlst file"
					DB=$(tail -n1 ${mlst_file} | cut -d'	' -f2)
					"${shareScript}/run_MLST.sh" "${sample}" "${project}" "-f" "${DB}"
				fi
				python3 "${shareScript}/check_and_fix_MLST.py" -i "${mlst_file}" -t "standard"
			fi
		fi # Done with main mlst file
	done
done < ${1}
