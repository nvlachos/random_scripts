#!/bin/sh -l

#$ -o abl-gracon.out
#$ -e abl-gracon.err
#$ -N abl-gracon
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
. ./config.sh
#Import the module file that loads all necessary mods
. "${mod_changers}/pipeline_mods"

#
# Usage ./act_by_list_graph_consolidater.sh path_to_list path_for_output_file
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./act_by_list_graph_consolidater.sh path_for_output_file"
	exit 0
elif [[ -z "${2}" ]]; then
	echo  "No output file input, exiting..."
	exit 1
fi

# Loop through and act on each sample name in the passed/provided list
counter=0
> ${2}
while IFS= read -r var; do
	sample_name=$(echo "${var}" | cut -d'/' -f2 | tr -d '[:space:]')
	project=$(echo "${var}" | cut -d'/' -f1 | tr -d '[:space:]')
	#echo "${counter}"
	if [[ -f "${processed}/${project}/${sample_name}/plasFlow/Unicycler_assemblies/${sample_name}_uni_assembly/${sample_name}_assembly.gfa" ]]; then
		while IFS= read plasFlow_graph; do
			line_type=$(echo "${plasFlow_graph}" | cut -d'	' -f1)
			if [[ "${line_type}" == "H" ]]; then
				echo "Don't know how H is used yet, skipping"
			elif [[ "${line_type}" == "S" ]]; then
				contig_num=$(echo "${plasFlow_graph}" | cut -d'	' -f2)
				seq=$(echo "${plasFlow_graph}" | cut -d'	' -f3)
				new_line="${line_type}	${sample_name}_${contig_num}	${seq}"
			elif [[ "${line_type}" == "L" ]]; then
				from=$(echo "${plasFlow_graph}" | cut -d'	' -f2)
				to=$(echo "${plasFlow_graph}" | cut -d'	' -f4)
				fromOrient=$(echo "${plasFlow_graph}" | cut -d'	' -f3)
				toOrient=$(echo "${plasFlow_graph}" | cut -d'	' -f5)
				overlap=$(echo "${plasFlow_graph}" | cut -d'	' -f6)
				new_line="${line_type}	${sample_name}_${from}	${fromOrient}	${sample_name}_${to}	${toOrient}	${overlap}"
			elif [[ "${line_type}" == "C" ]]; then
				container=$(echo "${plasFlow_graph}" | cut -d'	' -f2)
				contained=$(echo "${plasFlow_graph}" | cut -d'	' -f4)
				containerOrient=$(echo "${plasFlow_graph}" | cut -d'	' -f3)
				containedOrient=$(echo "${plasFlow_graph}" | cut -d'	' -f5)
				POS=$(echo "${plasFlow_graph}" | cut -d'	' -f6)
				overlap=$(echo "${plasFlow_graph}" | cut -d'	' -f7)
				new_line="${line_type}	${sample_name}_${container}	${containerOrient} ${sample_name}_${contained}	${containedOrient}	${POS}	${overlap}"
			elif [[ "${line_type}" == "P" ]]; then
				pathname=$(echo "${plasFlow_graph}" | cut -d'	' -f2)
				segment_names=$(echo "${plasFlow_graph}" | cut -d'	' -f3)
				new_segment_array=""
				IFS=', ' read -r -a segment_array <<< "${segment_names}"
				for segment in segment_array; do
					new_segment_array=new_segment_array+"${sample_name}_${segment},"
				done
				overlaps=$(echo "${plasFlow_graph}" | cut -d'	' -f4)
				new_line="${line_type}	${sample_name}_${pathname}	${segment::-1}	${overlaps}"
			else
				echo "Unknown line type - ${line_type}, skipping"
				continue
			fi
			echo "${new_line}" >> ${2}
		done < "${processed}/${project}/${sample_name}/plasFlow/Unicycler_assemblies/${sample_name}_uni_assembly/${sample_name}_assembly.gfa"
	fi
	#counter=$(( counter + 1 ))
done < "${1}"

echo "All isolates completed"
global_end_time=$(date "+%m-%d-%Y @ %Hh_%Mm_%Ss")
#Script exited gracefully (unless something else inside failed)
printf "%s %s" "act_by_list_graph_consolidater.sh has completed check of snd MLSTs" "${global_end_time}" | mail -s "act_by_list complete" nvx4@cdc.gov
exit 0
