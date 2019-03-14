#!/bin/sh -l

#$ -o act_by_list_barebones1.out
#$ -e act_by_list_barebones1.err
#$ -N ablb1
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
. ./config.sh
#Import the module file that loads all necessary mods
. "${mod_changers}/pipeline_mods"

#
# Usage ./act_by_list.sh path_to_list max_concurrent_submission output_directory_for_scripts
#
#
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to abl_mass_sub_node.sh, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./abl_mass_qsub_node.sh path_to_list_file(single sample ID per line, e.g. B8VHY/1700128 (it must include project id also)) max_submissions output_directory_for_scripts"
	echo "Output location varies depending on which tasks are performed but will be found somewhere under ${processed}"
	exit 0
fi

arr=()
while IFS= read -r line || [[ "$line" ]];  do
  arr+=("$line")
done < ${1}

arr_size="${#arr[@]}"
last_index=$(( arr_size -1 ))
echo "-${arr_size}:${arr[@]}-"

if [[ -z ${2} ]]; then
	max_subs=10
else
	max_subs=${2}
fi

# Loop through and act on each sample name in the passed/provided list
counter=0
main_dir="${3}/node_subs"
if [[ ! -d "${3}/node_subs" ]]; then
	mkdir "${3}/node_subs"
	mkdir "${3}/node_subs/complete"
elif [[ ! -d  "${3}/node_subs/complete" ]]; then
	mkdir "${3}/node_subs/complete"
fi

while [ ${counter} -lt ${arr_size} ] ; do
	sample=$(echo "${arr[${counter}]}" | cut -d'/' -f2)
	project=$(echo "${arr[${counter}]}" | cut -d'/' -f1)
	echo ${counter}
	if [ ${counter} -lt ${max_subs} ]; then
		if [[ -s "${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed.fasta" ]]; then
			header=$(head -n1 "${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed.fasta" | cut -d'_' -f1)
			if [[ "${header}" = ">NODE" ]]; then
				echo  "Index is below max submissions, submitting"
				echo -e "#!/bin/bash -l\n" > "${main_dir}/node_${sample}_${start_time}.sh"
				echo -e "#$ -o node_${sample}.out" >> "${main_dir}/node_${sample}_${start_time}.sh"
				echo -e "#$ -e node_${sample}.err" >> "${main_dir}/node_${sample}_${start_time}.sh"
				echo -e "#$ -N node_${sample}"   >> "${main_dir}/node_${sample}_${start_time}.sh"
				echo -e "#$ -cwd"  >> "${main_dir}/node_${sample}_${start_time}.sh"
				echo -e "#$ -q short.q\n"  >> "${main_dir}/node_${sample}_${start_time}.sh"
				echo -e "mv \"${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed.fasta\" \"${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed_original.fasta\"" >> "${main_dir}/node_${sample}_${start_time}.sh"
				echo -e "python3 \"${shareScript}/fasta_headers.py\" \"${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed_original.fasta\" \"${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed.fasta\"" >> "${main_dir}/node_${sample}_${start_time}.sh"
				echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_node_complete.txt\"" >> "${main_dir}/node_${sample}_${start_time}.sh"
				cd "${main_dir}"
				if [[ "${counter}" -lt "${last_index}" ]]; then
					qsub "${main_dir}/node_${sample}_${start_time}.sh"
				else
					qsub -sync y "${main_dir}/node_${sample}_${start_time}.sh"
				fi
			else
				echo "${project}/${sample} already had its nodes removed"
				echo "$(date)" > "${main_dir}/complete/${sample}_node_complete.txt"
			fi
		else
			echo "${project}/${sample} has no Assembly to fix"
			echo "$(date)" > "${main_dir}/complete/${sample}_node_complete.txt"
		fi
	else
		waiting_for_index=$(( counter - max_subs ))
		waiting_sample=$(echo "${arr[${waiting_for_index}]}" | cut -d'/' -f2)
		timer=0
		echo "Index is above max submissions, waiting for index ${waiting_for_index}:${waiting_sample} to complete"
		while :
		do
			if [[ ${timer} -gt 1800 ]]; then
				echo "Timer exceeded limit of 1800 seconds 30 minutes"
				break
			fi
			if [ -f "${main_dir}/complete/${waiting_sample}_node_complete.txt" ]; then
				if [[ -s "${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed.fasta" ]]; then
					header=$(head -n1 "${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed.fasta" | cut -d'_' -f1)
					if [[ "${header}" = ">NODE" ]]; then
						echo  "Index is below max submissions, submitting"
						echo -e "#!/bin/bash -l\n" > "${main_dir}/node_${sample}_${start_time}.sh"
						echo -e "#$ -o node_${sample}.out" >> "${main_dir}/node_${sample}_${start_time}.sh"
						echo -e "#$ -e node_${sample}.err" >> "${main_dir}/node_${sample}_${start_time}.sh"
						echo -e "#$ -N node_${sample}"   >> "${main_dir}/node_${sample}_${start_time}.sh"
						echo -e "#$ -cwd"  >> "${main_dir}/node_${sample}_${start_time}.sh"
						echo -e "#$ -q short.q\n"  >> "${main_dir}/node_${sample}_${start_time}.sh"
						echo -e "mv \"${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed.fasta\" \"${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed_original.fasta\"" >> "${main_dir}/node_${sample}_${start_time}.sh"
						echo -e "python3 \"${shareScript}/fasta_headers.py\" \"${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed_original.fasta\" \"${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed.fasta\"" >> "${main_dir}/node_${sample}_${start_time}.sh"
						echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_node_complete.txt\"" >> "${main_dir}/node_${sample}_${start_time}.sh"
						cd "${main_dir}"
						if [[ "${counter}" -lt "${last_index}" ]]; then
							qsub "${main_dir}/node_${sample}_${start_time}.sh"
						else
							qsub -sync y "${main_dir}/node_${sample}_${start_time}.sh"
						fi
						break
					else
						echo "${project}/${sample} already had its nodes removed"
						echo "$(date)" > "${main_dir}/complete/${sample}_node_complete.txt"
						break
					fi
				else
					echo "${project}/${sample} has no Assembly to fix"
					echo "$(date)" > "${main_dir}/complete/${sample}_node_complete.txt"
					break
				fi
			fi
		done
	fi
	counter=$(( counter + 1 ))
done

echo "All isolates completed"
global_end_time=$(date "+%m-%d-%Y @ %Hh_%Mm_%Ss")
#Script exited gracefully (unless something else inside failed)
printf "%s %s" "Act_by_list.sh has completed ${2}" "${global_end_time}" | mail -s "act_by_list complete" nvx4@cdc.gov
exit 0
