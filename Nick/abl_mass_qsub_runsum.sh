#!/bin/sh -l

#$ -o ablrunsum-cs.out
#$ -e ablrunsum-cs.err
#$ -N ablrunsum-cs
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
. ./config.sh
#Import the module file that loads all necessary mods
. "${mod_changers}/pipeline_mods"

#
# Usage ./abl_mass_qsub_srst2.sh path_to_list max_concurrent_submissions output_directory_for_scripts
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to abl_mass_qsub_runsum.sh, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./abl_mass_qsub_runsum.sh path_to_list_file(single runs per line) max_concurrent_submissions output_directory_for_scripts"
	exit 1
elif [[ ! -f "${1}" ]]; then
	echo "${1} (list) does not exist...exiting"
	exit 1
elif [[ ! -d "${3}" ]]; then
	echo "${3} location does not exist...exiting"
	exit 1
fi

# create an array of all samples in the list
arr=()
#while IFS= read -r line || [[ "$line" ]];  do
#	echo ${line}
#	arr+=$(echo "$line")
#done < ${1}

IFS=$'\n' read -r -a arr < ${1}

arr_size="${#arr[@]}"
last_index=$(( arr_size -1 ))
echo "-${arr_size}:${arr[@]}-"
exit
# Create counter and set max number of concurrent submissions
counter=0
max_subs=${2}

# Set script directory
main_dir="${3}/runsum_subs"
if [[ ! -d "${3}/runsum_subs" ]]; then
	mkdir "${3}/runsum_subs"
	mkdir "${3}/runsum_subs/complete"
elif [[ ! -d "${3}/runsum_subs/complete" ]]; then
	mkdir "${3}/runsum_subs/complete"
fi

time_run=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")

# Creates and submits qsub scripts to perform summaries of each run in the list
while [ ${counter} -lt ${arr_size} ] ; do
	project=${arr[${counter}]}
	#echo ${counter}"-${processed}/${project}/${project}/kraken/postAssembly/${project}_kraken_summary_assembled_BP_data.txt"
	if [[ ${counter} -lt ${max_subs} ]]; then
		echo  "Index is below max submissions, submitting"
		echo -e "#!/bin/bash -l\n" > "${main_dir}/runsum_${project}_${time_run}.sh"
		echo "Saving to ${main_dir}/runsum_${project}_${time_run}.sh"
		echo -e "#$ -o runsum_${project}.out" >> "${main_dir}/runsum_${project}_${time_run}.sh"
		echo -e "#$ -e runsum_${project}.err" >> "${main_dir}/runsum_${project}_${time_run}.sh"
		echo -e "#$ -N runsum_${project}"   >> "${main_dir}/runsum_${project}_${time_run}.sh"
		echo -e "#$ -cwd"  >> "${main_dir}/runsum_${project}_${time_run}.sh"
		echo -e "#$ -q short.q\n"  >> "${main_dir}/runsum_${project}_${time_run}.sh"
		echo -e "cd ${shareScript}" >> "${main_dir}/runsum_${project}_${time_run}.sh"
		echo -e "\"${shareScript}/run_sum.sh\" \"${project}\"" >> "${main_dir}/runsum_${project}_${time_run}.sh"
		echo -e "echo \"$(date)\" > \"${main_dir}/complete/${project}_runsum_complete.txt\"" >> "${main_dir}/runsum_${project}_${time_run}.sh"
		cd "${main_dir}"
		if [[ "${counter}" -lt "${last_index}" ]]; then
			qsub "${main_dir}/runsum_${project}_${time_run}.sh"
		else
			qsub -sync y "${main_dir}/runsum_${project}_${time_run}.sh"
		fi
		mv "${shareScript}/runsum_${project}.out" ${main_dir}
		mv "${shareScript}/runsum_${project}.err" ${main_dir}
	else
		waiting_for_index=$(( counter - max_subs ))
		waiting_project=$(echo "${arr[${waiting_for_index}]}" | cut -d'/' -f2)
		timer=0
		echo "Index is above max submissions, waiting for index ${waiting_for_index}:${waiting_project} to complete"
		while :
		do
			if [[ ${timer} -gt 1800 ]]; then
				echo "Timer exceeded limit of 1800 seconds 30 minutes"
				break
			fi
			if [[ -f "${main_dir}/complete/${waiting_project}_runsum_complete.txt" ]]; then
				echo  "Index is below max submissions, submitting"
				echo -e "#!/bin/bash -l\n" > "${main_dir}/runsum_${project}_${time_run}.sh"
				echo -e "#$ -o runsum_${project}.out" >> "${main_dir}/runsum_${project}_${time_run}.sh"
				echo -e "#$ -e runsum_${project}.err" >> "${main_dir}/runsum_${project}_${time_run}.sh"
				echo -e "#$ -N runsum_${project}"   >> "${main_dir}/runsum_${project}_${time_run}.sh"
				echo -e "#$ -cwd"  >> "${main_dir}/runsum_${project}_${time_run}.sh"
				echo -e "#$ -q short.q\n"  >> "${main_dir}/runsum_${project}_${time_run}.sh"
				echo -e "cd ${shareScript}" >> "${main_dir}/runsum_${project}_${time_run}.sh"
				echo -e "\"${shareScript}/run_runsum.sh\" \"${project}\"" >> "${main_dir}/runsum_${project}_${time_run}.sh"
				echo -e "echo \"$(date)\" > \"${main_dir}/complete/${project}_runsum_complete.txt\"" >> "${main_dir}/runsum_${project}_${time_run}.sh"
				cd "${main_dir}"
				if [[ "${counter}" -lt "${last_index}" ]]; then
					qsub "${main_dir}/runsum_${project}_${time_run}.sh"
				else
					qsub -sync y "${main_dir}/runsum_${project}_${time_run}.sh"
				fi
				mv "${shareScript}/runsum_${project}.out" ${main_dir}
				mv "${shareScript}/runsum_${project}.err" ${main_dir}
				break
			else
				timer=$(( timer + 5 ))
				echo "sleeping for 5 seconds, so far slept for ${timer}"
				sleep 5
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
