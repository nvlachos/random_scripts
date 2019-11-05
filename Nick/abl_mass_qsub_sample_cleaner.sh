#!/bin/sh -l

#$ -o ablclean.out
#$ -e ablclean.err
#$ -N cleaner
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
. ./config.sh

#
# Usage ./abl_mass_qsub_sample_cleaner.sh path_to_list max_concurrent_submissions output_directory_for_scripts
#

# Number regex to test max concurrent submission parameter
number='^[0-9]+$'

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./abl_mass_qsub_sclean.sh path_to_list_file(single runs per line) max_concurrent_submissions output_directory_for_scripts"
	exit 1
elif [[ ! -f "${1}" ]]; then
	echo "${1} (list) does not exist...exiting"
	exit 1
elif ! [[ ${2} =~ $number ]] || [[ -z "${2}" ]]; then
	echo "${2} is not a number or is empty. Please input max number of concurrent qsub submissions...exiting"
	exit 2
elif [[ ! -d "${3}" ]]; then
	echo "${3} location does not exist...exiting"
	exit 1
fi

# create an array of all samples in the list
arr=()
while IFS= read -r line || [ -n "$line" ];  do
	#echo "L:${line}"
	if [[ ! -z "${line}" ]]; then
		line=$(echo ${line} | tr -d '\n' | tr -d '\r')
		arr+=($line)
	fi
	#echo "A:${arr[@]}"
done < ${1}


#IFS=$'\n' read -d '' -r -a arr < ${1}

#readarray arr < ${1}

arr_size="${#arr[@]}"
last_index=$(( arr_size -1 ))
echo "-${arr_size}:${arr[@]}-"

# Create counter and set max number of concurrent submissions
counter=0
max_subs=${2}

# Set script directory
main_dir="${3}/sclean_subs"
if [[ ! -d "${3}/sclean_subs" ]]; then
	mkdir "${3}/sclean_subs"
	mkdir "${3}/sclean_subs/complete"
elif [[ ! -d "${3}/sclean_subs/complete" ]]; then
	mkdir "${3}/sclean_subs/complete"
fi

time_run=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")

# Creates and submits qsub scripts to perform summaries of each run in the list
while [ ${counter} -lt ${arr_size} ] ; do
	project=${arr[${counter}]}
	if [[ ${counter} -lt ${max_subs} ]]; then
		echo  "Index is below max submissions, submitting"
		echo -e "#!/bin/bash -l\n" > "${main_dir}/sclean_${project}_${time_run}.sh"
		echo "Saving to ${main_dir}/sclean_${project}_${time_run}.sh"
		echo -e "#$ -o sclean_${project}.out" >> "${main_dir}/sclean_${project}_${time_run}.sh"
		echo -e "#$ -e sclean_${project}.err" >> "${main_dir}/sclean_${project}_${time_run}.sh"
		echo -e "#$ -N sclean_${project}"   >> "${main_dir}/sclean_${project}_${time_run}.sh"
		echo -e "#$ -cwd"  >> "${main_dir}/sclean_${project}_${time_run}.sh"
		echo -e "#$ -q short.q\n"  >> "${main_dir}/sclean_${project}_${time_run}.sh"
		echo -e "cd ${shareScript}" >> "${main_dir}/sclean_${project}_${time_run}.sh"
		echo -e "\"${shareScript}/sample_cleaner.sh\" \"${sample_name}\" \"${project}\"" >> "${main_dir}/sclean_${project}_${time_run}.sh"
		echo -e "echo \"$(date)\" > \"${main_dir}/complete/${project}_sclean_complete.txt\"" >> "${main_dir}/sclean_${project}_${time_run}.sh"
		cd "${main_dir}"
		qsub "${main_dir}/sclean_${project}_${time_run}.sh"
		fi
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
			if [[ -f "${main_dir}/complete/${waiting_project}_sclean_complete.txt" ]]; then
				echo  "Index is below max submissions, submitting"
				echo -e "#!/bin/bash -l\n" > "${main_dir}/sclean_${project}_${time_run}.sh"
				echo -e "#$ -o sclean_${project}.out" >> "${main_dir}/sclean_${project}_${time_run}.sh"
				echo -e "#$ -e sclean_${project}.err" >> "${main_dir}/sclean_${project}_${time_run}.sh"
				echo -e "#$ -N sclean_${project}"   >> "${main_dir}/sclean_${project}_${time_run}.sh"
				echo -e "#$ -cwd"  >> "${main_dir}/sclean_${project}_${time_run}.sh"
				echo -e "#$ -q short.q\n"  >> "${main_dir}/sclean_${project}_${time_run}.sh"
				echo -e "cd ${shareScript}" >> "${main_dir}/sclean_${project}_${time_run}.sh"
				echo -e "\"${shareScript}/sample_cleaner.sh\" \"${sample_name}\" \"${project}\"" >> "${main_dir}/sclean_${project}_${time_run}.sh"
				echo -e "echo \"$(date)\" > \"${main_dir}/complete/${project}_sclean_complete.txt\"" >> "${main_dir}/sclean_${project}_${time_run}.sh"
				cd "${main_dir}"
				qsub "${main_dir}/sclean_${project}_${time_run}.sh"
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

# Check for completion is done right now, but the qsub fiules still need to be moved
timer=0
for item in "${arr[@]}"; do
	waiting_sample=$(echo "${item}" | cut -d'/' -f2)
	if [[ -f "${main_dir}/complete/${waiting_sample}_sclean_complete.txt" ]]; then
		echo "${item} is complete"
		if [[ -f "${shareScript}/sclean_${waiting_sample}.out" ]]; then
			mv "${shareScript}/sclean_${waiting_sample}.out" "${main_dir}"
		fi
		if [[ -f "${shareScript}/sclean_${waiting_sample}.err" ]]; then
			mv "${shareScript}/sclean_${waiting_sample}.err" "${main_dir}"
		fi
	else
		while :
		do
				if [[ ${timer} -gt 3600 ]]; then
					echo "Timer exceeded limit of 3600 seconds = 60 minutes"
					exit 1
				fi
				if [[ -f "${main_dir}/complete/${waiting_sample}_sclean_complete.txt" ]]; then
					echo "${item} is complete"
					if [[ -f "${shareScript}/sclean_${waiting_sample}.out" ]]; then
						mv "${shareScript}/sclean_${waiting_sample}.out" "${main_dir}"
					fi
					if [[ -f "${shareScript}/sclean_${waiting_sample}.err" ]]; then
						mv "${shareScript}/sclean_${waiting_sample}.err" "${main_dir}"
					fi
					break
				else
					timer=$(( timer + 5 ))
					echo "sleeping for 5 seconds, so far slept for ${timer}"
					sleep 5
				fi
		done
	fi
done

echo "All isolates completed"
global_end_time=$(date "+%m-%d-%Y @ %Hh_%Mm_%Ss")
#Script exited gracefully (unless something else inside failed)
printf "%s %s" "abl_mass_qsub_sclean.sh has completed ${2}" "${global_end_time}" | mail -s "abl_mass_qsub complete" nvx4@cdc.gov
exit 0
