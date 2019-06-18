#!/bin/sh -l

#$ -o abl_ID.out
#$ -e abl_ID.err
#$ -N ablID
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
. ./config.sh
#Import the module file that loads all necessary mods
. "${mod_changers}/pipeline_mods"

#List all currently loaded modules
#. ./module_changers/list_modules.sh

#
# Usage ./abl_masss_qsub_det_ID.sh path_to_list max_concurrent_submission output_directory_for_scripts clobberness[keep|clobber]
#

# Number regex to test max concurrent submission parametr
number='^[0-9]+$'

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./abl_mass_qsub_det_ID.sh path_to_list_file(single sample ID per line, e.g. B8VHY/1700128 (it must include project id also)) max_submissions output_directory_for_scripts clobberness[keep|clobber]"
	echo "Output location varies depending on which tasks are performed but will be found somewhere under ${processed}"
	exit 0
elif [[ ! -f "${1}" ]]; then
	echo "${1} (list) does not exist...exiting"
	exit 1
elif ! [[ ${2} =~ $number ]] || [[ -z "${2}" ]]; then
	echo "${2} is not a number or is empty. Please input max number of concurrent qsub submissions...exiting"
	exit 2
elif [[ -z "${3}" ]]; then
	echo "No script output directory given...exiting"
	exit 3
elif [[ -z "${4}" ]]; then
	echo "Clobberness was not input, be sure to add keep or clobber as 4th parameter...exiting"
	exit 4
fi

# Check that clobberness is a valid option
if [[ "${4}" != "keep" ]] && [[ "${4}" != "clobber" ]]; then
	echo "Clobberness was not input correctly, be sure to add keep or clobber as 5th parameter...exiting"
	exit 1
else
	clobberness="${4}"
fi

# Creates an array of all samples in the provided list
arr=()
while IFS= read -r line || [ "$line" ];  do
  arr+=("$line")
done < ${1}

arr_size="${#arr[@]}"
last_index=$(( arr_size -1 ))
echo "-${arr_size}:${arr[@]}-"

# Counter to track current isolate being worked on in list
counter=0
max_subs="${2}"

# Creates output directory for scripts if it does not exist yet
main_dir="${3}/detID_subs"
if [[ ! -d "${3}/detID_subs" ]]; then
	mkdir "${3}/detID_subs"
	mkdir "${3}/detID_subs/complete"
elif [[ ! -d  "${3}/detID_subs/complete" ]]; then
	mkdir "${3}/detID_subs/complete"
fi

# Go through whole list individually to analyze samples
while [ ${counter} -lt ${arr_size} ] ; do
	sample=$(echo "${arr[${counter}]}" | cut -d'/' -f2)
	project=$(echo "${arr[${counter}]}" | cut -d'/' -f1)
	if [[ "${clobberness}" == "clobber" ]]; then
		rm ${processed}/${project}/${sample}/${sample}.tax
	fi
	echo ${counter}
	# Check if counter is below max concurrent submissions limit
	if [ ${counter} -lt ${max_subs} ]; then
		# Check if old data exists and skip if so
		if [[ ! -s "${processed}/${project}/${sample}/${sample}.tax" ]]; then
			echo  "Index is below max submissions, submitting"
			echo -e "#!/bin/bash -l\n" > "${main_dir}/detID_${sample}_${start_time}.sh"
			echo -e "#$ -o detID_${sample}.out" >> "${main_dir}/detID_${sample}_${start_time}.sh"
			echo -e "#$ -e detID_${sample}.err" >> "${main_dir}/detID_${sample}_${start_time}.sh"
			echo -e "#$ -N detID_${sample}"   >> "${main_dir}/detID_${sample}_${start_time}.sh"
			echo -e "#$ -cwd"  >> "${main_dir}/detID_${sample}_${start_time}.sh"
			echo -e "#$ -q short.q\n"  >> "${main_dir}/detID_${sample}_${start_time}.sh"
			echo -e "cd ${shareScript}" >> "${main_dir}/detID_${sample}_${start_time}.sh"
			echo -e "${shareScript}/determine_taxID.sh ${sample} ${project}" >> "${main_dir}/detID_${sample}_${start_time}.sh"
			echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_detID_complete.txt\"" >> "${main_dir}/detID_${sample}_${start_time}.sh"
			cd "${main_dir}"
			#if [[ "${counter}" -lt "${last_index}" ]]; then
				qsub "${main_dir}/detID_${sample}_${start_time}.sh"
			#else
			#	qsub -sync y "${main_dir}/detID_${sample}_${start_time}.sh"
			#fi
		# Sample had old data, skipping
		else
			echo "${project}/${sample} already had its detIDs removed"
			echo "$(date)" > "${main_dir}/complete/${sample}_detID_complete.txt"
		fi
	# If counter is over max concurrent submissions limit
	else
		waiting_for_index=$(( counter - max_subs ))
		waiting_sample=$(echo "${arr[${waiting_for_index}]}" | cut -d'/' -f2)
		timer=0
		echo "Index is above max submissions, waiting for index ${waiting_for_index}:${waiting_sample} to complete"
		# Loop to check if "waiting" sample has completed
		while :
		do
			#	Check if max time limit has been reached, exit if so
			if [[ ${timer} -gt 1800 ]]; then
				echo "Timer exceeded limit of 1800 seconds 30 minutes"
				break
			fi
			# Check if "waiting" sample is complete, else wait longer
			if [ -f "${main_dir}/complete/${waiting_sample}_detID_complete.txt" ]; then
				# Check for old data, skip if present
				if [[ ! -s "${processed}/${project}/${sample}/${sample}.tax" ]]; then
					echo  "Index is below max submissions, submitting"
					echo -e "#!/bin/bash -l\n" > "${main_dir}/detID_${sample}_${start_time}.sh"
					echo -e "#$ -o detID_${sample}.out" >> "${main_dir}/detID_${sample}_${start_time}.sh"
					echo -e "#$ -e detID_${sample}.err" >> "${main_dir}/detID_${sample}_${start_time}.sh"
					echo -e "#$ -N detID_${sample}"   >> "${main_dir}/detID_${sample}_${start_time}.sh"
					echo -e "#$ -cwd"  >> "${main_dir}/detID_${sample}_${start_time}.sh"
					echo -e "#$ -q short.q\n"  >> "${main_dir}/detID_${sample}_${start_time}.sh"
					echo -e "cd ${shareScript}" >> "${main_dir}/detID_${sample}_${start_time}.sh"
					echo -e "${shareScript}/determine_taxID.sh ${sample} ${project}" >> "${main_dir}/detID_${sample}_${start_time}.sh"
					echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_detID_complete.txt\"" >> "${main_dir}/detID_${sample}_${start_time}.sh"
					cd "${main_dir}"
					#if [[ "${counter}" -lt "${last_index}" ]]; then
						qsub "${main_dir}/detID_${sample}_${start_time}.sh"
					#else
					#	qsub -sync y "${main_dir}/detID_${sample}_${start_time}.sh"
					#fi
					break
				# Has old data, skipping
				else
					echo "${project}/${sample} already had its detIDs removed"
					echo "$(date)" > "${main_dir}/complete/${sample}_detID_complete.txt"
					break
				fi
			fi
		done
	fi
	counter=$(( counter + 1 ))
done

# Check for completion of all samples
timer=0
for item in "${arr[@]}"; do
	waiting_sample=$(echo "${item}" | cut -d'/' -f2)
	if [[ -f "${main_dir}/complete/${waiting_sample}_detID_complete.txt" ]]; then
		echo "${item} is complete"
		if [[ -f "${shareScript}/detID_${waiting_sample}.out" ]]; then
			mv "${shareScript}/detID_${waiting_sample}.out" "${main_dir}"
		fi
		if [[ -f "${shareScript}/detID_${waiting_sample}.err" ]]; then
			mv "${shareScript}/detID_${waiting_sample}.err" "${main_dir}"
		fi
	else
		while :
		do
				if [[ ${timer} -gt 600 ]]; then
					echo "Timer exceeded limit of 600 seconds = 10 minutes"
					exit 1
				fi
				if [[ -f "${main_dir}/complete/${waiting_sample}_detID_complete.txt" ]]; then
					echo "${item} is complete"
					if [[ -f "${shareScript}/detID_${waiting_sample}.out" ]]; then
						mv "${shareScript}/detID_${waiting_sample}.out" "${main_dir}"
					fi
					if [[ -f "${shareScript}/detID_${waiting_sample}.err" ]]; then
						mv "${shareScript}/detID_${waiting_sample}.err" "${main_dir}"
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
printf "%s %s" "abl_mass_qsub_det_ID.sh has completed" "${global_end_time}" | mail -s "abl_mass_qsub_det_ID.sh complete" nvx4@cdc.gov

#Script exited gracefully (unless something else inside failed)
exit 0
