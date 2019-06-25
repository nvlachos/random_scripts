#!/bin/sh -l

#$ -o ablmq-16s.out
#$ -e ablmq-16s.err
#$ -N ablmq-16s
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi
. ./config.sh

#Import the module file that loads all necessary mods
. "${mod_changers}/pipeline_mods"

#List all currently loaded modules
#. ./module_changers/list_modules.sh

#
# Usage ./abl_mass_qsub_16s.sh path_to_list max_concurrent_submission output_directory_for_scripts clobberness[clobber/keep]
#


# Number regex to test max concurrent submission parametr
number='^[0-9]+$'

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./abl_mass_qsub_16s.sh path_to_list_file(single sample ID per line, e.g. B8VHY/1700128 (it must include project id also) max_concurrent_submissions output_directory_for_scripts clobberness[keep|clobber]"
	exit 0
elif ! [[ ${2} =~ $number ]] || [[ -z "${2}" ]]; then
	echo "${2} is not a number or is empty. Please input max number of concurrent qsub submissions...exiting"
	exit 2
elif [[ -z "${4}" ]] ; then
	echo "clobberness is empty please choose clobber or keep...exiting"
	exit 1
fi

# Check that clobberness is a valid option
if [[ "${4}" != "keep" ]] && [[ "${4}" != "clobber" ]]; then
	echo "Clobberness was not input, be sure to add keep or clobber as 5th parameter...exiting"
	exit 1
else
	clobberness="${4}"
fi

# Creates an array of all samples in the provided list
arr=()
while IFS= read -r line || [ "$line" ];  do
  arr+=("$line")
done < ${1}

# Gets size and prints out all array elements
arr_size="${#arr[@]}"
echo "-${arr_size}:${arr[@]}-"

 # Set known variables
counter=0
max_subs=${2}

# Create all sub output folders
main_dir="${3}/blast16s_subs"
if [[ ! -d "${3}/blast16s_subs" ]]; then
	mkdir "${3}/blast16s_subs"
	mkdir "${3}/blast16s_subs/complete"
elif [[ ! -d "${3}/blast16s_subs/complete" ]]; then
	mkdir "${3}/blast16s_subs/complete"
fi

start_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")

# Loops through every item in the array. Will wait for one to finish is maximum submission number is hit
while [ ${counter} -lt ${arr_size} ] ; do
	sample=$(echo "${arr[${counter}]}" | cut -d'/' -f2)
	project=$(echo "${arr[${counter}]}" | cut -d'/' -f1)
	# Delete the old folder
	if [[ "${clobberness}" = "clobber" ]]; then
		rm -r "${processed}/${project}/${sample}/16s"
	fi
	# Check to see if there is an assembly file to look in before starting
	if [[ -s "${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed.fasta" ]]; then
		# Check if max submissions has been hit
		if [[ ${counter} -lt ${max_subs} ]]; then
			# Checks if there is already an output file for this sample
			if [[ ! -f "${processed}/${project}/${sample}/16s/${sample}_16s_blast_id.txt" ]]; then
				# Write all lines of temp script to file"
				echo  "Index is below max submissions, submitting"
				echo -e "#!/bin/bash -l\n" > "${main_dir}/blast16s_${sample}_${start_time}.sh"
				echo -e "#$ -o blast16s_${sample}.out" >> "${main_dir}/blast16s_${sample}_${start_time}.sh"
				echo -e "#$ -e blast16s_${sample}.err" >> "${main_dir}/blast16s_${sample}_${start_time}.sh"
				echo -e "#$ -N blast16s_${sample}" >> "${main_dir}/blast16s_${sample}_${start_time}.sh"
				echo -e "#$ -cwd"  >> "${main_dir}/blast16s_${sample}_${start_time}.sh"
				echo -e "#$ -q short.q\n" >> "${main_dir}/blast16s_${sample}_${start_time}.sh"
				echo -e "cd ${shareScript}" >> "${main_dir}/blast16s_${sample}_${start_time}.sh"
				echo -e "\"${shareScript}/16s_blast.sh\" \"-n\" \"${sample}\" \"-p\" \"${project}\"" >> "${main_dir}/blast16s_${sample}_${start_time}.sh"
				echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_blast16s_complete.txt\"" >> "${main_dir}/blast16s_${sample}_${start_time}.sh"
				# Moves to working dir of mass qsub scripts
				cd "${main_dir}"
				# Checks o see if it is last entry to sync it, not allowing script to end
				#if [[ "${counter}" -lt "${last_index}" ]]; then
					qsub "${main_dir}/blast16s_${sample}_${start_time}.sh"
				#else
				#	qsub -sync y "${main_dir}/blast16s_${sample}_${start_time}.sh"
				#fi
			# Case that data already exists for this sample, creates a completed txt notifier for if there are more samples than max submissions allowed
			else
				echo "${project}/${sample} already has 16s summary"
				echo "$(date)" > "${main_dir}/complete/${sample}_blast16s_complete.txt"
			fi
		# Once max submissions are submitted script will wait until the isolate - max submissions is complete
		else
			waiting_for_index=$(( counter - max_subs ))
			waiting_sample=$(echo "${arr[${waiting_for_index}]}" | cut -d'/' -f2)
			timer=0
			echo "Index is above max submissions, waiting for index ${waiting_for_index}:${waiting_sample} to complete"
			# Loops at most 10 mins while waiting for "waiting_sample" completes
			while :
			do
				# Exit loop due to time limit exceeded
				if [[ ${timer} -gt 1800 ]]; then
					echo "Timer exceeded limit of 1800 seconds 30 minutes"
					break
				fi
				# Check to see if the pwaiting sample is finished by checking for the text notification
				if [[ -f "${main_dir}/complete/${waiting_sample}_blast16s_complete.txt" ]]; then
					# Checks if there is already an output file for this sample
					if [[ ! -f "${processed}/${project}/${sample}/16s/${sample}_16s_blast_id.txt" ]]; then
						echo  "Index is below max submissions, submitting"
						echo -e "#!/bin/bash -l\n" > "${main_dir}/blast16s_${sample}_${start_time}.sh"
						echo -e "#$ -o blast16s_${sample}.out" >> "${main_dir}/blast16s_${sample}_${start_time}.sh"
						echo -e "#$ -e blast16s_${sample}.err" >> "${main_dir}/blast16s_${sample}_${start_time}.sh"
						echo -e "#$ -N blast16s_${sample}" >> "${main_dir}/blast16s_${sample}_${start_time}.sh"
						echo -e "#$ -cwd"  >> "${main_dir}/blast16s_${sample}_${start_time}.sh"
						echo -e "#$ -q short.q\n" >> "${main_dir}/blast16s_${sample}_${start_time}.sh"
						echo -e "cd ${shareScript}" >> "${main_dir}/blast16s_${sample}_${start_time}.sh"
						echo -e "\"${shareScript}/16s_blast.sh\" \"-n\" \"${sample}\" \"-p\" \"${project}\"" >> "${main_dir}/blast16s_${sample}_${start_time}.sh"
						echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_blast16s_complete.txt\"" >> "${main_dir}/blast16s_${sample}_${start_time}.sh"
						# Moves to working dir of mass qsub scripts
						cd "${main_dir}"
						qsub "${main_dir}/blast16s_${sample}_${start_time}.sh"
					# Case that data already exists for this sample, creates a completed txt notifier for if there are more samples than max submissions allowed
					else
						echo "${project}/${sample} already has 16s summary"
						echo "$(date)" > "${main_dir}/complete/${sample}_blast16s_complete.txt"
					fi
					break
				# If waiting sample is not complete yet, then wait 5 seconds and try again
				else
					timer=$(( timer + 5 ))
					echo "sleeping for 5 seconds, so far slept for ${timer}"
					sleep 5
				fi
			done
		fi
	# Skip isolate due to no assembly being available. Mark sample as being done so that if another sample is waiting it doesnt get hung up
	else
		echo "${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed.fasta not found"
		echo "$(date)" > "${main_dir}/complete/${sample}_blast16s_complete.txt"
	fi
	counter=$(( counter + 1 ))
done

# Do a final check that all samples are completed before exiting script
for item in "${arr[@]}"; do
	# Get current sample ID
	waiting_sample=$(echo "${item}" | cut -d'/' -f2)
	# Checks if the sample_complete file is existant or if not that an assembly file is also absent, indicating nothing can be done
	if [[ -f "${main_dir}/complete/${waiting_sample}_blast16s_complete.txt" ]] || [[ ! -s "${processed}/${project}/${waiting_sample}/Assembly/${waiting_sample}_scaffolds_trimmed.fasta" ]]; then
		echo "${item} is complete"
		if [[ -f "${shareScript}/blast16s_${sample}.out" ]]; then
			mv "${shareScript}/blast16s_${sample}.out" "${main_dir}"
		fi
		if [[ -f "${shareScript}/blast16s_${sample}.err" ]]; then
			mv "${shareScript}/blast16s_${sample}.err" "${main_dir}"
		fi
	# If sample is not complete wait until time limit
	else
		while :
		do
				# Check if time elapsed is greater than max time limit to wait for all isolates to finish
				if [[ ${timer} -gt 3600 ]]; then
					echo "Timer exceeded limit of 3600 seconds = 60 minutes"
					exit 1
				fi
				# Move to next sample if current isolate finishes
				if [[ -f "${main_dir}/complete/${waiting_sample}_blast16s_complete.txt" ]]; then
					echo "${item} is complete"
					if [[ -f "${shareScript}/blast16s_${sample}.out" ]]; then
						mv "${shareScript}/blast16s_${sample}.out" "${main_dir}"
					fi
					if [[ -f "${shareScript}/blast16s_${sample}.err" ]]; then
						mv "${shareScript}/blast16s_${sample}.err" "${main_dir}"
					fi
					break
				# Wait set time before checking for completion again
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
printf "%s %s" "abl_mass_qsub_16s.sh has completed" "${global_end_time}" | mail -s "abl_mass_qsub_16s.sh complete" nvx4@cdc.gov

#Script exited gracefully (unless something else inside failed)
exit 0
