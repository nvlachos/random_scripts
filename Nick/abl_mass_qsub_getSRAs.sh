#!/bin/sh -l

#$ -o ablmq_getSRAs.out
#$ -e ablmq_getSRAs.err
#$ -N ablmq_getSRAs
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
# Usage ./abl_mass_qsub_getSRAs.sh path_to_list max_concurrent_submissions output_directory_for_scripts clobberness[keep|clobber]
#

# Number regex to test max concurrent submission parametr
number='^[0-9]+$'

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./abl_mass_qsub_getSRAs.sh path_to_list_file(single sample ID per line, e.g. B8VHY/1700128 (it must include project id also)) max_concurrent_submissions path_to_alt_database output_directory_for_scripts clobberness[keep|clobber]"
	exit 1
elif [[ ! -f "${1}" ]]; then
	echo "${1} (list) does not exist...exiting"
	exit 1
elif ! [[ ${2} =~ $number ]] || [[ -z "${2}" ]]; then
	echo "${2} is not a number or is empty. Please input max number of concurrent qsub submissions...exiting"
	exit 2
elif [[ -z "${3}" ]]; then
	echo "Output directory parameter is empty...exiting"
	exit 1
elif [[ -z "${4}" ]]; then
	echo "Clobberness was not input, be sure to add keep or clobber as 4th parameter...exiting"
	exit 1
fi

# Check that clobberness is a valid option
if [[ "${4}" != "keep" ]] && [[ "${4}" != "clobber" ]]; then
	echo "Clobberness was not input, be sure to add keep or clobber as 5th parameter...exiting"
	exit 1
else
	clobberness="${4}"
fi

# create an array of all samples in the list
arr=()
while IFS= read -r line || [ "$line" ];  do
  arr+=("$line")
done < ${1}

arr_size="${#arr[@]}"
last_index=$(( arr_size -1 ))
echo "-${arr_size}:${arr[@]}-"

# Create counter and set max number of concurrent submissions
counter=0
max_subs=${2}

# Set script directory
main_dir="${3}/getSRAs_subs"
if [[ ! -d "${3}/getSRAs_subs" ]]; then
	mkdir "${3}/getSRAs_subs"
	mkdir "${3}/getSRAs_subs/complete"
elif [[ ! -d  "${3}/getSRAs_subs/complete" ]]; then
	mkdir "${3}/getSRAs_subs/complete"
fi

start_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")

# Creates and submits qsub scripts to check all isolates on the list against the newest ResGANNCBI DB
while [ ${counter} -lt ${arr_size} ] ; do
	sample=$(echo "${arr[${counter}]}" | cut -d'/' -f2)
	project=$(echo "${arr[${counter}]}" | cut -d'/' -f1)
	if [[ "${clobberness}" = "clobber" ]]; then
		rm ${processed}/${project}/${sample}
	fi
	echo ${counter}
	# Check if counter is below max number of concurrent submissions
	if [ ${counter} -lt ${max_subs} ]; then
		# Check if either one of the concatted output files of trimming reads files exist, skip submission if so
		if [[ ! -f "${processed}/${project}/${sample}/trimmed/${sample}.paired.fq" ]]; then
			echo  "Index is below max submissions, submitting"
			echo -e "#!/bin/bash -l\n" > "${main_dir}/getSRAsAR_${sample}_${start_time}.sh"
			echo -e "#$ -o getSRAs_${sample}.out" >> "${main_dir}/getSRAs_${sample}_${start_time}.sh"
			echo -e "#$ -e getSRAs_${sample}.err" >> "${main_dir}/getSRAs_${sample}_${start_time}.sh"
			echo -e "#$ -N getSRAs_${sample}"   >> "${main_dir}/getSRAs_${sample}_${start_time}.sh"
			echo -e "#$ -cwd"  >> "${main_dir}/getSRAs_${sample}_${start_time}.sh"
			echo -e "#$ -q short.q\n"  >> "${main_dir}/getSRAs_${sample}_${start_time}.sh"
			# Can we somehow consolidate into one srst2 analysis to do MLST/AR/SEROTYPE
			echo -e "cd ${shareScript}" >> "${main_dir}/getSRAs_${sample}_${start_time}.sh"
			echo -e "\"${shareScript}/get_SRA_reads.sh\" \"${sample}\" \"${project}\"" >> "${main_dir}/getSRAs_${sample}_${start_time}.sh"
			echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_getSRAs_complete.txt\"" >> "${main_dir}/getSRAs_${sample}_${start_time}.sh"

			cd "${main_dir}"
			qsub "${main_dir}/getSRAs_${sample}_${start_time}.sh"
		# Old data existed, skipping
		else
			echo -e $(date) > "${main_dir}/complete/${sample}_getSRAs_complete.txt"
			echo "${project}/${sample} already has trimmed paired.fq"
		fi
	# Counter is above max submission, must wait for previous ones to finish before moving on
	else
		waiting_for_index=$(( counter - max_subs ))
		waiting_sample=$(echo "${arr[${waiting_for_index}]}" | cut -d'/' -f2)
		timer=0
		echo "Index is above max submissions, waiting for index ${waiting_for_index}:${waiting_sample} to complete"
		while :
		do
			# Check if timer is above max time allowed
			if [[ ${timer} -gt 1800 ]]; then
				echo "Timer exceeded limit of 1800 seconds 30 minutes"
				break
			fi
			# Check if waiting sample is finished
			if [ -f "${main_dir}/complete/${waiting_sample}_getSRAs_complete.txt" ]; then
				# Check if current sample has etiher one of the output files from srst2, skip analysis if so
				if [[ ! -f "${processed}/${project}/${sample}/trimmed/${sample}.paired.fq" ]]; then
					echo  "Index is below max submissions, submitting"
					echo -e "#!/bin/bash -l\n" > "${main_dir}/getSRAsAR_${sample}_${start_time}.sh"
					echo -e "#$ -o getSRAs_${sample}.out" >> "${main_dir}/getSRAs_${sample}_${start_time}.sh"
					echo -e "#$ -e getSRAs_${sample}.err" >> "${main_dir}/getSRAs_${sample}_${start_time}.sh"
					echo -e "#$ -N getSRAs_${sample}"   >> "${main_dir}/getSRAs_${sample}_${start_time}.sh"
					echo -e "#$ -cwd"  >> "${main_dir}/getSRAs_${sample}_${start_time}.sh"
					echo -e "#$ -q short.q\n"  >> "${main_dir}/getSRAs_${sample}_${start_time}.sh"
					# Can we somehow consolidate into one srst2 analysis to do MLST/AR/SEROTYPE
					echo -e "cd ${shareScript}" >> "${main_dir}/getSRAs_${sample}_${start_time}.sh"
					echo -e "\"${shareScript}/get_SRA_reads.sh\" \"${sample}\" \"${project}\"" >> "${main_dir}/getSRAs_${sample}_${start_time}.sh"
					echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_getSRAs_complete.txt\"" >> "${main_dir}/getSRAs_${sample}_${start_time}.sh"

					cd "${main_dir}"
					qsub "${main_dir}/getSRAs_${sample}_${start_time}.sh"
				# Old data existed, skipping
				else
					echo -e $(date) > "${main_dir}/complete/${sample}_getSRAs_complete.txt"
					echo "${project}/${sample} already has trimmed paired.fq"
				fi
				break
			# Wait 5 seconds and then check if "waiting" sample is complete
			else
				timer=$(( timer + 5 ))
				echo "sleeping for 5 seconds, so far slept for ${timer}"
				sleep 5
			fi
		done
	fi
	counter=$(( counter + 1 ))
done

# Check for completion of all samples
timer=0
for item in "${arr[@]}"; do
	waiting_sample=$(echo "${item}" | cut -d'/' -f2)
	if [[ -f "${main_dir}/complete/${waiting_sample}_getSRAs_complete.txt" ]]; then
		echo "${item} is complete"
		if [[ -f "${shareScript}/getSRAs_${waiting_sample}.out" ]]; then
			mv "${shareScript}/getSRAs_${waiting_sample}.out" "${main_dir}"
		fi
		if [[ -f "${shareScript}/getSRAs_${waiting_sample}.err" ]]; then
			mv "${shareScript}/getSRAs_${waiting_sample}.err" "${main_dir}"
		fi
	else
		while :
		do
				if [[ ${timer} -gt 3600 ]]; then
					echo "Timer exceeded limit of 3600 seconds = 60 minutes"
					exit 1
				fi
				if [[ -f "${main_dir}/complete/${waiting_sample}_getSRAs_complete.txt" ]]; then
					echo "${item} is complete"
					if [[ -f "${shareScript}/getSRAs_${waiting_sample}.out" ]]; then
						mv "${shareScript}/getSRAs_${waiting_sample}.out" "${main_dir}"
					fi
					if [[ -f "${shareScript}/getSRAs_${waiting_sample}.err" ]]; then
						mv "${shareScript}/getSRAs_${waiting_sample}.err" "${main_dir}"
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

exit 0
