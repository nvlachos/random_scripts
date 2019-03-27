#!/bin/sh -l

#$ -o ablmq-altm.out
#$ -e ablmq-altm.err
#$ -N ablmq-altm
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi
. ./config

#Import the module file that loads all necessary mods
. "${mod_changers}/pipeline_mods"

#List all currently loaded modules
#. ./module_changers/list_modules.sh

#
# Usage ./abl_mass_qsub_MLST.sh path_to_list max_concurrent_submissions  alt_mlst_DB output_directory_for_scripts clobberness[keep/clobber]
#

# Number regex to test max concurrent submission parametr
number='^[0-9]+$'

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to ./abl_mass_qsub_MLST.sh, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./abl_mass_qsub_MLST.sh path_to_list_file(single sample ID per line, e.g. B8VHY/1700128 (it must include project id also)) max_concurrent_submissions alt_databse output_directory_for_scripts clobberness[keep|clobber]"
	exit 1
elif [[ ! -f "${1}" ]]; then
	echo "${1} (list) does not exist...exiting"
	exit 1
elif ! [[ ${2} =~ $number ]] || [[ -z "${2}" ]]; then
	echo "${2} is not a number or is empty. Please input max number of concurrent qsub submissions...exiting"
	exit 2
elif [[ -z "${3}" ]]; then
	echo "alt_db is empty, exiting...."
elif [[ -z "${4}" ]]; then
	echo "No output directory for scripts given, exiting..."
elif [[ -z "${5}" ]]; then
	echo "Clobberness is empty, exiting...."
fi

# Check that clobberness is a valid option
if [[ "${5}" != "keep" ]] && [[ "${5}" != "clobber" ]]; then
	echo "Clobberness was not input correctly [keep|clobber]...exiting"
	exit 1
else
	clobberness="${5}"
fi

# create an array of all samples in the list
arr=()
while IFS= read -r line || [[ "$line" ]];  do
  arr+=("$line")
done < ${1}

arr_size="${#arr[@]}"
last_index=$(( arr_size -1 ))
echo "-${arr_size}:${arr[@]}-"

# Sets location of alt MLST database
mlst_dbs=$(mlst -list)
echo "${mlst_dbs}"
exit

# Create counter and set max concurrent submissions
counter=0
max_subs=${2}
alt_DB="${3}"

# Create direcory to hold all temporary qsub scripts
main_dir="${4}/mlst_subs"
if [[ ! -d "${4}/mlst_subs" ]]; then
	mkdir "${4}/mlst_subs"
	mkdir "${4}/mlst_subs/complete"
elif [[ ! -d "${4}/mlst_subs/complete" ]]; then
	mkdir "${4}/mlst_subs/complete"
fi

start_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")

# Creates and submits qsub script to determine MlST of all samples on the list
while [ ${counter} -lt ${arr_size} ] ; do
	sample=$(echo "${arr[${counter}]}" | cut -d'/' -f2)
	project=$(echo "${arr[${counter}]}" | cut -d'/' -f1)
	# Removes old data if clobbering is set
	if [[ "${clobberness}" = "clobber'" ]]; then
		rm "${processed}/${project}/${sample}/MLST/${sample}_${alt_DB}.mlst"
	fi
	# Check to see if the sample has usable assembly data
	if [[ -s "${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed.fasta" ]]; then
		# Check if current submission is less than max
		if [[ ${counter} -lt ${max_subs} ]]; then
			# Check if precious result exists and skip if so
			if [[ ! -f "${processed}/${project}/${sample}/MLST/${sample}_${alt_DB}.mlst" ]]; then
				echo  "Index is below max submissions, submitting"
				echo -e "#!/bin/bash -l\n" > "${main_dir}/alt_mlst_${sample}_${start_time}.sh"
				echo -e "#$ -o alt_mlst_${sample}.out" >> "${main_dir}/alt_mlst_${sample}_${start_time}.sh"
				echo -e "#$ -e alt_mlst_${sample}.err" >> "${main_dir}/alt_mlst_${sample}_${start_time}.sh"
				echo -e "#$ -N alt_mlst_${sample}"   >> "${main_dir}/alt_mlst_${sample}_${start_time}.sh"
				echo -e "#$ -cwd"  >> "${main_dir}/alt_mlst_${sample}_${start_time}.sh"
				echo -e "#$ -q short.q\n"  >> "${main_dir}/alt_mlst_${sample}_${start_time}.sh"
				echo -e "\"${shareScript}/run_MLST.sh\" \"${sample}\" \"${project}\" \"-f\" \"${alt_DB}\"" >> "${main_dir}/alt_mlst_${sample}_${start_time}.sh"
				echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_alt_mlst_complete.txt\"" >> "${main_dir}/alt_mlst_${sample}_${start_time}.sh"
				cd "${main_dir}"
				if [[ "${counter}" -lt "${last_index}" ]]; then
					qsub "${main_dir}/alt_mlst_${sample}_${start_time}.sh"
				else
					qsub -sync y "${main_dir}/alt_mlst_${sample}_${start_time}.sh"
				fi
			else
				echo "${project}/${sample} already has mlst summary"
				echo "$(date)" > "${main_dir}/complete/${sample}_alt_mlst_complete.txt"
			fi
			# If counter is above maximum concurrent submissions
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
				# Check if waiting sample has finished
				if [[ -f "${main_dir}/complete/${waiting_sample}_alt_mlst_complete.txt" ]]; then
					# Check if precious data exists and skip if so
					if [[ ! -f "${processed}/${project}/${sample}/MLST/${sample}_${alt_DB}.mlst" ]]; then
						echo  "Index is below max submissions, submitting"
						echo -e "#!/bin/bash -l\n" > "${main_dir}/alt_mlst_${sample}_${start_time}.sh"
						echo -e "#$ -o alt_mlst_${sample}.out" >> "${main_dir}/alt_mlst_${sample}_${start_time}.sh"
						echo -e "#$ -e alt_mlst_${sample}.err" >> "${main_dir}/alt_mlst_${sample}_${start_time}.sh"
						echo -e "#$ -N alt_mlst_${sample}"   >> "${main_dir}/alt_mlst_${sample}_${start_time}.sh"
						echo -e "#$ -cwd"  >> "${main_dir}/alt_mlst_${sample}_${start_time}.sh"
						echo -e "#$ -q short.q\n"  >> "${main_dir}/alt_mlst_${sample}_${start_time}.sh"
						echo -e "\"${shareScript}/run_MLST.sh\" \"${sample}\" \"${project}\" \"-f\" \"${alt_DB}\"" >> "${main_dir}/alt_mlst_${sample}_${start_time}.sh"
						echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_alt_mlst_complete.txt\"" >> "${main_dir}/alt_mlst_${sample}_${start_time}.sh"
						cd "${main_dir}"
						if [[ "${counter}" -lt "${last_index}" ]]; then
							qsub "${main_dir}/alt_mlst_${sample}_${start_time}.sh"
						else
							qsub -sync y "${main_dir}/alt_mlst_${sample}_${start_time}.sh"
						fi
					# Previous data exists, skip
					else
						echo "${project}/${sample} already has mlst summary"
						echo "$(date)" > "${main_dir}/complete/${sample}_alt_mlst_complete.txt"
					fi
					break
				# waiting sample has not completed, then wait 5 seconds and try again
				else
					timer=$(( timer + 5 ))
					echo "sleeping for 5 seconds, so far slept for ${timer}"
					sleep 5
				fi
			done
		fi
	# Sample does not have a usable assembly file
	else
		echo "${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed.fasta not found"
		echo "$(date)" > "${main_dir}/complete/${sample}_alt_mlst_complete.txt"
	fi
	counter=$(( counter + 1 ))
done

# Loop to ensure all samples are complete (or time runs) before allowing the script to exit
timer=0
for item in "${arr[@]}"; do
	waiting_sample=$(echo "${item}" | cut -d'/' -f2)
	if [[ -f "${main_dir}/complete/${waiting_sample}_csstarn_complete.txt" ]] || [[ ! -s "${processed}/${project}/${waiting_sample}/Assembly/${waiting_sample}_scaffolds_trimmed.fasta" ]]; then
		echo "${item} is complete normal"
	else
		# Check every 5 seconds to see if the sample has completed normal csstar analysis
		while :
		do
				if [[ ${timer} -gt 3600 ]]; then
					echo "Timer exceeded limit of 3600 seconds = 60 minutes"
					exit 1
				fi
				if [[ -f "${main_dir}/complete/${waiting_sample}_csstarn_complete.txt" ]]; then
					echo "${item} is complete"
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
exit 0
