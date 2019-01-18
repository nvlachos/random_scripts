#!/bin/sh -l

#$ -o ablmq-cs.out
#$ -e ablmq-cs.err
#$ -N ablmq-cs
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
. ./config.sh
#Import the module file that loads all necessary mods
. "${mod_changers}/pipeline_mods"

#
# Usage ./abl_mass_qsub_MLST.sh path_to_list max_concurrent_submissions
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to ./abl_mass_qsub_MLST.sh, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./abl_mass_qsub_MLST.sh path_to_list_file(single sample ID per line, e.g. B8VHY/1700128 (it must include project id also)) max_concurrent_submissions alt_databse(optional)"
	exit 1
elif [[ ! -f "${1}" ]]; then
	echo "${1} (list) does not exist...exiting"
	exit 1
fi

# create an array of all samples in the list
arr=()
while IFS= read -r line || [[ "$line" ]];  do
  arr+=("$line")
done < ${1}

arr_size="${#arr[@]}"
last_index=$(( arr_size -1 ))
echo "-${arr_size}:${arr[@]}-"

alt_DB=""
if [[ ! -z ${3} ]]; then
	alt_DB="${3}"
fi

# Create direcory to hold all temporary qsub scripts
counter=0
max_subs=${2}

# Set script directory
main_dir="${share}/mass_subs/mlst_subs"
if [[ ! -d "${share}/mass_subs/mlst_subs" ]]; then
	mkdir "${share}/mass_subs/mlst_subs"
	mkdir "${share}/mass_subs/mlst_subs/complete"
elif [[ ! -d "${share}/mass_subs/mlst_subs/complete" ]]; then
	mkdir "${share}/mass_subs/mlst_subs/complete"
fi

start_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")

# Creates and submits qsub script to determine MlST of all samples on the list
while [ ${counter} -lt ${arr_size} ] ; do
	sample=$(echo "${arr[${counter}]}" | cut -d'/' -f2)
	project=$(echo "${arr[${counter}]}" | cut -d'/' -f1)
	if [[ -s "${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed.fasta" ]]; then
		if [[ ${counter} -lt ${max_subs} ]]; then
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
				if [[ "${counter}" -lt "${last_index}" ]]; then
					qsub "${main_dir}/alt_mlst_${sample}_${start_time}.sh"
				else
					qsub -sync y "${main_dir}/alt_mlst_${sample}_${start_time}.sh"
				fi
			else
				echo "${project}/${sample} already has mlst summary"
				echo "$(date)" > "${main_dir}/complete/${sample}_alt_mlst_complete.txt"
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
				if [[ -f "${main_dir}/complete/${waiting_sample}_alt_mlst_complete.txt" ]]; then
					if [[ ! -f "${processed}/${project}/${sample}/MLST/${sample}.mlst" ]]; then
						echo  "Index is below max submissions, submitting"
						echo -e "#!/bin/bash -l\n" > "${main_dir}/alt_mlst_${sample}_${start_time}.sh"
						echo -e "#$ -o alt_mlst_${sample}.out" >> "${main_dir}/alt_mlst_${sample}_${start_time}.sh"
						echo -e "#$ -e alt_mlst_${sample}.err" >> "${main_dir}/alt_mlst_${sample}_${start_time}.sh"
						echo -e "#$ -N alt_mlst_${sample}"   >> "${main_dir}/alt_mlst_${sample}_${start_time}.sh"
						echo -e "#$ -cwd"  >> "${main_dir}/alt_mlst_${sample}_${start_time}.sh"
						echo -e "#$ -q short.q\n"  >> "${main_dir}/alt_mlst_${sample}_${start_time}.sh"
						echo -e "\"${shareScript}/run_MLST.sh\" \"${sample}\" \"${project}\" \"-f\" \"${alt_DB}\"" >> "${main_dir}/alt_mlst_${sample}_${start_time}.sh"
						echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_alt_mlst_complete.txt\"" >> "${main_dir}/alt_mlst_${sample}_${start_time}.sh"
						if [[ "${counter}" -lt "${last_index}" ]]; then
							qsub "${main_dir}/alt_mlst_${sample}_${start_time}.sh"
						else
							qsub -sync y "${main_dir}/alt_mlst_${sample}_${start_time}.sh"
						fi
					else
						echo "${project}/${sample} already has mlst summary"
						echo "$(date)" > "${main_dir}/complete/${sample}_alt_mlst_complete.txt"
					fi
					break
				else
					timer=$(( timer + 5 ))
					echo "sleeping for 5 seconds, so far slept for ${timer}"
					sleep 5
				fi
			done
		fi
	else
		echo "${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed.fasta not found"
		echo "$(date)" > "${main_dir}/complete/${sample}_alt_mlst_complete.txt"
	fi
	counter=$(( counter + 1 ))
done

echo "All isolates completed"
exit 0
