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
# Usage ./abl_mass_qsub_ANI.sh path_to_list max_concurrent_submissions
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to ./abl_mass_qsub_ANI.sh, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./abl_mass_qsub_ANI.sh path_to_list_file(single sample ID per line, e.g. B8VHY/1700128 (it must include project id also)) max_concurrent_submissions"
	exit 1
elif [[ ! -f "${1}" ]]; then
	echo "${1} (list) does not exist...exiting"
	exit 1
fi

# Create an array of all samples in the list
arr=()
while IFS= read -r line || [[ "$line" ]];  do
  arr+=("$line")
done < ${1}

arr_size="${#arr[@]}"
last_index=$(( arr_size -1 ))
echo "-${arr_size}:${arr[@]}-"

# Create counter and set max number of concurrent submissions
counter=0
max_subs=${2}

# Set script directory
main_dir="${share}/mass_subs/ANI_subs"
if [[ ! -d "${share}/mass_subs/ANI_subs" ]]; then
	mkdir "${share}/mass_subs/ANI_subs"
	mkdir "${share}/mass_subs/ANI_subs/complete"
elif [[ ! -d "${share}/mass_subs/ANI_subs/complete" ]]; then
	mkdir "${share}/mass_subs/ANI_subs/complete"
fi

start_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")

# Create and submit qsub scripts to get ANI for all isolates
while [ ${counter} -lt ${arr_size} ] ; do
	sample=$(echo "${arr[${counter}]}" | cut -d'/' -f2)
	short_sample=${sample:0:20}
	project=$(echo "${arr[${counter}]}" | cut -d'/' -f1)
	echo "${sample} and ${project}:"
	echo "${counter}-${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed.fasta"
	# If sample has assembly, then delete old ANI folder to allow rerun
	if [[ -s "${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed.fasta" ]]; then
		mv -r "${processed}/${project}/${sample}/ANI_original/"
	fi
	if [[ -s "${processed}/${project}/${sample}/${sample}.tax" ]]; then
		while IFS= read -r line;
		do
			# Grab first letter of line (indicating taxonomic level)
			first=${line:0:1}
			# Assign taxonomic level value from 4th value in line (1st-classification level,2nd-% by kraken, 3rd-true % of total reads, 4th-identifier)
			if [ "${first}" = "s" ]
			then
				species=$(echo "${line}" | awk -F ' ' '{print $2}')
			elif [ "${first}" = "G" ]
			then
				genus=$(echo "${line}" | awk -F ' ' '{print $2}')
			fi
		done < "${processed}/${project}/${sample}/${sample}.tax"
		if [[ ! -f "${processed}/${sample}/ANI/best_ANI_hits_ordered(${sample}_vs_${genus,})" ]]; then
			rm -r "${processed}/${sample}/ANI/best_ANI_hits_ordered(${sample}_vs_${genus,})"
		fi
		if [[ ! -d "${processed}/${sample}/ANI/aniM" ]]; then
			rm -r "${processed}/${sample}/ANI/aniM"
		fi
	 	if [[ ${counter} -lt ${max_subs} ]]; then
			if [[ ! -f "${processed}/${sample}/ANI/best_ANI_hits_ordered(${sample}_vs_${genus,})" ]]; then
				echo  "Index is below max submissions, submitting"
				echo "Going to make ${main_dir}/ani_${short_sample}_${start_time}.sh"
				echo -e "#!/bin/bash -l\n" > "${main_dir}/ani_${short_sample}_${start_time}.sh"
				echo -e "#$ -o ani_${sample}.out" >> "${main_dir}/ani_${short_sample}_${start_time}.sh"
				echo -e "#$ -e ani_${sample}.err" >> "${main_dir}/ani_${short_sample}_${start_time}.sh"
				echo -e "#$ -N ani_${sample}"   >> "${main_dir}/ani_${short_sample}_${start_time}.sh"
				echo -e "#$ -cwd"  >> "${main_dir}/ani_${short_sample}_${start_time}.sh"
				echo -e "#$ -q short.q\n"  >> "${main_dir}/ani_${short_sample}_${start_time}.sh"
				echo -e "\"${shareScript}/run_ANI.sh\" \"${sample}\" \"${genus}\" \"${species}\" \"${project}\"" >> "${main_dir}/ani_${short_sample}_${start_time}.sh"
				echo -e "echo \"$(date)\" > \"${main_dir}/complete/${short_sample}_ani_complete.txt\"" >> "${main_dir}/ani_${short_sample}_${start_time}.sh"
				if [[ "${counter}" -lt "${last_index}" ]]; then
					qsub "${main_dir}/ani_${short_sample}_${start_time}.sh"
				else
					qsub -sync y "${main_dir}/ani_${short_sample}_${start_time}.sh"
				fi
			else
				echo "${project}/${sample} already has ANI summary"
				echo "$(date)" > "${main_dir}/complete/${sample}_ani_complete.txt"
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
				if [[ -f "${main_dir}/complete/${waiting_sample}_ani_complete.txt" ]]; then
					if [[ ! -f "${processed}/${project}/${sample}/ANI/best_ANI_hits_ordered(${sample}_vs_${genus,})" ]]; then
						echo  "${waiting_sample}(${waiting_for_index}) is not complete, submitting ${sample} ($counter)"
						echo -e "#!/bin/bash -l\n" > "${main_dir}/ani_${sample}_${start_time}.sh"
						echo -e "#$ -o ani_${sample}.out" >> "${main_dir}/ani_${sample}_${start_time}.sh"
						echo -e "#$ -e ani_${sample}.err" >> "${main_dir}/ani_${sample}_${start_time}.sh"
						echo -e "#$ -N ani_${sample}"   >> "${main_dir}/ani_${sample}_${start_time}.sh"
						echo -e "#$ -cwd"  >> "${main_dir}/ani_${sample}_${start_time}.sh"
						echo -e "#$ -q short.q\n"  >> "${main_dir}/ani_${sample}_${start_time}.sh"
						echo -e "\"${shareScript}/run_ANI.sh\" \"${sample}\" \"${genus}\" \"${species}\" \"${project}\"" >> "${main_dir}/ani_${sample}_${start_time}.sh"
						echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_ani_complete.txt\"" >> "${main_dir}/ani_${sample}_${start_time}.sh"
						if [[ "${counter}" -lt "${last_index}" ]]; then
							qsub "${main_dir}/ani_${short_sample}_${start_time}.sh"
						else
							qsub -sync y "${main_dir}/ani_${short_sample}_${start_time}.sh"
						fi
					else
						echo "${project}/${sample} already has ANI summary"
						echo "$(date)" > "${main_dir}/complete/${sample}_ani_complete.txt"
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
	fi
	counter=$(( counter + 1 ))
done

echo "All isolates completed"
exit 0
