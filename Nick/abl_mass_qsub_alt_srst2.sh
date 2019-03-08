#!/bin/sh -l

#$ -o ablmq_srst2_alt.out
#$ -e ablmq_srst2_alt.err
#$ -N ablmq_srst2_alt
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
. ./config.sh
#Import the module file that loads all necessary mods
. "${mod_changers}/pipeline_mods"

#
# Usage ./abl_mass_qsub_alt_srst2.sh path_to_list max_concurrent_submissions path_to_alt_database output_directory_for_scripts
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to ./abl_mass_qsub_alt_srst2.sh, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./abl_mass_qsub_alt_srst2.sh path_to_list_file(single sample ID per line, e.g. B8VHY/1700128 (it must include project id also)) max_concurrent_submissions path_to_alt_database output_directory_for_scripts"
	exit 1
elif [[ ! -f "${1}" ]]; then
	echo "${1} (list) does not exist...exiting"
	exit 1
elif [[ ! -f "${3}" ]]; then
	echo "${3} (alt_db) does not exist...exiting"
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

# Create counter and set max number of concurrent submissions
counter=0
max_subs=${2}

# Set script directory
main_dir="${4}/srst2_alt_subs"
if [[ ! -d "${4}/srst2_alt_subs" ]]; then
	mkdir "${4}/srst2_alt_subs"
	mkdir "${4}/srst2_alt_subs/complete"
elif [[ ! -d  "${4}/srst2_alt_subs/complete" ]]; then
	mkdir "${4}/srst2_alt_subs/complete"
fi

# Format name of DB used from filrname for srst2
alt_DB_path=${3}
alt_DB=$(echo ${alt_DB_path##*/} | cut -d'.' -f1)
alt_DB=${alt_DB//_srst2/}
echo ${alt_DB}
start_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")

# Create and submit srst2 qsub scripts for each isolate in the list
while [ ${counter} -lt ${arr_size} ] ; do
	sample=$(echo "${arr[${counter}]}" | cut -d'/' -f2)
	project=$(echo "${arr[${counter}]}" | cut -d'/' -f1)
	echo ${counter}
	if [ ${counter} -lt ${max_subs} ]; then
		# Fix any improperly named (older style) files
		if [[ -f "${processed}/${project}/${sample_name}/srst2/${sample_name}_${alt_DB}__genes__${alt_DB}_srst2__results.txt" ]]; then
				mv "${processed}/${project}/${sample_name}/srst2/${sample_name}_${alt_DB}__genes__${alt_DB}_srst2__results.txt" "${processed}/${project}/${sample_name}/srst2/${sample_name}__genes__${alt_DB}_srst2__results.txt"
				if [[ -f "${processed}/${project}/${sample_name}/srst2/${sample_name}_${alt_DB}__fullgenes__${alt_DB}_srst2__results.txt" ]]; then
					mv "${processed}/${project}/${sample_name}/srst2/${sample_name}_${alt_DB}__fullgenes__${alt_DB}_srst2__results.txt" "${processed}/${project}/${sample_name}/srst2/${sample_name}__genes__${alt_DB}_srst2__results.txt"
				fi
				continue
		elif [[ -f "${processed}/${project}/${sample_name}/srst2/${sample_name}_${alt_DB}__fullgenes__${alt_DB}_srst2__results.txt" ]]; then
			mv "${processed}/${project}/${sample_name}/srst2/${sample_name}_${alt_DB}__fullgenes__${alt_DB}_srst2__results.txt" "${processed}/${project}/${sample_name}/srst2/${sample_name}__genes__${alt_DB}_srst2__results.txt"
			continue
		fi

		if [[ ! -f "${processed}/${project}/${sample_name}/srst2/${sample_name}__genes__${alt_DB}_srst2__results.txt" ]] || [[ ! -f "${processed}/${project}/${sample_name}/srst2/${sample_name}__fullgenes__${alt_DB}_srst2__results.txt" ]]; then
			echo  "Index is below max submissions, submitting"
			echo -e "#!/bin/bash -l\n" > "${main_dir}/srst2AR_${sample}_${start_time}.sh"
			echo -e "#$ -o srst2AR_${sample}.out" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
			echo -e "#$ -e srst2AR_${sample}.err" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
			echo -e "#$ -N srst2AR_${sample}"   >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
			echo -e "#$ -cwd"  >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
			echo -e "#$ -q short.q\n"  >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
			echo -e "module unload Python/2.7" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
			echo -e "module unload Python/3.5.2" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
			echo -e "module unload perl/5.22.1" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
			echo -e "module load Python/2.7.15" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
			echo -e "module load bowtie2/2.2.4" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
			echo -e "module load samtools/0.1.18" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
			echo -e "module load perl/5.16.1-MT" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
			echo -e "module load srst2" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
			# Can we somehow consolidate into one srst2 analysis to do MLST/AR/SEROTYPE
			echo -e "\"${shareScript}/run_srst2_on_singleDB_alternateDB.sh\" \"${sample}\" \"${project}\" \"${alt_DB_path}\"" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
			echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_srst2_complete.txt\"" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
			cd "${main_dir}"
			if [[ "${counter}" -lt "${last_index}" ]]; then
				qsub "${main_dir}/srst2AR_${sample}_${start_time}.sh"
			else
				qsub -sync y "${main_dir}/srst2AR_${sample}_${start_time}.sh"
			fi
		else
			echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_srst2_complete.txt\"" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
			echo "${project}/${sample} already has ${alt_DB}"
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
			if [ -f "${main_dir}/complete/${waiting_sample}_srst2_complete.txt" ]; then
				if [[ ! -f "${processed}/${project}/${sample_name}/srst2/${sample_name}__genes__${alt_DB}_srst2__results.txt" ]] || [[ ! -f "${processed}/${project}/${sample_name}/srst2/${sample_name}__fullgenes__${alt_DB}_srst2__results.txt" ]]; then
					echo "${waiting_sample} has completed, starting ${sample}"
					echo -e "#!/bin/bash -l\n" > "${main_dir}/srst2AR_${sample}_${start_time}.sh"
					echo -e "#$ -o srst2AR_${sample}.out" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
					echo -e "#$ -e srst2AR_${sample}.err" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
					echo -e "#$ -N srst2AR_${sample}"   >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
					echo -e "#$ -cwd"  >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
					echo -e "#$ -q short.q\n"  >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
					echo -e "module unload Python/2.7" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
					echo -e "module unload Python/3.5.2" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
					echo -e "module unload perl/5.22.1" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
					echo -e "module load Python/2.7.15" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
					echo -e "module load bowtie2/2.2.4" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
					echo -e "module load samtools/0.1.18" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
					echo -e "module load perl/5.16.1-MT" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
					echo -e "module load srst2" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
					# Can we somehow consolidate into one srst2 analysis to do MLST/AR/SEROTYPE
					echo -e "\"${shareScript}/run_srst2_on_singleDB_alternateDB.sh\" \"${sample}\" \"${project}\" \"${alt_DB_path}\"" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
					echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_srst2_complete.txt\"" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
					cd "${main_dir}"
					if [[ "${counter}" -lt "${last_index}" ]]; then
						qsub "${main_dir}/srst2AR_${sample}_${start_time}.sh"
					else
						qsub -sync y "${main_dir}/srst2AR_${sample}_${start_time}.sh"
					fi
				else
					echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_srst2_complete.txt\"" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
					echo "${project}/${sample} already has 20180608"
				fi
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
exit 0
