#!/bin/sh -l

#$ -o ablmq_tax.out
#$ -e ablmq_tax..err
#$ -N ablb1
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
. ./config.sh
#Import the module file that loads all necessary mods
. "${mod_changers}/pipeline_mods"

#
# Usage ./abl_mass_qsub_taxonomy.sh path_to_list max_concurrent_submission output_directory_for_scripts
#

# Number regex to test max concurrent submission parametr
number='^[0-9]+$'

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./abl_mass_qsub_taxonomy.sh path_to_list_file(single sample ID per line, e.g. B8VHY/1700128 (it must include project id also)) max_submissions output_directory_for_scripts"
	echo "Output location varies depending on which tasks are performed but will be found somewhere under ${processed}"
	exit 0
elif ! [[ ${2} =~ $number ]] || [[ -z "${2}" ]]; then
	echo "${2} is not a number or is empty. Please input max number of concurrent qsub submissions...exiting"
	exit 2
fi

arr=()
while IFS= read -r line || [ "$line" ];  do
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
main_dir="${3}/taxonomy_subs"
if [[ ! -d "${3}/taxonomy_subs" ]]; then
	mkdir "${3}/taxonomy_subs"
	mkdir "${3}/taxonomy_subs/complete"
elif [[ ! -d  "${3}/taxonomy_subs/complete" ]]; then
	mkdir "${3}/taxonomy_subs/complete"
fi

while [ ${counter} -lt ${arr_size} ] ; do
	sample=$(echo "${arr[${counter}]}" | cut -d'/' -f2)
	project=$(echo "${arr[${counter}]}" | cut -d'/' -f1)
	echo "${counter} of ${arr_size}"
	if [[ -f "${processed}/${project}/${sample}/${sample}_no_calc.tax" ]]; then
		rm "${processed}/${project}/${sample}/${sample}_no_calc.tax"
	fi
	if [[ -d "${processed}/${project}/${sample}/ANI_NO_CALCO)" ]]; then
		rm -r "${processed}/${project}/${sample}/ANI_NO_CALCO)"
	fi
	if [[ -f "${processed}/${project}/${sample}/${sample}_preANI.tax" ]]; then
		rm "${processed}/${project}/${sample}/${sample}_preANI.tax"
	fi
	if [ ${counter} -lt ${max_subs} ]; then
		echo  "Index is below max submissions, submitting"
		echo -e "#!/bin/bash -l\n" > "${main_dir}/tax_${sample}_${start_time}.sh"
		echo -e "#$ -o tax_${sample}.out" >> "${main_dir}/tax_${sample}_${start_time}.sh"
		echo -e "#$ -e tax_${sample}.err" >> "${main_dir}/tax_${sample}_${start_time}.sh"
		echo -e "#$ -N tax_${sample}"   >> "${main_dir}/tax_${sample}_${start_time}.sh"
		echo -e "#$ -cwd"  >> "${main_dir}/tax_${sample}_${start_time}.sh"
		echo -e "#$ -q short.q\n"  >> "${main_dir}/tax_${sample}_${start_time}.sh"
		echo -e "\"${shareScript}/determine_taxID.sh\" \"${sample}\" \"${project}\"" >> "${main_dir}/tax_${sample}_${start_time}.sh"
		echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_taxonomy_complete.txt\"" >> "${main_dir}/tax_${sample}_${start_time}.sh"
		cd "${main_dir}"
		if [[ "${counter}" -lt "${last_index}" ]]; then
			qsub "${main_dir}/tax_${sample}_${start_time}.sh"
		else
			qsub -sync y "${main_dir}/tax_${sample}_${start_time}.sh"
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
			if [ -f "${main_dir}/complete/${waiting_sample}_taxonomy_complete.txt" ]; then
				echo  "Index is below max submissions, submitting"
				echo -e "#!/bin/bash -l\n" > "${main_dir}/tax_${sample}_${start_time}.sh"
				echo -e "#$ -o tax_${sample}.out" >> "${main_dir}/tax_${sample}_${start_time}.sh"
				echo -e "#$ -e tax_${sample}.err" >> "${main_dir}/tax_${sample}_${start_time}.sh"
				echo -e "#$ -N tax_${sample}"   >> "${main_dir}/tax_${sample}_${start_time}.sh"
				echo -e "#$ -cwd"  >> "${main_dir}/tax_${sample}_${start_time}.sh"
				echo -e "#$ -q short.q\n"  >> "${main_dir}/tax_${sample}_${start_time}.sh"
				echo -e "\"${shareScript}/determine_taxID.sh\" \"${sample}\" \"${project}\"" >> "${main_dir}/tax_${sample}_${start_time}.sh"
				echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_taxonomy_complete.txt\"" >> "${main_dir}/tax_${sample}_${start_time}.sh"
				cd "${main_dir}"
				if [[ "${counter}" -lt "${last_index}" ]]; then
					qsub "${main_dir}/tax_${sample}_${start_time}.sh"
				else
					qsub -sync y "${main_dir}/tax_${sample}_${start_time}.sh"
				fi
				break
			fi
		done
	fi
	counter=$(( counter + 1 ))
done

# Check for completion of all samples
timer=0
for item in "${arr[@]}"; do
	waiting_sample=$(echo "${item}" | cut -d'/' -f2)
	if [[ -f "${main_dir}/complete/${waiting_sample}_taxonomy_complete.txt" ]]; then
		echo "${item} is complete"
		if [[ -f "${shareScript}/tax_${waiting_sample}.out" ]]; then
			mv "${shareScript}/tax_${waiting_sample}.out" "${main_dir}"
		fi
		if [[ -f "${shareScript}/tax_${waiting_sample}.err" ]]; then
			mv "${shareScript}/tax_${waiting_sample}.err" "${main_dir}"
		fi
	else
		while :
		do
				if [[ ${timer} -gt 3600 ]]; then
					echo "Timer exceeded limit of 3600 seconds = 60 minutes"
					exit 1
				fi
				if [[ -f "${main_dir}/complete/${waiting_sample}_taxonomy_complete.txt" ]]; then
					echo "${item} is complete"
					if [[ -f "${shareScript}/tax_${waiting_sample}.out" ]]; then
						mv "${shareScript}/tax_${waiting_sample}.out" "${main_dir}"
					fi
					if [[ -f "${shareScript}/tax_${waiting_sample}.err" ]]; then
						mv "${shareScript}/tax_${waiting_sample}.err" "${main_dir}"
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
printf "%s %s" "abl_mass_qsub_taxonomy.sh has completed ${2}" "${global_end_time}" | mail -s "abl_mass_qsub complete" nvx4@cdc.gov
exit 0
