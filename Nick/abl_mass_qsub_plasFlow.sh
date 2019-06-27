#!/bin/sh -l

#$ -o ablmq_pFlow.out
#$ -e ablmq_pFlow.err
#$ -N ablmq_pFlow
#$ -cwd
#$ -q short.q

echo "1"
echo $(pwd)
echo "2"
echo $(ls -l)
echo "3"

#Import the config file with shortcuts and settings
if [[ ! -f ./config.sh ]]; then
	cp ./config_template.sh ./config.sh
fi
. ./config.sh
#Import the module file that loads all necessary mods
. "${mod_changers}/pipeline_mods"

#Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi
. ./config.sh
#Import the module file that loads all necessary mods
. "${mod_changers}/pipeline_mods"



#
# Usage ./abl_mass_qsub_plasFlow.sh path_to_list max_concurrent_submission output_directory_for_scripts clobberness[keep|clobber]
#

# Number regex to test max concurrent submission parametr
number='^[0-9]+$'

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./abl_mass_qsub_plasFlow.sh path_to_list_file(single sample ID per line, e.g. B8VHY/1700128 (it must include project id also)) max_submissions output_directory_for_scripts clobberness[keep|clobber]"
	echo "Output location varies depending on which tasks are performed but will be found somewhere under ${processed}"
	exit 0
elif ! [[ ${2} =~ $number ]] || [[ -z "${2}" ]]; then
	echo "${2} is not a number or is empty. Please input max number of concurrent qsub submissions...exiting"
	exit 2
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

arr=()
while IFS= read -r line || [ "$line" ];  do
  arr+=("$line")
done < ${1}

arr_size="${#arr[@]}"
echo "-${arr_size}:${arr[@]}-"

if [[ -z ${2} ]]; then
	max_subs=10
fi

# Loop through and act on each sample name in the passed/provided list
counter=0
max_subs=${2}
main_dir="${3}/plasFlow_subs"
if [[ ! -d "${3}/plasFlow_subs" ]]; then
	mkdir "${3}/plasFlow_subs"
	mkdir "${3}/plasFlow_subs/complete"
elif [[ ! -d  "${3}/plasFlow_subs/complete" ]]; then
	mkdir "${3}/plasFlow_subs/complete"
fi

start_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")

while [ ${counter} -lt ${arr_size} ] ; do
	sample=$(echo "${arr[${counter}]}" | cut -d'/' -f2)
	project=$(echo "${arr[${counter}]}" | cut -d'/' -f1)
	echo ${counter}
	if [[ "${clobberness}" == "clobber" ]]; then
		if [[ -d "${processed}/${project}/${sample}/plasFlow" ]]; then
			rm -r "${processed}/${project}/${sample}/plasFlow"
		fi
	fi
	if [[ ! -f "${processed}/${project}/${sample}/${sample}.tax" ]]; then
		${shareScript}/determine_taxID.sh "${sample}" "${project}"
	fi
	family=""
	while IFS= read -r line  || [ -n "$line" ]; do
		# Grab first letter of line (indicating taxonomic level)
		first=${line:0:1}
		# Assign taxonomic level value from 4th value in line (1st-classification level, 2nd-% by kraken, 3rd-true % of total reads, 4th-identifier)
		if [ "${first}" = "F" ]; then
			family=$(echo "${line}" | awk -F ' ' '{print $2}')
		fi
	done < "${processed}/${project}/${sample}/${sample}.tax"

	if [[ "${family}" == "Enterobacteriaceae" ]]; then
		if [ ${counter} -lt ${max_subs} ]; then
			if [[ ! -f "${processed}/${project}/${sample_name}/plasFlow/Unicycler_assemblies/${sample_name}_uni_assembly/${sample_name}_plasmid_assembly_trimmed.fasta" ]]; then
				echo  "Index is below max submissions, submitting"
				echo -e "#!/bin/bash -l\n" > "${main_dir}/pFlow_${sample}_${start_time}.sh"
				echo -e "#$ -o pFlow_${sample}.out" >> "${main_dir}/pFlow_${sample}_${start_time}.sh"
				echo -e "#$ -e pFlow_${sample}.err" >> "${main_dir}/pFlow_${sample}_${start_time}.sh"
				echo -e "#$ -N pFlow_${sample}"   >> "${main_dir}/pFlow_${sample}_${start_time}.sh"
				echo -e "#$ -cwd"  >> "${main_dir}/pFlow_${sample}_${start_time}.sh"
				echo -e "#$ -q short.q\n"  >> "${main_dir}/pFlow_${sample}_${start_time}.sh"
				# Add all necessary modules
				### echo -e "module load XXX" >> "${main_dir}/pFlow_${sample}_${start_time}.sh"
				echo -e "cd ${shareScript}" >> "${main_dir}/pFlow_${sample}_${start_time}.sh"
				echo -e "\"${shareScript}/run_plasFlow.sh\" \"${sample}\" \"${project}\"" >> "${main_dir}/pFlow_${sample}_${start_time}.sh"
				echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_pFlow_complete.txt\"" >> "${main_dir}/pFlow_${sample}_${start_time}.sh"
				cd "${main_dir}"

				qsub "${main_dir}/pFlow_${sample}_${start_time}.sh"
			else
				echo -e "$(date)" > "${main_dir}/complete/${sample}_pFlow_complete.txt"
				echo "${project}/${sample} already has plasFlow completed"
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
				if [ -f "${main_dir}/complete/${waiting_sample}_pFlow_complete.txt" ]; then
					if [[ ! -f "${processed}/${project}/${sample_name}/plasFlow/Unicycler_assemblies/${sample_name}_uni_assembly/${sample_name}_plasmid_assembly_trimmed.fasta" ]]; then
						echo  "Index is below max submissions, submitting"
						echo -e "#!/bin/bash -l\n" > "${main_dir}/pFlow_${sample}_${start_time}.sh"
						echo -e "#$ -o pFlow_${sample}.out" >> "${main_dir}/pFlow_${sample}_${start_time}.sh"
						echo -e "#$ -e pFlow_${sample}.err" >> "${main_dir}/pFlow_${sample}_${start_time}.sh"
						echo -e "#$ -N pFlow_${sample}"   >> "${main_dir}/pFlow_${sample}_${start_time}.sh"
						echo -e "#$ -cwd"  >> "${main_dir}/pFlow_${sample}_${start_time}.sh"
						echo -e "#$ -q short.q\n"  >> "${main_dir}/pFlow_${sample}_${start_time}.sh"
						# Add all necessary modules
						### echo -e "module load XXX" >> "${main_dir}/pFlow_${sample}_${start_time}.sh"
						echo -e "cd ${shareScript}" >> "${main_dir}/pFlow_${sample}_${start_time}.sh"
						echo -e "\"${shareScript}/run_plasFlow.sh\" \"${sample}\" \"${project}\"" >> "${main_dir}/pFlow_${sample}_${start_time}.sh"
						echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_pFlow_complete.txt\"" >> "${main_dir}/pFlow_${sample}_${start_time}.sh"
						cd "${main_dir}"

						qsub "${main_dir}/pFlow_${sample}_${start_time}.sh"
					else
						echo -e "$(date)" > "${main_dir}/complete/${sample}_pFlow_complete.txt"
						echo "${project}/${sample} already has plasFlow completed"
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
		echo -e "$(date)" > "${main_dir}/complete/${sample}_pFlow_complete.txt"
		echo "${project}/${sample} no plasFlow - not in Enterobacteriaceae family"
	fi
	counter=$(( counter + 1 ))
done

# Check for completion of all samples
timer=0
for item in "${arr[@]}"; do
	waiting_sample=$(echo "${item}" | cut -d'/' -f2)
	if [[ -f "${main_dir}/complete/${waiting_sample}_pFlow_complete.txt" ]]; then
		echo "${item} is complete"
		if [[ -f "${shareScript}/pFlow_${waiting_sample}.out" ]]; then
			mv "${shareScript}/pFlow_${waiting_sample}.out" "${main_dir}"
		fi
		if [[ -f "${shareScript}/pFlow_${waiting_sample}.err" ]]; then
			mv "${shareScript}/pFlow_${waiting_sample}.err" "${main_dir}"
		fi
	else
		while :
		do
				if [[ ${timer} -gt 3600 ]]; then
					echo "Timer exceeded limit of 3600 seconds = 60 minutes"
					exit 1
				fi
				if [[ -f "${main_dir}/complete/${waiting_sample}_pFlow_complete.txt" ]]; then
					echo "${item} is complete"
					if [[ -f "${shareScript}/pFlow_${waiting_sample}.out" ]]; then
						mv "${shareScript}/pFlow_${waiting_sample}.out" "${main_dir}"
					fi
					if [[ -f "${shareScript}/pFlow_${waiting_sample}.err" ]]; then
						mv "${shareScript}/pFlow_${waiting_sample}.err" "${main_dir}"
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
printf "%s %s" "abl_mass_qsub_template.sh has completed ${2}" "${global_end_time}" | mail -s "abl_mass_qsub complete" nvx4@cdc.gov
exit 0
