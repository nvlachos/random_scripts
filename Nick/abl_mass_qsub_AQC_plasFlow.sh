#!/bin/sh -l

#$ -o ablmq-AQC.out
#$ -e ablmq-AQC.err
#$ -N ablmq-AQC
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
# Usage ./abl_mass_qsub_AQC.sh path_to_list max_concurrent_submissions output_directory_for_scripts clobberness[keep|clobber]
#

# Number regex to test max concurrent submission parametr
number='^[0-9]+$'

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./abl_mass_qsub_AQC.sh path_to_list_file(single sample ID per line, e.g. B8VHY/1700128 (it must include project id also)) max_concurrent_submissions output_directory_for_scripts clobberness[keep|clobber]"
	exit 1
elif [[ ! -f "${1}" ]]; then
	echo "${1} (list) does not exist...exiting"
	exit 1
elif ! [[ ${2} =~ $number ]] || [[ -z "${2}" ]]; then
	echo "${2} is not a number or is empty. Please input max number of concurrent qsub submissions...exiting"
	exit 2
elif [[ -z "${3}" ]]; then
	echo "Output directory parameter empty...exiting"
	exit 1
elif [[ -z "${4}" ]]; then
	echo "Clobberness is empty...exiting"
	exit 1
fi

# Check that clobberness is a valid option
if [[ "${4}" != "keep" ]] && [[ "${4}" != "clobber" ]]; then
	echo "Clobberness was not input correctly [keep|clobber]...exiting"
	exit 1
else
	clobberness="${4}"
fi

# Create an array of all samples in the list
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
main_dir="${3}/ACS_plasFlow_subs"
if [[ ! -d "${3}/ACS_plasFlow_subs" ]]; then
	mkdir "${3}/ACS_plasFlow_subs"
	mkdir "${3}/ACS_plasFlow_subs/complete"
elif [[ ! -d "${3}/ACS_plasFlow_subs/complete" ]]; then
	mkdir "${3}/ACS_plasFlow_subs/complete"
fi

start_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")

# Create and submit qsub scripts to do Assembly quality checks of all isolates on the list
while [ ${counter} -lt ${arr_size} ] ; do
	sample=$(echo "${arr[${counter}]}" | cut -d'/' -f2)
	project=$(echo "${arr[${counter}]}" | cut -d'/' -f1)
	# Delete old files if clobber is set
	if [[ "${clobberness}" == "clobber" ]]; then
		if [[ "${processed}/${project}/${sample}/Assembly_Stats_plasFlow" ]]; then
			rm -r "${processed}/${project}/${sample}/Assembly_Stats_plasFlow"
		fi
	fi
	# Check if there is a usable assembly to work with
	if [[ -s "${processed}/${project}/${sample}//plasFlow/Unicycler_assemblies/${sample}_uni_assembly/${sample}_plasmid_assembly_trimmed.fasta" ]]; then
		# Check if counter is below max submission limit
		if [[ ${counter} -lt ${max_subs} ]]; then
			# Check if old results exist, skip if so
			if [[ ! -f "${processed}/${project}/${sample}/Assembly_Stats_plasFlow/${sample}_report.txt" ]] ; then
				echo  "Index is below max submissions, submitting"
				echo -e "#!/bin/bash -l\n" > "${main_dir}/ACSP_${sample}_${start_time}.sh"
				echo -e "#$ -o ACSP_${sample}.out" >> "${main_dir}/ACSP_${sample}_${start_time}.sh"
				echo -e "#$ -e ACSP_${sample}.err" >> "${main_dir}/ACSP_${sample}_${start_time}.sh"
				echo -e "#$ -N ACSP_${sample}"   >> "${main_dir}/ACSP_${sample}_${start_time}.sh"
				echo -e "#$ -cwd"  >> "${main_dir}/ACSP_${sample}_${start_time}.sh"
				echo -e "#$ -q short.q\n"  >> "${main_dir}/ACSP_${sample}_${start_time}.sh"
				echo -e "cd ${shareScript}" >> "${main_dir}/ACSP_${sample}_${start_time}.sh"
				echo -e "\"${shareScript}/run_Assembly_Quality_Check.sh\" \"${sample}\" \"${project}\" \"-p\"" >> "${main_dir}/ACSP_${sample}_${start_time}.sh"
				echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_ACSP_complete.txt\"" >> "${main_dir}/ACSP_${sample}_${start_time}.sh"
				cd "${main_dir}"
				qsub "${main_dir}/ACSP_${sample}_${start_time}.sh"
			# Old data exists, skipping
			else
				echo "${project}/${sample} already has QUAST summaries"
				echo "$(date)" > "${main_dir}/complete/${sample}_ACSP_complete.txt"
			fi
		# Counter above max submissions, must wait for slot to optn up
		else
			waiting_for_index=$(( counter - max_subs ))
			waiting_sample=$(echo "${arr[${waiting_for_index}]}" | cut -d'/' -f2)
			timer=0
			echo "Index is above max submissions, waiting for index ${waiting_for_index}:${waiting_sample} to complete"
			# Loop to wait until slot opens up
			while :
			do
				# Check if max time limit has been reached
				if [[ ${timer} -gt 1800 ]]; then
					echo "Timer exceeded limit of 1800 seconds 30 minutes"
					break
				fi
				# Check if "waiting" sample has completed
				if [[ -f "${main_dir}/complete/${waiting_sample}_ACSP_complete.txt" ]] || [[ ! -f "${processed}/${project}/${waiting_sample}/plasFlow/Unicycler_assemblies/${waiting_sample}_uni_assembly/${waiting_sample}_plasmid_assembly_trimmed.fasta" ]]; then
					# Check if old data exists, skip if so
					if [[ ! -f "${processed}/${project}/${sample}/Assembly_Stats_plasFlow/${sample}_report.txt" ]] ; then
						echo  "Index is below max submissions, submitting"
						echo -e "#!/bin/bash -l\n" > "${main_dir}/ACSP_${sample}_${start_time}.sh"
						echo -e "#$ -o ACSP_${sample}.out" >> "${main_dir}/ACSP_${sample}_${start_time}.sh"
						echo -e "#$ -e ACSP_${sample}.err" >> "${main_dir}/ACSP_${sample}_${start_time}.sh"
						echo -e "#$ -N ACSP_${sample}"   >> "${main_dir}/ACSP_${sample}_${start_time}.sh"
						echo -e "#$ -cwd"  >> "${main_dir}/ACSP_${sample}_${start_time}.sh"
						echo -e "#$ -q short.q\n"  >> "${main_dir}/ACSP_${sample}_${start_time}.sh"
						echo -e "cd ${shareScript}" >> "${main_dir}/ACSP_${sample}_${start_time}.sh"
						echo -e "\"${shareScript}/run_Assembly_Quality_Check.sh\" \"${sample}\" \"${project}\" \"-p\"" >> "${main_dir}/ACSP_${sample}_${start_time}.sh"
						echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_ACSP_complete.txt\"" >> "${main_dir}/ACSP_${sample}_${start_time}.sh"
						cd "${main_dir}"
						qsub "${main_dir}/ACSP_${sample}_${start_time}.sh"
					# Old data exists, skipping
					else
						echo "${project}/${sample} already has QUAST summaries"
						echo "$(date)" > "${main_dir}/complete/${sample}_ACSP_complete.txt"
					fi
					break
				# Wait 5 seconds before checking if :waiting: sample has completed yet
				else
					timer=$(( timer + 5 ))
					echo "sleeping for 5 seconds, so far slept for ${timer}"
					sleep 5
				fi
			done
		fi
	# No Assembly file to work with
	else
		echo "${project}/${sample} does not have assembly"
		echo "$(date)" > "${main_dir}/complete/${sample}_ACSP_complete.txt"
	fi
	counter=$(( counter + 1 ))
done

# Loop to ensure all samples are complete (or time runs) before allowing the script to exit
timer=0
for item in "${arr[@]}"; do
	waiting_sample=$(echo "${item}" | cut -d'/' -f2)
	if [[ -f "${main_dir}/complete/${waiting_sample}_ACSP_complete.txt" ]] || [[ ! -f "${processed}/${project}/${waiting_sample}/plasFlow/Unicycler_assemblies/${waiting_sample}_uni_assembly/${waiting_sample}_plasmid_assembly_trimmed.fasta" ]]; then
		echo "${item} is complete"
		if [[ -f "${shareScript}/ACSP_${waiting_sample}.out" ]]; then
			mv "${shareScript}/ACSP_${waiting_sample}.out" ${main_dir}
		fi
		if [[ -f "${shareScript}/ACSP_${waiting_sample}.err" ]]; then
			mv "${shareScript}/ACSP_${waiting_sample}.err" ${main_dir}
		fi
	else
		# Check every 5 seconds to see if the sample has completed normal csstar analysis
		while :
		do
				if [[ ${timer} -gt 1800 ]]; then
					echo "Timer exceeded limit of 1800 seconds = 30 minutes"
					exit 1
				fi
				if [[ -f "${main_dir}/complete/${waiting_sample}_ACSP_complete.txt" ]]; then
					echo "${item} is complete"
					if [[ -f "${shareScript}/ACSP_${waiting_sample}.out" ]]; then
						mv "${shareScript}/ACSP_${waiting_sample}.out" ${main_dir}
					fi
					if [[ -f "${shareScript}/ACSP_${waiting_sample}.err" ]]; then
						mv "${shareScript}/ACSP_${waiting_sample}.err" ${main_dir}
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
printf "%s %s" "abl_mass_qsub_AQC.sh has completed" "${global_end_time}" | mail -s "abl_mass_qsub_AQC.sh complete" nvx4@cdc.gov

#Script exited gracefully (unless something else inside failed)
exit 0
