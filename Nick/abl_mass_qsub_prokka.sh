#!/bin/sh -l

#$ -o ablmq-prok.out
#$ -e ablmq-prok.err
#$ -N ablmq-prok
#$ -cwd
#$ -q short.q

##Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi

. ./config.sh
#Import the module file that loads all necessary mods
. "${mod_changers}/pipeline_mods"

#List all currently loaded modules
#. ./module_changers/list_modules.sh

#
# Usage ./abl_mass_qsub_prokka.sh path_to_list max_concurrent_submissions output_directory_for_scripts clobberness[keep|clobber]
#

# Number regex to test max concurrent submission parametr
number='^[0-9]+$'

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./abl_mass_qsub_prokka.sh path_to_list_file(single sample ID per line, e.g. B8VHY/1700128 (it must include project id also)) max_concurrent_submissions output_directory_for_scripts"
	exit 1
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

# create an array of all samples in the list
arr=()
while IFS= read -r line || [ "$line" ];  do
  arr+=("$line")
done < ${1}

arr_size="${#arr[@]}"
last_index=$(( arr_size -1 ))
echo "-${arr_size}:${arr[@]}-"


# Create direcory to hold all temporary qsub scripts
counter=0
max_subs=${2}

# Set script directory
main_dir="${3}/prokka_subs"
if [[ ! -d "${3}/prokka_subs" ]]; then
	mkdir "${3}/prokka_subs"
	mkdir "${3}/prokka_subs/complete"
elif [[ ! -d "${3}/prokka_subs/complete" ]]; then
	mkdir "${3}/prokka_subs/complete"
fi

start_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")

# Creates and submits qsub scripts to try to find plasmid replicons in all the samples in the list
while [ ${counter} -lt ${arr_size} ] ; do
	sample=$(echo "${arr[${counter}]}" | cut -d'/' -f2)
	project=$(echo "${arr[${counter}]}" | cut -d'/' -f1)
	if [[ "${clobberness}" == "clobber" ]]; then
		if [[ -d "${processed}/${project}/${sample}/prokka" ]]; then
			rm -r "${processed}/${project}/${sample}/prokka"
		fi
		# if [[ -d "${processed}/${project}/${sample}/plasmid_on_plasFlow" ]]; then
		# 	rm -r "${processed}/${project}/${sample}/plasmid_on_plasFlow"
		# fi
	fi
	#echo ${counter}
	# Check for assembly file to run plasmidFinder on
	if [[ -s "${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed.fasta" ]]; then
		# Check if counter is below max sub limit
		if [[ ${counter} -lt ${max_subs} ]]; then
			# Check for old data, skip if present
			if [[ ! -f "${processed}/${project}/${sample}/prokka/${sample}_PROKKA.gbf" ]]; then
				echo  "Index is below max submissions, submitting"
				echo -e "#!/bin/bash -l\n" > "${main_dir}/prok_${sample}_${start_time}.sh"
				echo -e "#$ -o prok_${sample}.out" >> "${main_dir}/prok_${sample}_${start_time}.sh"
				echo -e "#$ -e prok_${sample}.err" >> "${main_dir}/prok_${sample}_${start_time}.sh"
				echo -e "#$ -N prok_${sample}"   >> "${main_dir}/prok_${sample}_${start_time}.sh"
				echo -e "#$ -cwd"  >> "${main_dir}/prok_${sample}_${start_time}.sh"
				echo -e "#$ -q short.q\n"  >> "${main_dir}/prok_${sample}_${start_time}.sh"
				echo -e "cd ${shareScript}" >> "${main_dir}/prok_${sample}_${start_time}.sh"
				echo -e "mv ${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed.fasta ${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed.fasta_backup" >> "${main_dir}/prok_${sample}_${start_time}.sh"
				echo -e "python3 \"${shareScript}/fasta_headers.sh\" -i \"${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed.fasta_backup\" -o \"${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed.fasta\" --reverse" >> "${main_dir}/prok_${sample}_${start_time}.sh"
				echo -e "\"${shareScript}/run_prokka.sh\" \"${sample}\" \"${project}\"" >> "${main_dir}/prok_${sample}_${start_time}.sh"
				echo -e "rm ${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed.fasta" >> "${main_dir}/prok_${sample}_${start_time}.sh"
				echo -e "mv ${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed.fasta_backup ${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed.fasta" >> "${main_dir}/prok_${sample}_${start_time}.sh"
				echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_prok_complete.txt\"" >> "${main_dir}/prok_${sample}_${start_time}.sh"
				cd "${main_dir}"
				qsub "${main_dir}/prok_${sample}_${start_time}.sh"
			# Old data exists, skipping plasmidFinder on normal Assembly
			else
				echo "${project}/${sample} already has prokka output"
				echo "$(date)" > "${main_dir}/complete/${sample}_prok_complete.txt"
			fi
		# Counter is above max sub limit, must wait until slot becomes available
		else
			waiting_for_index=$(( counter - max_subs ))
			waiting_sample=$(echo "${arr[${waiting_for_index}]}" | cut -d'/' -f2)
			timer=0
			echo "Index is above max submissions, waiting for index ${waiting_for_index}:${waiting_sample} to complete"
			# Loop to check if "waiting" sample is complete
			while :
			do
				# Check if timer has exceeded maximum elapsed allowed time
				if [[ ${timer} -gt 1800 ]]; then
					echo "Timer exceeded limit of 1800 seconds 30 minutes"
					break
				fi
				# Check if "waiting" sample has completed
				if [[ -f "${main_dir}/complete/${waiting_sample}_prok_complete.txt" ]] || [[ ! -s "${processed}/${project}/${waiting_sample}/Assembly/${waiting_sample}_scaffolds_trimmed.fasta" ]]; then
					# Check if old data exists
					if [[ ! -f "${processed}/${project}/${sample}/prokka/${sample}_PROKKA.gbf" ]]; then
						echo  "Index is below max submissions, submitting"
						echo -e "#!/bin/bash -l\n" > "${main_dir}/prok_${sample}_${start_time}.sh"
						echo -e "#$ -o prok_${sample}.out" >> "${main_dir}/prok_${sample}_${start_time}.sh"
						echo -e "#$ -e prok_${sample}.err" >> "${main_dir}/prok_${sample}_${start_time}.sh"
						echo -e "#$ -N prok_${sample}"   >> "${main_dir}/prok_${sample}_${start_time}.sh"
						echo -e "#$ -cwd"  >> "${main_dir}/prok_${sample}_${start_time}.sh"
						echo -e "#$ -q short.q\n"  >> "${main_dir}/prok_${sample}_${start_time}.sh"
						echo -e "cd ${shareScript}" >> "${main_dir}/prok_${sample}_${start_time}.sh"
						echo -e "mv ${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed.fasta ${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed.fasta_backup" >> "${main_dir}/prok_${sample}_${start_time}.sh"
						echo -e "python3 \"${shareScript}/fasta_headers.sh\" -i \"${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed.fasta_backup\" -o \"${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed.fasta\" --reverse" >> "${main_dir}/prok_${sample}_${start_time}.sh"
						echo -e "\"${shareScript}/run_prokka.sh\" \"${sample}\" \"${project}\"" >> "${main_dir}/prok_${sample}_${start_time}.sh"
						echo -e "rm ${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed.fasta" >> "${main_dir}/prok_${sample}_${start_time}.sh"
						echo -e "mv ${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed.fasta_backup ${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed.fasta" >> "${main_dir}/prok_${sample}_${start_time}.sh"
						echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_prok_complete.txt\"" >> "${main_dir}/prok_${sample}_${start_time}.sh"
						cd "${main_dir}"
						qsub "${main_dir}/prok_${sample}_${start_time}.sh"
					# Old data exists, skipping plasmid Finder on normal Assembly
					else
						echo "${project}/${sample} already has prokka output"
						echo "$(date)" > "${main_dir}/complete/${sample}_prok_complete.txt"
					fi
					break
				# Wait 5 seconds before seeing if the "waiting" sample is done yet
				else
					timer=$(( timer + 5 ))
					echo "sleeping for 5 seconds, so far slept for ${timer}"
					sleep 5
				fi
			done
		fi
	# No assembly exists to run prokka on
	else
		echo "${project}/${sample} has no assembly to run prokka on"
		echo "$(date)" > "${main_dir}/complete/${sample}_prok_complete.txt"
	fi
	counter=$(( counter + 1 ))
done

# Check for completion of all samples
timer=0
for item in "${arr[@]}"; do
	waiting_sample=$(echo "${item}" | cut -d'/' -f2)
	if [[ -f "${main_dir}/complete/${waiting_sample}_prok_complete.txt" ]] || [[ ! -s "${processed}/${project}/${waiting_sample}/Assembly/${waiting_sample}_scaffolds_trimmed.fasta" ]]; then
		echo "${item} is complete"
		if [[ -f "${shareScript}/prok_${waiting_sample}.out" ]]; then
			mv "${shareScript}/prok_${waiting_sample}.out" ${main_dir}
		fi
		if [[ -f "${shareScript}/prok_${waiting_sample}.err" ]]; then
 			mv "${shareScript}/prok_${waiting_sample}.err" ${main_dir}
		fi
	else
		while :
		do
				if [[ ${timer} -gt 3600 ]]; then
					echo "Timer exceeded limit of 3600 seconds = 60 minutes"
					exit 1
				fi
				if [[ -f "${main_dir}/complete/${waiting_sample}_prok_complete.txt" ]]; then
					echo "${item} is complete"
					if [[ -f "${shareScript}/prok_${waiting_sample}.out" ]]; then
						mv "${shareScript}/prok_${waiting_sample}.out" ${main_dir}
					fi
					if [[ -f "${shareScript}/prok_${waiting_sample}.err" ]]; then
			 			mv "${shareScript}/prok_${waiting_sample}.err" ${main_dir}
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
printf "%s %s" "abl_mass_qsub_prokka.sh has completed" "${global_end_time}" | mail -s "abl_mass_qsub_prokka.sh complete" nvx4@cdc.gov

#Script exited gracefully (unless something else inside failed)
exit 0
