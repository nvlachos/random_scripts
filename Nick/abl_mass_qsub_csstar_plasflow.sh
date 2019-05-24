#!/bin/sh -l

#$ -o ablmq-cs.out
#$ -e ablmq-cs.err
#$ -N ablmq-cs
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
# Usage ./abl_mass_qsub_csstar.sh path_to_list max_concurrent_submissions output_folder_for_scripts clobberness[keep|clobber]
#

# Number regex to test max concurrent submission parametr
number='^[0-9]+$'

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./abl_mass_qsub_csstar.sh path_to_list_file(single sample ID per line, e.g. B8VHY/1700128 (it must include project id also)) max_concurrent_submissions output_directory_for_scripts clobberness[keep|clobber]"
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
while IFS= read -r line || [[ "$line" ]];  do
	if [[ ! -z "${line}" ]]; then
		line=$(echo ${line} | tr -d '\n' | tr -d '\r')
		arr+=($line)
	fi
done < ${1}

arr_size="${#arr[@]}"
last_index=$(( arr_size -1 ))
echo "-${arr_size}:${arr[@]}-"


# Create direcory to hold all temporary qsub scripts
counter=0
max_subs=${2}

# format name being extracted from alt database
main_dir="${3}/csstar_plasflow_subs"
#cp ./config ${main_dir}
if [[ ! -d "${3}/csstar_plasflow_subs" ]]; then
	mkdir -p "${3}/csstar_plasflow_subs/complete"
elif [[ ! -d  "${3}/csstar_plasflow_subs/complete" ]]; then
	mkdir "${3}/csstar_plasflow_subs/complete"
fi

start_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")

# Create and submit scripts to run default csstar on all samples on the list
while [ ${counter} -lt ${arr_size} ] ; do
	sample=$(echo "${arr[${counter}]}" | cut -d'/' -f2)
	project=$(echo "${arr[${counter}]}" | cut -d'/' -f1)
	if [[ "${clobberness}" == "clobber" ]]; then
		rm -r ${processed}/${project}/${sample}/c-sstar_plasFlow
	fi
	#echo ${counter}-${project}-${sample}
	# Check if sample has a usable assembly file
	if [[ -s "${processed}/${project}/${sample}/plasFlow/Unicycler_assemblies/${sample}_uni_assembly/${sample}_plasmid_assembly_trimmed.fasta" ]]; then
		#echo "Test"
		# Check if counter is below max number of concurrent submissions
		if [[ ${counter} -lt ${max_subs} ]]; then
			# Check if old data exists, skip if so
			if [[ ! -f "${processed}/${project}/${sample}/c-sstar_plasFlow/${sample}.${resGANNOT_srst2_filename}.gapped_40_sstar_summary.txt" ]]; then
				echo  "Index is below max submissions, submitting"
				echo -e "#!/bin/bash -l\n" > "${main_dir}/csstp_${sample}_${start_time}.sh"
				echo -e "#$ -o csstp_${sample}.out" >> "${main_dir}/csstp_${sample}_${start_time}.sh"
				echo -e "#$ -e csstp_${sample}.err" >> "${main_dir}/csstp_${sample}_${start_time}.sh"
				echo -e "#$ -N csstp_${sample}"   >> "${main_dir}/csstp_${sample}_${start_time}.sh"
				echo -e "#$ -cwd"  >> "${main_dir}/csstp_${sample}_${start_time}.sh"
				echo -e "#$ -q short.q\n"  >> "${main_dir}/csstp_${sample}_${start_time}.sh"
				echo -e "module load Python/3.6.1\n" >> "${main_dir}/csstp_${sample}_${start_time}.sh"
				# Defaulting to gapped/98, change if you want to include user preferences
				echo -e "cd ${shareScript}" >> "${main_dir}/csstp_${sample}_${start_time}.sh"
				echo -e "\"${shareScript}/run_c-sstar_on_single_plasFlow.sh\" \"${sample}\" g o \"${project}\" \"-p\"" >> "${main_dir}/csstp_${sample}_${start_time}.sh"
				echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_csstarp_complete.txt\"" >> "${main_dir}/csstp_${sample}_${start_time}.sh"
				cd "${main_dir}"
				echo "submitting ${main_dir}/csstp_${sample}_${start_time}.sh"
				qsub "${main_dir}/csstp_${sample}_${start_time}.sh"
			# Old data exists
			else
				echo "${project}/${sample} already has the newest ResGANNOT for csstar plasFlow(${resGANNOT_srst2_filename})"
				echo -e "$(date)" > "${main_dir}/complete/${sample}_csstarp_complete.txt"
			fi

		# Counter is above max number of submissions
		else
			waiting_for_index=$(( counter - max_subs ))
			waiting_sample=$(echo "${arr[${waiting_for_index}]}" | cut -d'/' -f2)
			timer=0
			echo "Index is above max submissions, waiting for index ${waiting_for_index}:${waiting_sample} to complete"
			while :
			do
				# Check that the timer has not exceeded max amount of time to wait
				if [[ ${timer} -gt 1800 ]]; then
					echo "Timer exceeded limit of 1800 seconds 30 minutes"
					break
				fi
				# Check if usable assembly exists for current sample or that one does not exist for the waiting sample (therefore there would be no need to wait on it)
				if [[ -f "${main_dir}/complete/${waiting_sample}_csstarp_complete.txt" ]] || [[ ! -s "${processed}/${project}/${waiting_sample}/plasFlow/Unicycler_assemblies/${sample}_uni_assembly/${sample}_plasmid_assembly_trimmed.fasta" ]]; then
					# Check if old data exists, skip if so
					if [[ ! -f "${processed}/${project}/${sample}/c-sstar_plasFlow/${sample}.${resGANNOT_srst2_filename}.gapped_40_sstar_summary.txt" ]]; then
						echo  "Index is below max submissions, submitting"
						echo -e "#!/bin/bash -l\n" > "${main_dir}/csstp_${sample}_${start_time}.sh"
						echo -e "#$ -o csstp_${sample}.out" >> "${main_dir}/csstp_${sample}_${start_time}.sh"
						echo -e "#$ -e csstp_${sample}.err" >> "${main_dir}/csstp_${sample}_${start_time}.sh"
						echo -e "#$ -N csstp_${sample}"   >> "${main_dir}/csstp_${sample}_${start_time}.sh"
						echo -e "#$ -cwd"  >> "${main_dir}/csstp_${sample}_${start_time}.sh"
						echo -e "#$ -q short.q\n"  >> "${main_dir}/csstp_${sample}_${start_time}.sh"
						echo -e "module load Python/3.6.1\n" >> "${main_dir}/csstp_${sample}_${start_time}.sh"
						# Defaulting to gapped/98, change if you want to include user preferences
						echo -e "cd ${shareScript}" >> "${main_dir}/csstp_${sample}_${start_time}.sh"
						echo -e "\"${shareScript}/run_c-sstar_on_single_plasFlow.sh\" \"${sample}\" g o \"${project}\" \"-p\"" >> "${main_dir}/csstp_${sample}_${start_time}.sh"
						echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_csstarp_complete.txt\"" >> "${main_dir}/csstp_${sample}_${start_time}.sh"
						cd ${main_dir}
						qsub "${main_dir}/csstp_${sample}_${start_time}.sh"
					# Skipping because old data exists
					else
						echo "${project}/${sample} already has the newest ResGANNOT on csstar plasFlow (${resGANNOT_srst2_filename})"
						echo -e "$(date)" > "${main_dir}/complete/${sample}_csstarp_complete.txt"
					fi
					break
				# If waiting sample has not completed, wait 5 more seconds and try again
				else
					timer=$(( timer + 5 ))
					echo "${main_dir}/complete/${waiting_sample}_csstarp_complete.txt not ready, sleeping for 5 seconds, so far slept for ${timer}"
					sleep 5
				fi
			done
		fi
	else
		echo "${project}/${sample} has no plasflow assembly to run csstar on (! ${processed}/${project}/${sample}/plasFlow/Unicycler_assemblies/${sample}_uni_assembly/${sample}_plasmid_assembly_trimmed.fasta)"
		echo "$(date)" > "${main_dir}/complete/${sample}_csstarp_complete.txt"
	fi
	counter=$(( counter + 1 ))
done

# Loop to ensure all samples are complete (or time runs) before allowing the script to exit
timer=0
ptimer=0
for item in "${arr[@]}"; do
	waiting_sample=$(echo "${item}" | cut -d'/' -f2)
	if [[ -f "${main_dir}/complete/${waiting_sample}_csstarp_complete.txt" ]] || [[ ! -s "${processed}/${project}/${waiting_sample}/plasFlow/Unicycler_assemblies/${waiting_sample}_uni_assembly/${waiting_sample}_plasmid_assembly_trimmed.fasta" ]]; then
		echo "${item} is complete normal"
		if [[ -f "${shareScript}/csstp_${sample}.out" ]]; then
			mv "${shareScript}/csstp_${sample}.out" "${main_dir}"
		fi
		if [[ -f "${shareScript}/csstp_${sample}.err" ]]; then
			mv "${shareScript}/csstp_${sample}.err" "${main_dir}"
		fi
	else
		# Check every 5 seconds to see if the sample has completed normal csstar analysis
		while :
		do
				if [[ ${timer} -gt 1800 ]]; then
					echo "Timer exceeded limit of 1800 seconds = 30 minutes"
					exit 1
				fi
				if [[ -f "${main_dir}/complete/${waiting_sample}_csstarp_complete.txt" ]]; then
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
global_end_time=$(date "+%m-%d-%Y @ %Hh_%Mm_%Ss")
printf "%s %s" "abl_mass_qsub_csstar.sh has completed" "${global_end_time}" | mail -s "abl_mass_qsub_csstar.sh complete" nvx4@cdc.gov

#Script exited gracefully (unless something else inside failed)
exit 0
