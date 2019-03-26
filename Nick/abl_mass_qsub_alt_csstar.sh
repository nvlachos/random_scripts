#!/bin/sh -l

#$ -o amq-altc.out
#$ -e amq-altc.err
#$ -N amq-altc
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi
. ./config
#Import the module file that loads all necessary mods
. "${mod_changers}/pipeline_mods"

#
# Usage ./abl_mass_qsub_alt_csstar.sh path_to_list max_concurrent_submissions path_to_alt_database output_directory_for_scripts cloberness (keep|clobber)
#

# Number regex to test max concurrent submission parametr
number='^[0-9]+$'

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to ./abl_mass_qsub_alt_csstar.sh, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./abl_mass_qsub_alt_csstar.sh path_to_list_file(single sample ID per line, e.g. B8VHY/1700128 (it must include project id also)) max_concurrent_submissions path_to_alt_database output_directory_for_scripts"
	exit 1
elif [[ ! -f "${1}" ]]; then
	echo "${1} (list) does not exist...exiting"
	exit 1
elif ! [[ ${2} =~ $number ]] || [[ -z "${2}" ]]; then
	echo "${2} is not a number or is empty. Please input max number of concurrent qsub submissions...exiting"
	exit 2
elif [[ -z "${3}" ]]; then
	echo "${3} (alt_db) does not exist or parameter is empty...exiting"
	exit 1
elif [[ -z "${4}" ]]; then
	echo "No script output directory given...exiting"
	exit 4
elif [[ -z "${5}" ]] ; then
	echo "Clobberness was not input, be sure to add keep or clobber as 4th parameter...exiting"
	exit 1
fi

# Check that clobberness is a valid option
if [[ "${5}" != "keep" ]] && [[ "${5}" != "clobber" ]]; then
	echo "Clobberness was not input, be sure to add keep or clobber as 5th parameter...exiting"
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


# Create direcory to hold all temporary qsub scripts
counter=0
max_subs=${2}

# Set script directory
main_dir="${4}/csstar_alt_subs"
if [[ ! -d "${4}/csstar_alt_subs" ]]; then
	mkdir "${4}/csstar_alt_subs"
	mkdir "${4}/csstar_alt_subs/complete"
elif [[ ! -d  "${4}/csstar_alt_subs/complete" ]]; then
	mkdir "${4}/csstar_alt_subs/complete"
fi

# format name being extracted from alt database
alt_database_path=$(basename -- "${3}")
echo "${alt_database_path}"
alt_database=$(echo ${alt_database_path##*/} | cut -d'.' -f1)
echo "${alt_database}"
alt_database=${alt_database//_srst2/}
echo "${alt_database}"


start_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")

# Create and submit qsub file for each isolate in the list, pausing if sunbissions reach max_concurrent_submissions until space frees up
while [ ${counter} -lt ${arr_size} ] ; do
	echo ${counter}
	sample=$(echo "${arr[${counter}]}" | cut -d'/' -f2)
	project=$(echo "${arr[${counter}]}" | cut -d'/' -f1)
	# Delete old results if clobbering
	if [[ "${5}" == "clobber" ]]; then
		rm ${processed}/${project}/${sample}/c-sstar/${sample}.${alt_database}.gapped_98_sstar_summary
		rm -r ${processed}/${project}/${sample}/c-sstar/${alt_database}_gapped/
	fi

	# Check if there is an assembly to use
	if [[ -s "${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed.fasta" ]]; then
		# Check if counter has reached max submissions yet
		if [[ ${counter} -lt ${max_subs} ]]; then
			echo  "Index is below max submissions, submitting"
			# Check if old results exist, and skip if so
			if [[ ! -f "${processed}/${project}/${sample}/c-sstar/${sample}.${alt_database}.gapped_98_sstar_summary.txt" ]]; then
				echo  "Index is below max submissions, submitting"
				echo -e "#!/bin/bash -l\n" > "${main_dir}/csstn_${sample}_${start_time}.sh"
				echo -e "#$ -o csstn_${sample}.out" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
				echo -e "#$ -e csstn_${sample}.err" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
				echo -e "#$ -N csstn_${sample}"   >> "${main_dir}/csstn_${sample}_${start_time}.sh"
				echo -e "#$ -cwd"  >> "${main_dir}/csstn_${sample}_${start_time}.sh"
				echo -e "#$ -q short.q\n"  >> "${main_dir}/csstn_${sample}_${start_time}.sh"
				echo -e "module load Python/3.6.1\n" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
				# Defaulting to gapped/98, change if you want to include user preferences
				echo -e "\"${shareScript}/run_c-sstar_on_single_alternate_DB.sh\" \"${sample}\" g h \"${project}\" \"${3}\"" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
				echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_csstarn_complete.txt\"" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
				#if [[ "${counter}" -lt "${last_index}" ]]; then
					qsub "${main_dir}/csstn_${sample}_${start_time}.sh"
				#else
				#	if [[ -d "${processed}/${project}/${sample}/c-sstar_plasmid" ]]; then
				#		qsub "${main_dir}/csstn_${sample}_${start_time}.sh"
				#	else
				#		qsub -sync y "${main_dir}/csstn_${sample}_${start_time}.sh"
				#	fi
				#fi

			else
				echo "${project}/${sample} already has 0608"
				echo "$(date)" > "${main_dir}/complete/${sample}_csstarn_complete.txt"
			fi

			# Check if there is a plasmid folder to run also
			if [[ -d "${processed}/${project}/${sample}/c-sstar_plasmid" ]]; then
				# Delete old results if clobbering
				if [[ "${5}" == "clobber" ]]; then
					rm ${processed}/${project}/${sample}/c-sstar_plasmid/${sample}.${alt_database}.gapped_40_sstar_summary
					rm -r ${processed}/${project}/${sample}/c-sstar_plasmid/${alt_database}_gapped/
				fi
				# Check if old results exist, and skip if so
				if [[ ! -f "${processed}/${project}/${sample}/c-sstar_plasmid/${sample}.${alt_database}.gapped_40_sstar_summary.txt" ]]; then
					echo -e "#!/bin/bash -l\n" > "${main_dir}/csstp_${sample}_${start_time}.sh"
					echo -e "#$ -o csstp_${sample}.out" >> "${main_dir}/csstp_${sample}_${start_time}.sh"
					echo -e "#$ -e csstp_${sample}.err" >> "${main_dir}/csstp_${sample}_${start_time}.sh"
					echo -e "#$ -N csstp_${sample}"   >> "${main_dir}/csstp_${sample}_${start_time}.sh"
					echo -e "#$ -cwd"  >> "${main_dir}/csstp_${sample}_${start_time}.sh"
					echo -e "#$ -q short.q\n"  >> "${main_dir}/csstp_${sample}_${start_time}.sh"
					echo -e "module load Python/3.6.1\n" >> "${main_dir}/csstp_${sample}_${start_time}.sh"
					# Defaulting to gapped/98, change if you want to include user preferences
					echo -e "\"${shareScript}/run_c-sstar_on_single_alternate_DB.sh\" \"${sample}\" g o \"${project}\" \"${3}\" \"--plasmid\"" >> "${main_dir}/csstp_${sample}_${start_time}.sh"
					echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_csstarp_complete.txt\"" >> "${main_dir}/csstp_${sample}_${start_time}.sh"
					#if [[ "${counter}" -lt "${last_index}" ]]; then
						qsub "${main_dir}/csstp_${sample}_${start_time}.sh"
					#else
					#	qsub -sync y "${main_dir}/csstp_${sample}_${start_time}.sh"
					#fi
					if [[ -f "${processed}/${project}/${sample}/c-sstar_plasmid/${sample}.${alt_database}.gapped_98_sstar_summary.txt" ]]; then
						rm "${processed}/${project}/${sample}/c-sstar_plasmid/${sample}.${alt_database}.gapped_98_sstar_summary.txt"
					fi
				else
					echo "${project}/${sample} already has 0608 PLASMID"
					echo "$(date)" > "${main_dir}/complete/${sample}_csstarp_complete.txt"
				fi
			else
				echo "${project}/${sample} doesnt have a plasmid folder, so no further actions required"
				echo "$(date)" > "${main_dir}/complete/${sample}_csstarp_complete.txt"
			fi
		# Counter has exceeded max submissions limit
		else
			waiting_for_index=$(( counter - max_subs ))
			waiting_sample=$(echo "${arr[${waiting_for_index}]}" | cut -d'/' -f2)
			timer=0
			echo "Index is above max submissions, waiting for index ${waiting_for_index}:${waiting_sample} to complete"
			while :
			do
				# Check to see if max time limit has been reached, exit if so.
				if [[ ${timer} -gt 1800 ]]; then
					echo "Timer exceeded limit of 1800 seconds 30 minutes"
					break
				fi
				# Check if the sample it is waiting on has been completed or if it has a usable assembly
				if [[ -f "${main_dir}/complete/${waiting_sample}_csstarn_complete.txt" ]] || [[ ! -s "${processed}/${project}/${waiting_sample}/Assembly/${waiting_sample}_scaffolds_trimmed.fasta" ]]; then
					# Check if ol results exist and skip if so
					if [[ ! -f "${processed}/${project}/${sample}/c-sstar/${sample}.${alt_database}.gapped_98_sstar_summary.txt" ]]; then
						echo  "Index is below max submissions, submitting"
						echo -e "#!/bin/bash -l\n" > "${main_dir}/csstn_${sample}_${start_time}.sh"
						echo -e "#$ -o csstn_${sample}.out" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
						echo -e "#$ -e csstn_${sample}.err" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
						echo -e "#$ -N csstn_${sample}"   >> "${main_dir}/csstn_${sample}_${start_time}.sh"
						echo -e "#$ -cwd"  >> "${main_dir}/csstn_${sample}_${start_time}.sh"
						echo -e "#$ -q short.q\n"  >> "${main_dir}/csstn_${sample}_${start_time}.sh"
						echo -e "module load Python/3.6.1\n" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
						# Defaulting to gapped/98, change if you want to include user preferences
						echo -e "\"${shareScript}/run_c-sstar_on_single_alternate_DB.sh\" \"${sample}\" g h \"${project}\" \"${3}\"" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
						echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_csstarn_complete.txt\"" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
						#if [[ "${counter}" -lt "${last_index}" ]]; then
							qsub "${main_dir}/csstn_${sample}_${start_time}.sh"
						#else
						#	if [[ -d "${processed}/${project}/${sample}/c-sstar_plasmid" ]]; then
						#		qsub "${main_dir}/csstn_${sample}_${start_time}.sh"
						#	else
						#		qsub -sync y "${main_dir}/csstn_${sample}_${start_time}.sh"
						#	fi
						#fi
					else
						echo "${project}/${sample} already has 0608"
						echo "$(date)" > "${main_dir}/complete/${sample}_csstarn_complete.txt"
					fi
					# Check if c-sstar plasmid folder exists
					if [[ -d "${processed}/${project}/${sample}/c-sstar_plasmid" ]]; then
						# Check if data already exists and skip if so
						if [[ ! -f "${processed}/${project}/${sample}/c-sstar_plasmid/${sample}.${alt_database}.gapped_40_sstar_summary.txt" ]]; then
							echo -e "#!/bin/bash -l\n" > "${main_dir}/csstp_${sample}_${start_time}.sh"
							echo -e "#$ -o csstp_${sample}.out" >> "${main_dir}/csstp_${sample}_${start_time}.sh"
							echo -e "#$ -e csstp_${sample}.err" >> "${main_dir}/csstp_${sample}_${start_time}.sh"
							echo -e "#$ -N csstp_${sample}"   >> "${main_dir}/csstp_${sample}_${start_time}.sh"
							echo -e "#$ -cwd"  >> "${main_dir}/csstp_${sample}_${start_time}.sh"
							echo -e "#$ -q short.q\n"  >> "${main_dir}/csstp_${sample}_${start_time}.sh"
							echo -e "module load Python/3.6.1\n" >> "${main_dir}/csstp_${sample}_${start_time}.sh"
							# Defaulting to gapped/98, change if you want to include user preferences
							echo -e "\"${shareScript}/run_c-sstar_on_single_alternate_DB.sh\" \"${sample}\" g o \"${project}\" \"${3}\" \"--plasmid\"" >> "${main_dir}/csstp_${sample}_${start_time}.sh"
							echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_csstarp_complete.txt\"" >> "${main_dir}/csstp_${sample}_${start_time}.sh"
							#if [[ "${counter}" -lt "${last_index}" ]]; then
								qsub "${main_dir}/csstp_${sample}_${start_time}.sh"
							#else
							#	qsub -sync y "${main_dir}/csstp_${sample}_${start_time}.sh"
							#fi
							if [[ -f "${processed}/${project}/${sample}/c-sstar_plasmid/${sample}.${alt_database}.gapped_98_sstar_summary.txt" ]]; then
								rm "${processed}/${project}/${sample}/c-sstar_plasmid/${sample}.${alt_database}.gapped_98_sstar_summary.txt"
							fi
						else
							echo "${project}/${sample} already has 0608 PLASMID"
							echo "$(date)" > "${main_dir}/complete/${sample}_csstarp_complete.txt"
						fi
					else
						echo "${project}/${sample} doesnt have a plasmid folder, so no further actions required"
						echo "$(date)" > "${main_dir}/complete/${sample}_csstarp_complete.txt"
					fi
					break
				# Wait 5 seconds and check if waiting sample is complete yet
				else
					timer=$(( timer + 5 ))
					echo "sleeping for 5 seconds, so far slept for ${timer}"
					sleep 5
				fi
			done
		fi
	# If no assembly file exists
	else
		echo "${project}/${sample} does not have an assembly file???"
	fi
	counter=$(( counter + 1 ))
done

# Loop to ensure all samples are complete (or time runs) before allowing the script to exit
timer=0
ptimer=0
for item in "${arr[@]}"; do
	waiting_sample=$(echo "${item}" | cut -d'/' -f2)
	if [[ -f "${main_dir}/complete/${waiting_sample}_csstarn_complete.txt" ]] || [[ ! -s "${processed}/${project}/${waiting_sample}/Assembly/${waiting_sample}_scaffolds_trimmed.fasta" ]]; then
		echo "${item} is complete normal"
		if [[ -f "${shareScript}/csstp_${sample}.out" ]]; then
			mv "${shareScript}/csstp_${sample}.out" "${main_dir}"
		fi
		if [[ -f "${shareScript}/csstp_${sample}.err" ]]; then
			mv "${shareScript}/csstp_${sample}.err" "${main_dir}"
		fi
		# Check if plasmid csstar is complete also and wait a total of 30 minutes for all samples to be checked
		if [[ -f "${main_dir}/complete/${waiting_sample}_csstarp_complete.txt" ]] || [[ ! -s "${processed}/${project}/${waiting_sample}/plasmidAssembly/${waiting_sample}_plasmid_scaffolds_trimmed.fasta" ]]; then
			while :
			do
					if [[ ${ptimer} -gt 1800 ]]; then
						echo "Timer exceeded limit of 1800 seconds = 30 minutes"
						exit 1
					fi
					if [[ -f "${main_dir}/complete/${waiting_sample}_csstarp_complete.txt" ]]; then
						echo "${item} is complete plasmid"
						if [[ -f "${shareScript}/csstp_${sample}.out" ]]; then
							mv "${shareScript}/csstp_${sample}.out" "${main_dir}"
						fi
						if [[ -f "${shareScript}/csstp_${sample}.err" ]]; then
							mv "${shareScript}/csstp_${sample}.err" "${main_dir}"
						fi
						break
					else
						ptimer=$(( ptimer + 5 ))
						echo "sleeping for 5 seconds, so far slept for ${ptimer}"
						sleep 5
					fi
			done
		fi
	else
		# Check every 5 seconds to see if the sample has completed normal csstar analysis
		while :
		do
				if [[ ${timer} -gt 1800 ]]; then
					echo "Timer exceeded limit of 3600 seconds = 30 minutes"
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
