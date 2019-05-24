#!/bin/sh -l

#$ -o ablmq-bus.out
#$ -e ablmq-bus.err
#$ -N ablmq-bus
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
# Usage ./abl_mass_qsub_BUSCO.sh path_to_list max_concurrent_submissions output_directory_for_scripts clobberness[keep|clobber]
#

# Number regex to test max concurrent submission parametr
number='^[0-9]+$'

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./abl_mass_qsub_BUSCO.sh path_to_list_file(single sample ID per line, e.g. B8VHY/1700128 (it must include project id also)) max_concurrent_submissions output_directory_for_scripts"
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
main_dir="${3}/busco_subs"
if [[ ! -d "${3}/busco_subs" ]]; then
	mkdir "${3}/busco_subs"
	mkdir "${3}/busco_subs/complete"
elif [[ ! -d "${3}/busco_subs/complete" ]]; then
	mkdir "${3}/busco_subs/complete"
fi

start_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")

# Create and submit qsub scripts to do BUSCO analysis on all isolates on the list
while [ ${counter} -lt ${arr_size} ] ; do
	sample=$(echo "${arr[${counter}]}" | cut -d'/' -f2)
	project=$(echo "${arr[${counter}]}" | cut -d'/' -f1)
	if [[ "${clobberness}" == "clobber" ]]; then
		rm -r "${processed}/${project}/${sample}/BUSCO"
	fi
	# Get taxonomy from tax file
	if [[ ! -f "${processed}/${project}/${sample}/${sample}.tax" ]]; then
		"${shareScript}/determine_taxID.sh" "${sample}" "${project}"
	fi
	if [[ -s "${processed}/${project}/${sample}/${sample}.tax" ]]; then
		while IFS= read -r line;
		do
			# Grab first letter of line (indicating taxonomic level)
			first=${line::1}
			# Assign taxonomic level value from 4th value in line (1st-classification level,2nd-% by kraken, 3rd-true % of total reads, 4th-identifier)
			if [ "${first}" = "s" ]
			then
				species=$(echo "${line}" | awk -F ' ' '{print $2}')
			elif [ "${first}" = "G" ]
			then
				genus=$(echo "${line}" | awk -F ' ' '{print $2}')
			elif [ "${first}" = "F" ]
			then
				family=$(echo "${line}" | awk -F ' ' '{print $2}')
			elif [ "${first}" = "O" ]
			then
				order=$(echo "${line}" | awk -F ' ' '{print $2}')
			elif [ "${first}" = "C" ]
			then
				class=$(echo "${line}" | awk -F ' ' '{print $2}')
			elif [ "${first}" = "P" ]
			then
				phylum=$(echo "${line}" | awk -F ' ' '{print $2}')
			elif [ "${first}" = "K" ]
			then
				kingdom=$(echo "${line}" | awk -F ' ' '{print $2}')
			elif [ "${first}" = "D" ]
			then
				domain=$(echo "${line}" | awk -F ' ' '{print $2}')
			fi
		done < "${processed}/${project}/${sample}/${sample}.tax"
		buscoDB="bacteria_odb9"
		# Iterate through taxon levels (species to domain) and test if a match occurs to entry in database. If so, compare against it
		for tax in $species $genus $family $order $class $phylum $kingdom $domain
		do
			echo "$tax"
			if [ -d "${local_DBs}/BUSCO/${tax,}_odb9" ]
			then
				buscoDB="${tax}_odb9"
				break
			fi
		done
		echo "Using $buscoDB fo ${sample}"
		# Check if prokka output exists and is usable
		if [[ -s "${processed}/${project}/${sample}/prokka/${sample}_PROKKA.gbf" ]] || [[ -s "${processed}/${project}/${sample}/prokka/${sample}_PROKKA.gff" ]]; then
			# Check if counter is below max sub limit
			if [[ ${counter} -lt ${max_subs} ]]; then
				# Check if old data exists, skip if so
				if [[ ! -f "${processed}/${project}/${sample}/BUSCO/short_summary_${sample}.txt" ]]; then
					echo  "Index is below max submissions, submitting"
					echo -e "#!/bin/bash -l\n" > "${main_dir}/busco_${sample}_${start_time}.sh"
					echo -e "#$ -o busco_${sample}.out" >> "${main_dir}/busco_${sample}_${start_time}.sh"
					echo -e "#$ -e busco_${sample}.err" >> "${main_dir}/busco_${sample}_${start_time}.sh"
					echo -e "#$ -N busco_${sample}"   >> "${main_dir}/busco_${sample}_${start_time}.sh"
					echo -e "#$ -cwd"  >> "${main_dir}/busco_${sample}_${start_time}.sh"
					echo -e "#$ -q short.q\n"  >> "${main_dir}/busco_${sample}_${start_time}.sh"
					echo -e "module load busco/3.0.1" >> "${main_dir}/busco_${sample}_${start_time}.sh"
					echo -e "cd ${shareScript}" >> "${main_dir}/busco_${sample}_${start_time}.sh"
					echo -e "\"${shareScript}/do_busco.sh\" \"${sample}\" \"${buscoDB}\" \"${project}\"" >> "${main_dir}/busco_${sample}_${start_time}.sh"
					echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_busco_complete.txt\"" >> "${main_dir}/busco_${sample}_${start_time}.sh"
					cd "${main_dir}"
					#if [[ "${counter}" -lt "${last_index}" ]]; then
						qsub "${main_dir}/busco_${sample}_${start_time}.sh"
					#else
					#	qsub -sync y "${main_dir}/busco_${sample}_${start_time}.sh"
					#fi
				# Old data exists, skipping
				else
					echo "${project}/${sample} already has busco summary"
					echo "$(date)" > "${main_dir}/complete/${sample}_busco_complete.txt"
				fi
			# Counter is ablove max subs, must wait for slot to open
			else
				waiting_for_index=$(( counter - max_subs ))
				waiting_sample=$(echo "${arr[${waiting_for_index}]}" | cut -d'/' -f2)
				timer=0
				echo "Index is above max submissions, waiting for index ${waiting_for_index}:${waiting_sample} to complete"
				# Loop to wait for :waiting: sample to complete
				while :
				do
						# Check if timer is above max limit
					if [[ ${timer} -gt 1800 ]]; then
						echo "Timer exceeded limit of 1800 seconds 30 minutes"
						break
					fi
					# Check if "waiting" sample is complete
					if [[ -f "${main_dir}/complete/${waiting_sample}_busco_complete.txt" ]]; then
						# Check if old results exist
						if [[ ! -f "${processed}/${project}/${sample}/BUSCO/short_summary_${sample}.txt" ]]; then
							echo  "Index is below max submissions, submitting"
							echo -e "#!/bin/bash -l\n" > "${main_dir}/busco_${sample}_${start_time}.sh"
							echo -e "#$ -o busco_${sample}.out" >> "${main_dir}/busco_${sample}_${start_time}.sh"
							echo -e "#$ -e busco_${sample}.err" >> "${main_dir}/busco_${sample}_${start_time}.sh"
							echo -e "#$ -N busco_${sample}"   >> "${main_dir}/busco_${sample}_${start_time}.sh"
							echo -e "#$ -cwd"  >> "${main_dir}/busco_${sample}_${start_time}.sh"
							echo -e "#$ -q short.q\n"  >> "${main_dir}/busco_${sample}_${start_time}.sh"
							echo -e "module load busco/3.0.1" >> "${main_dir}/busco_${sample}_${start_time}.sh"
							echo -e "cd ${shareScript}" >> "${main_dir}/busco_${sample}_${start_time}.sh"
							echo -e "\"${shareScript}/do_busco.sh\" \"${sample}\" \"${buscoDB}\" \"${project}\"" >> "${main_dir}/busco_${sample}_${start_time}.sh"
							echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_busco_complete.txt\"" >> "${main_dir}/busco_${sample}_${start_time}.sh"
							cd "${main_dir}"
							#if [[ "${counter}" -lt "${last_index}" ]]; then
								qsub "${main_dir}/busco_${sample}_${start_time}.sh"
							#else
							#	qsub -sync y "${main_dir}/busco_${sample}_${start_time}.sh"
							#fi
						# Old data exists, skipping
						else
							echo "${project}/${sample} already has busco summary"
							echo "$(date)" > "${main_dir}/complete/${sample}_busco_complete.txt"
						fi
						break
					# Wait 5 seconds before checking if "waiting" sample is done again
					else
						timer=$(( timer + 5 ))
						echo "sleeping for 5 seconds, so far slept for ${timer}"
						sleep 5
					fi
				done
			fi
		# No prokka output to run BUSCO on
		else
			echo "${project}/${sample} does not have prokka output to ru nbusco on"
			echo "$(date)" > "${main_dir}/complete/${sample}_busco_complete.txt"
		fi
	counter=$(( counter + 1 ))
done

# Loop to ensure all samples are complete (or time runs) before allowing the script to exit
timer=0
for item in "${arr[@]}"; do
	waiting_sample=$(echo "${item}" | cut -d'/' -f2)
	if [[ -f "${main_dir}/complete/${waiting_sample}_busco_complete.txt" ]]; then
		echo "${item} is complete"
		if [[ -f "${shareScript}/busco_${sample}.out" ]]; then
			mv "${shareScript}/busco_${sample}.out" ${main_dir}
		fi
		if [[ -f "${shareScript}/busco_${sample}.err" ]]; then
			mv "${shareScript}/busco_${sample}.err" ${main_dir}
		fi
	else
		# Check every 5 seconds to see if the sample has completed normal csstar analysis
		while :
		do
				if [[ ${timer} -gt 1800 ]]; then
					echo "Timer exceeded limit of 1800 seconds = 30 minutes"
					exit 1
				fi
				if [[ -f "${main_dir}/complete/${waiting_sample}_busco_complete.txt" ]]; then
					echo "${item} is complete"
					if [[ -f "${shareScript}/busco_${sample}.out" ]]; then
						mv "${shareScript}/busco_${sample}.out" ${main_dir}
					fi
					if [[ -f "${shareScript}/busco_${sample}.err" ]]; then
						mv "${shareScript}/busco_${sample}.err" ${main_dir}
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
printf "%s %s" "abl_mass_qsub_BUSCO.sh has completed" "${global_end_time}" | mail -s "abl_mass_qsub_BUSCO.sh complete" nvx4@cdc.gov

#Script exited gracefully (unless something else inside failed)
exit 0
