#!/bin/sh -l

#$ -o ablmq-anim.out
#$ -e ablmq-anim.err
#$ -N ablmq-anim
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
# Usage ./abl_mass_qsub_ANIm.sh path_to_list max_concurrent_submissions output_directory_for_scripts clobberness[keep|clobber]
#

# Number regex to test max concurrent submission parametr
number='^[0-9]+$'

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./abl_mass_qsub_ANIm.sh path_to_list_file(single sample ID per line, e.g. B8VHY/1700128 (it must include project id also)) max_concurrent_submissions output_directory_for_scripts clobberness[keep|clobber]"
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
main_dir="${3}/ANIm_scaffold_subs"

if [[ ! -d "${main_dir}" ]]; then
	mkdir -p "${main_dir}/complete"
elif [[ ! -d "${main_dir}/complete" ]]; then
	mkdir "${main_dir}/complete"
fi
cp ./config.sh ${main_dir}

start_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")

# Create and submit qsub scripts to get ANI for all isolates
while [ ${counter} -lt ${arr_size} ] ; do
	sample=$(echo "${arr[${counter}]}" | cut -d'/' -f2 | cut -d':' -f1)
	#sample=${sample:0:20}
	project=$(echo "${arr[${counter}]}" | cut -d'/' -f1)
	#genus=${5}
	#species=${6}



	if [[ "${clobberness}" == "clobber" ]]; then
		echo "Trying to remove ${processed}/${project}/${sample}/ANI"
		rm -r "${processed}/${project}/${sample}/ANI"
		rm -r "${processed}/${project}/${sample}/${sample}.tax"
		"${shareScript}/determine_taxID.sh" "${sample}" "${project}"
	fi
	#echo "${sample} and ${project}:"
	#echo "${counter}-${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed.fasta"

	# Ensure tax file exists to get proper DB to run ANI against
	if [[ -s "${processed}/${project}/${sample}/${sample}.tax" ]]; then
		# Parse tax file
		while IFS= read -r line; do
			# Grab first letter of line (indicating taxonomic level)
			first=${line:0:1}
			# Assign taxonomic level value from 4th value in line (1st-classification level, 2nd-% by kraken, 3rd-true % of total reads, 4th-identifier)
			if [ "${first}" = "s" ]
			then
				species=$(echo "${line}" | awk -F ' ' '{print $2}')
			elif [ "${first}" = "G" ]
			then
				genus=$(echo "${line}" | awk -F ' ' '{print $2}')
				# Only until ANI gets fixed
				if [[ ${genus} == "Clostridioides" ]]; then
					genus="Clostridium"
				fi
			fi
		done < "${processed}/${project}/${sample}/${sample}.tax"


		sleep 15
		#Temp assignment, if specific DB is necessary
		#genus="Acinetobacter"
		#species="baumannii"

		echo "G-${genus},S-${species}"

		# Check if counter is below max submission limit
	 	if [[ ${counter} -lt ${max_subs} ]]; then
			# Check if old data exists, skip if so
			if [[ -s "${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed.fasta" ]]; then
				if [[ ! -f "${processed}/${sample}/ANI/best_anim_hits_ordered(${sample}_vs_${genus,})" ]]; then
					echo  "Index is below max submissions, submitting"
					echo "Going to make ${main_dir}/anim_${sample}_${start_time}.sh"
					echo -e "#!/bin/bash -l\n" > "${main_dir}/anim_${sample}_${start_time}.sh"
					echo -e "#$ -o anim_${sample}.out" >> "${main_dir}/anim_${sample}_${start_time}.sh"
					echo -e "#$ -e anim_${sample}.err" >> "${main_dir}/anim_${sample}_${start_time}.sh"
					echo -e "#$ -N anim_${sample}"   >> "${main_dir}/anim_${sample}_${start_time}.sh"
					echo -e "#$ -cwd"  >> "${main_dir}/anim_${sample}_${start_time}.sh"
					echo -e "#$ -q short.q\n"  >> "${main_dir}/anim_${sample}_${start_time}.sh"
					echo -e "cd ${shareScript}" >> "${main_dir}/anim_${sample}_${start_time}.sh"
					echo -e "\"${shareScript}/run_ANIm.sh\" \"${sample}\" \"${genus}\" \"${species}\" \"${project}\"" >> "${main_dir}/anim_${sample}_${start_time}.sh"
					echo -e "\"${shareScript}/determine_taxID.sh\" \"${sample}\" \"${project}\"" >> "${main_dir}/anim_${sample}_${start_time}.sh"
					echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_anim_complete.txt\"" >> "${main_dir}/anim_${sample}_${start_time}.sh"
					cd "${main_dir}"
					#if [[ "${counter}" -lt "${last_index}" ]]; then
						qsub "${main_dir}/anim_${sample}_${start_time}.sh"
					#else
					#	qsub -sync y "${main_dir}/anim_${sample}_${start_time}.sh"
					#fi
				# Old data exists, skipping
				else
					echo "${project}/${sample} already has ANI summary"
					echo "$(date)" > "${main_dir}/complete/${sample}_anim_complete.txt"
				fi
			# No Assembly file to run ANI on, skipping
			else
				echo "${project}/${sample} does not have assembly"
				echo "$(date)" > "${main_dir}/complete/${sample}_anim_complete.txt"
			fi
		# Counter is above limit, wait until "slot" opens up"
		else
			waiting_for_index=$(( counter - max_subs ))
			waiting_sample=$(echo "${arr[${waiting_for_index}]}" | cut -d'/' -f2)
			timer=0
			echo "Index is above max submissions, waiting for index ${waiting_for_index}:${waiting_sample} to complete"
			# Loop to wait until "waiting" sample is complete
			while :
			do
				# If timer exceeeds limit then exit
				if [[ ${timer} -gt 1800 ]]; then
					echo "Timer exceeded limit of 1800 seconds 30 minutes"
					break
				fi
				# Check if "waiting" sample is complete
				if [[ -f "${main_dir}/complete/${waiting_sample}_anim_complete.txt" ]]; then
					# Check if an assembly exists to run ANI on
					if [[ -s "${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed.fasta" ]]; then
						# Check if old data exists, skip if so
						if [[ ! -f "${processed}/${project}/${sample}/ANI/best_anim_hits_ordered(${sample}_vs_${genus,})" ]]; then
							echo  "${waiting_sample}(${waiting_for_index}) is not complete, submitting ${sample} ($counter)"
							echo -e "#!/bin/bash -l\n" > "${main_dir}/anim_${sample}_${start_time}.sh"
							echo -e "#$ -o anim_${sample}.out" >> "${main_dir}/anim_${sample}_${start_time}.sh"
							echo -e "#$ -e anim_${sample}.err" >> "${main_dir}/anim_${sample}_${start_time}.sh"
							echo -e "#$ -N anim_${sample}"   >> "${main_dir}/anim_${sample}_${start_time}.sh"
							echo -e "#$ -cwd"  >> "${main_dir}/anim_${sample}_${start_time}.sh"
							echo -e "#$ -q short.q\n"  >> "${main_dir}/anim_${sample}_${start_time}.sh"
							echo -e "cd ${shareScript}" >> "${main_dir}/anim_${sample}_${start_time}.sh"
							echo -e "\"${shareScript}/run_ANIm.sh\" \"${sample}\" \"${genus}\" \"${species}\" \"${project}\"" >> "${main_dir}/anim_${sample}_${start_time}.sh"
							echo -e "\"${shareScript}/determine_taxID.sh\" \"${sample}\" \"${project}\"" >> "${main_dir}/anim_${sample}_${start_time}.sh"
							echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_anim_complete.txt\"" >> "${main_dir}/anim_${sample}_${start_time}.sh"
							cd "${main_dir}"
							#if [[ "${counter}" -lt "${last_index}" ]]; then
								qsub "${main_dir}/anim_${sample}_${start_time}.sh"
							#else
							#	qsub -sync y "${main_dir}/anim_${sample}_${start_time}.sh"
							#fi
						# Old data exists, skipping
						else
							echo "${project}/${sample} already has ANI summary"
							echo "$(date)" > "${main_dir}/complete/${sample}_anim_complete.txt"
						fi
						# No Assembly file to run ANI on, skipping
					else
						echo "${project}/${sample} does not have assembly"
						echo "$(date)" > "${main_dir}/complete/${sample}_anim_complete.txt"
					fi
					break
				# Wait 5 seconds before checking if "waiting" sample is complete
				else
					timer=$(( timer + 5 ))
					echo "sleeping for 5 seconds, so far slept for ${timer}"
					sleep 5
				fi
			done
		fi
	else
		echo "No TAX file available to use for ANI determination"
	fi
	counter=$(( counter + 1 ))
done

# Loop to ensure all samples are complete (or time runs) before allowing the script to exit
timer=0
for item in "${arr[@]}"; do
	waiting_sample=$(echo "${item}" | cut -d'/' -f2)
	if [[ -f "${main_dir}/complete/${waiting_sample}_anim_complete.txt" ]]; then
		echo "${item} is complete"
		if [[ -f "${shareScript}/anim_${waiting_sample}.out" ]]; then
			mv "${shareScript}/anim_${waiting_sample}.out" ${main_dir}
		fi
		if [[ -f "${shareScript}/anim_${waiting_sample}.err" ]]; then
			mv "${shareScript}/anim_${waiting_sample}.err" ${main_dir}
		fi
	else
		# Check every 5 seconds to see if the sample has completed normal csstar analysis
		while :
		do
				if [[ ${timer} -gt 3600 ]]; then
					echo "Timer exceeded limit of 3600 seconds = 60 minutes"
					exit 1
				fi
				if [[ -f "${main_dir}/complete/${waiting_sample}_anim_complete.txt" ]]; then
					echo "${item} is complete"
					if [[ -f "${shareScript}/anim_$waiting_{sample}.out" ]]; then
						mv "${shareScript}/anim_${waiting_sample}.out" ${main_dir}
					fi
					if [[ -f "${shareScript}/anim_${waiting_sample}.err" ]]; then
						mv "${shareScript}/anim_${waiting_sample}.err" ${main_dir}
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
printf "%s %s" "abl_mass_qsub_ANIm.sh has completed" "${global_end_time}" | mail -s "abl_mass_qsub_ANIm.sh complete" nvx4@cdc.gov

#Script exited gracefully (unless something else inside failed)
exit 0
