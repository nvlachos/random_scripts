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
shareScript="/scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/scripts/Quaisar-H"
#
# Usage ./act_by_list.sh list_name(currently has to be placed in /scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR folder) description of list function
#
# script changes depending on what needs to be run through the list
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to act_by_list.sh, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./abl_mass_qsub_template.sh path_to_list_file(single sample ID per line, e.g. B8VHY/1700128 (it must include project id also)) max_submissions"
	echo "Output location varies depending on which tasks are performed but will be found somewhere under ${processed}"
	exit 0
fi

arr=()
while IFS= read -r line || [[ "$line" ]];  do
echo "${line}"
	line=$(echo "${line}" | tr -d [:space:] )
	arr+=("${line}")
done < ${1}

arr_size="${#arr[@]}"
echo "-${arr_size}:${arr[@]}-"

if [[ ! -z "${2}" ]]; then
	max_subs="${2}"
else
	max_subs=10
fi

# Loop through and act on each sample name in the passed/provided list
counter=0
max_subs=${2}
main_dir="${share}/mass_subs/ANI_subs"
if [[ ! -d "${share}/mass_subs/ANI_subs" ]]; then
	mkdir "${share}/mass_subs/ANI_subs"
	mkdir "${share}/mass_subs/ANI_subs/complete"
elif [[ ! -d "${share}/mass_subs/ANI_subs/complete" ]]; then
	mkdir "${share}/mass_subs/ANI_subs/complete"
fi

start_time=$(date)

while [ ${counter} -lt ${arr_size} ] ; do
	sample=$(echo "${arr[${counter}]}" | cut -d'/' -f2)
	project=$(echo "${arr[${counter}]}" | cut -d'/' -f1)
	#genus=$(echo "${arr[${counter}]}" | cut -d'/' -f3 | cut -d'_' -f1)
	#species=$(echo "${arr[${counter}]}" | cut -d'/' -f3 | cut -d'_' -f1)
	echo "${sample} and ${project}:"
	echo "${counter}-${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed.fasta"
	if [[ -s "${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed.fasta" ]]; then
		rm -r "${processed}/${project}/${sample}/ANI/"
		genus="Acinetobacter"
		species="baumannii"
	# if [[ -s "${processed}/${project}/${sample}/${sample}.tax" ]]; then
	# 	while IFS= read -r line;
	# 	do
	# 		# Grab first letter of line (indicating taxonomic level)
	# 		first=${line::1}
	# 		# Assign taxonomic level value from 4th value in line (1st-classification level,2nd-% by kraken, 3rd-true % of total reads, 4th-identifier)
	# 		if [ "${first}" = "S" ]
	# 		then
	# 			species=$(echo "${line}" | awk -F ' ' '{print $4}')
	# 		elif [ "${first}" = "G" ]
	# 		then
	# 			genus=$(echo "${line}" | awk -F ' ' '{print $4}')
	# 		#elif [ "${first}" = "F" ]
	# 		#then
	# 		#	family=$(echo "${line}" | awk -F ' ' '{print $4}')
	# 		#elif [ "${first}" = "O" ]
	# 		#then
	# 		#	order=$(echo "${line}" | awk -F ' ' '{print $4}')
	# 		#elif [ "${first}" = "C" ]
	# 		#then
	# 		#	class=$(echo "${line}" | awk -F ' ' '{print $4}')
	# 		#elif [ "${first}" = "P" ]
	# 		#then
	# 		#	phylum=$(echo "${line}" | awk -F ' ' '{print $4}')
	# 		#elif [ "${first}" = "K" ]
	# 		#then
	# 		#	kingdom=$(echo "${line}" | awk -F ' ' '{print $4}')
	# 		#elif [ "${first}" = "D" ]
	# 		#then
	# 		#	domain=$(echo "${line}" | awk -F ' ' '{print $4}')
	# 		fi
	# 	done < "${processed}/${project}/${sample}/${sample}.tax"
	 	if [[ ${counter} -lt ${max_subs} ]]; then
			if [[ ! -f "${processed}/${sample}/ANI/best_ANI_hits_ordered(${sample}_vs_${genus,})" ]]; then
				echo  "Index is below max submissions, submitting"
				echo -e "#!/bin/bash -l\n" > "${main_dir}/ani_${sample}_${start_time}.sh"
				echo -e "#$ -o ani_${sample}.out" >> "${main_dir}/ani_${sample}_${start_time}.sh"
				echo -e "#$ -e ani_${sample}.err" >> "${main_dir}/ani_${sample}_${start_time}.sh"
				echo -e "#$ -N ani_${sample}"   >> "${main_dir}/ani_${sample}_${start_time}.sh"
				echo -e "#$ -cwd"  >> "${main_dir}/ani_${sample}_${start_time}.sh"
				echo -e "#$ -q short.q\n"  >> "${main_dir}/ani_${sample}_${start_time}.sh"
				# Defaulting to gapped/98, change if you want to include user preferences
				echo -e "\"${shareScript}/run_ANI.sh\" \"${sample}\" \"${genus}\" \"${species}\" \"${project}\"" >> "${main_dir}/ani_${sample}_${start_time}.sh"
				echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_ani_complete.txt\"" >> "${main_dir}/ani_${sample}_${start_time}.sh"
				qsub "${main_dir}/ani_${sample}_${start_time}.sh"
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
						echo  "Index is below max submissions, submitting"
						echo -e "#!/bin/bash -l\n" > "${main_dir}/ani_${sample}_${start_time}.sh"
						echo -e "#$ -o ani_${sample}.out" >> "${main_dir}/ani_${sample}_${start_time}.sh"
						echo -e "#$ -e ani_${sample}.err" >> "${main_dir}/ani_${sample}_${start_time}.sh"
						echo -e "#$ -N ani_${sample}"   >> "${main_dir}/ani_${sample}_${start_time}.sh"
						echo -e "#$ -cwd"  >> "${main_dir}/ani_${sample}_${start_time}.sh"
						echo -e "#$ -q short.q\n"  >> "${main_dir}/ani_${sample}_${start_time}.sh"
						# Defaulting to gapped/98, change if you want to include user preferences
						echo -e "\"${shareScript}/run_ANI.sh\" \"${sample}\" \"${genus}\" \"${species}\" \"${project}\"" >> "${main_dir}/ani_${sample}_${start_time}.sh"
						echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_ani_complete.txt\"" >> "${main_dir}/ani_${sample}_${start_time}.sh"
						qsub "${main_dir}/ani_${sample}_${start_time}.sh"
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
global_end_time=$(date "+%m-%d-%Y @ %Hh_%Mm_%Ss")
#Script exited gracefully (unless something else inside failed)
printf "%s %s" "Act_by_list.sh has completed ${2}" "${global_end_time}" | mail -s "act_by_list complete" nvx4@cdc.gov
exit 0
