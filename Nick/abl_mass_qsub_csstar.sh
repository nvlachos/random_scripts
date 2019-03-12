#!/bin/sh -l

#$ -o ablmq-cs.out
#$ -e ablmq-cs.err
#$ -N ablmq-cs
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
# Import the config file with shortcuts and settings
pwd
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi
. ./config.sh
#Import the module file that loads all necessary mods
. "${mod_changers}/pipeline_mods"

#
# Usage ./abl_mass_qsub_csstar.sh path_to_list max_concurrent_submissions
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to ./abl_mass_qsub_csstar.sh, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./abl_mass_qsub_csstar.sh path_to_list_file(single sample ID per line, e.g. B8VHY/1700128 (it must include project id also)) max_submissions output_directory_for_scripts"
	exit 1
elif [[ ! -f "${1}" ]]; then
	echo "${1} (list) does not exist...exiting"
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


# Create direcory to hold all temporary qsub scripts
counter=0
max_subs=${2}

# format name being extracted from alt database
main_dir="${3}/csstar_subs"
#cp ./config ${main_dir}
if [[ ! -d "${3}/csstar_subs" ]]; then
	mkdir -p "${3}/csstar_subs/complete"
elif [[ ! -d  "${3}/csstar_subs/complete" ]]; then
	mkdir "${3}/csstar_subs/complete"
fi

start_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")

# Create and submit scripts to run default csstar on all samples on the list
while [ ${counter} -lt ${arr_size} ] ; do
	sample=$(echo "${arr[${counter}]}" | cut -d'/' -f2)
	project=$(echo "${arr[${counter}]}" | cut -d'/' -f1)
	echo ${counter}-${project}-${sample}
	if [[ -s "${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed.fasta" ]]; then
		echo "Test"
		if [[ ${counter} -lt ${max_subs} ]]; then
			rm "${processed}/${project}/${sample}/c-sstar/${sample}.${resGANNOT_srst2_filename}.gapped_98_sstar_summary.txt"
			if [[ ! -f "${processed}/${project}/${sample}/c-sstar/${sample}.${resGANNOT_srst2_filename}.gapped_98_sstar_summary.txt" ]]; then
				echo  "Index is below max submissions, submitting"
				echo -e "#!/bin/bash -l\n" > "${main_dir}/csstn_${sample}_${start_time}.sh"
				echo -e "#$ -o csstn_${sample}.out" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
				echo -e "#$ -e csstn_${sample}.err" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
				echo -e "#$ -N csstn_${sample}"   >> "${main_dir}/csstn_${sample}_${start_time}.sh"
				echo -e "#$ -cwd"  >> "${main_dir}/csstn_${sample}_${start_time}.sh"
				echo -e "#$ -q short.q\n"  >> "${main_dir}/csstn_${sample}_${start_time}.sh"
				echo -e "module load Python/3.6.1\n" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
				# Defaulting to gapped/98, change if you want to include user preferences
				echo -e "cd ${shareScript}" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
				echo -e "\"${shareScript}/run_c-sstar_on_single.sh\" \"${sample}\" g h \"${project}\"" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
				echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_csstarn_complete.txt\"" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
				cd "${main_dir}"
				echo "submitting ${main_dir}/csstn_${sample}_${start_time}.sh"
				if [[ "${counter}" -lt "${last_index}" ]]; then
					qsub "${main_dir}/csstn_${sample}_${start_time}.sh"
				else
					if [[ -d "${processed}/${project}/${sample}/c-sstar_plasmid" ]]; then
						qsub "${main_dir}/csstn_${sample}_${start_time}.sh"
					else
						qsub -sync y "${main_dir}/csstn_${sample}_${start_time}.sh"
					fi
				fi
				mv "${shareScript}/csstn_${sample}.out" ${main_dir}
				mv "${shareScript}/csstn_${sample}.err" ${main_dir}
			else
				echo "${project}/${sample} already has the newest ResGANNOT (${resGANNOT_srst2_filename})"
				echo -e "$(date)" > "${main_dir}/complete/${sample}_csstarn_complete.txt"
			fi
			if [[ -d "${processed}/${project}/${sample}/c-sstar_plasmid" ]]; then
				if [[ ! -f "${processed}/${project}/${sample}/c-sstar_plasmid/${sample}.${resGANNOT_srst2_filename}.gapped_40_sstar_summary.txt" ]]; then
					echo "Index below max submissions, submitting plasmid"
					echo -e "#!/bin/bash -l\n" > "${main_dir}/csstp_${sample}_${start_time}.sh"
					echo -e "#$ -o csstp_${sample}.out" >> "${main_dir}/csstp_${sample}_${start_time}.sh"
					echo -e "#$ -e csstp_${sample}.err" >> "${main_dir}/csstp_${sample}_${start_time}.sh"
					echo -e "#$ -N csstp_${sample}"   >> "${main_dir}/csstp_${sample}_${start_time}.sh"
					echo -e "#$ -cwd"  >> "${main_dir}/csstp_${sample}_${start_time}.sh"
					echo -e "#$ -q short.q\n"  >> "${main_dir}/csstp_${sample}_${start_time}.sh"
					echo -e "module load Python/3.6.1\n" >> "${main_dir}/csstp_${sample}_${start_time}.sh"
					# Defaulting to gapped/98, change if you want to include user preferences
					echo -e "cd ${shareScript}" >> "${main_dir}/csstp_${sample}_${start_time}.sh"
					echo -e "\"${shareScript}/run_c-sstar_on_single.sh\" \"${sample}\" g o \"${project}\" \"--plasmid\"" >> "${main_dir}/csstp_${sample}_${start_time}.sh"
					echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_csstarp_complete.txt\"" >> "${main_dir}/csstp_${sample}_${start_time}.sh"
					cd "${main_dir}"
					if [[ "${counter}" -lt "${last_index}" ]]; then
						qsub "${main_dir}/csstp_${sample}_${start_time}.sh"
					else
						qsub -sync y "${main_dir}/csstp_${sample}_${start_time}.sh"
					fi
					mv "${shareScript}/csstp_${sample}.out" ${main_dir}
					mv "${shareScript}/csstp_${sample}.err" ${main_dir}
				else
					echo "${project}/${sample} plasmid already has the newest ResGANNOT (${resGANNOT_srst2_filename})"
					echo -e "$(date)" > "${main_dir}/complete/${sample}_csstarp_complete.txt"
				fi
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
				if [[ -f "${main_dir}/complete/${waiting_sample}_csstarn_complete.txt" ]] || [[ ! -s "${processed}/${project}/${waiting_sample}/Assembly/${waiting_sample}_scaffolds_trimmed.fasta" ]]; then
					if [[ ! -f "${processed}/${project}/${sample}/c-sstar/${sample}.${resGANNOT_srst2_filename}.gapped_98_sstar_summary.txt" ]]; then
						echo  "Index is below max submissions, submitting"
						echo -e "#!/bin/bash -l\n" > "${main_dir}/csstn_${sample}_${start_time}.sh"
						echo -e "#$ -o csstn_${sample}.out" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
						echo -e "#$ -e csstn_${sample}.err" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
						echo -e "#$ -N csstn_${sample}"   >> "${main_dir}/csstn_${sample}_${start_time}.sh"
						echo -e "#$ -cwd"  >> "${main_dir}/csstn_${sample}_${start_time}.sh"
						echo -e "#$ -q short.q\n"  >> "${main_dir}/csstn_${sample}_${start_time}.sh"
						echo -e "module load Python/3.6.1\n" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
						# Defaulting to gapped/98, change if you want to include user preferences
						echo -e "cd ${shareScript}" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
						echo -e "\"${shareScript}/run_c-sstar_on_single.sh\" \"${sample}\" g h \"${project}\"" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
						echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_csstarn_complete.txt\"" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
						cd ${main_dir}
						if [[ "${counter}" -lt "${last_index}" ]]; then
							qsub "${main_dir}/csstn_${sample}_${start_time}.sh"
						else
							if [[ -d "${processed}/${project}/${sample}/c-sstar_plasmid" ]]; then
								qsub "${main_dir}/csstn_${sample}_${start_time}.sh"
							else
								qsub -sync y "${main_dir}/csstn_${sample}_${start_time}.sh"
							fi
						fi
						mv "${shareScript}/csstn_${sample}.out" ${main_dir}
						mv "${shareScript}/csstn_${sample}.err" ${main_dir}
					else
						echo "${project}/${sample} already has the newest ResGANNOT (${resGANNOT_srst2_filename})"
						echo -e "$(date)" > "${main_dir}/complete/${sample}_csstarn_complete.txt"
					fi
					if [[ -d "${processed}/${project}/${sample}/c-sstar_plasmid" ]]; then
						if [[ ! -f "${processed}/${project}/${sample}/c-sstar_plasmid/${sample}.${resGANNOT_srst2_filename}.gapped_40_sstar_summary.txt" ]]; then
							echo -e "#!/bin/bash -l\n" > "${main_dir}/csstp_${sample}_${start_time}.sh"
							echo -e "#$ -o csstp_${sample}.out" >> "${main_dir}/csstp_${sample}_${start_time}.sh"
							echo -e "#$ -e csstp_${sample}.err" >> "${main_dir}/csstp_${sample}_${start_time}.sh"
							echo -e "#$ -N csstp_${sample}"   >> "${main_dir}/csstp_${sample}_${start_time}.sh"
							echo -e "#$ -cwd"  >> "${main_dir}/csstp_${sample}_${start_time}.sh"
							echo -e "#$ -q short.q\n"  >> "${main_dir}/csstp_${sample}_${start_time}.sh"
							echo -e "module load Python/3.6.1\n" >> "${main_dir}/csstp_${sample}_${start_time}.sh"
							# Defaulting to gapped/98, change if you want to include user preferences
							echo -e "cd ${shareScript}" >> "${main_dir}/csstp_${sample}_${start_time}.sh"
							echo -e "\"${shareScript}/run_c-sstar_on_single.sh\" \"${sample}\" g o \"${project}\" \"--plasmid\"" >> "${main_dir}/csstp_${sample}_${start_time}.sh"
							echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_csstarp_complete.txt\"" >> "${main_dir}/csstp_${sample}_${start_time}.sh"
							cd ${main_dir}
							if [[ "${counter}" -lt "${last_index}" ]]; then
								qsub "${main_dir}/csstp_${sample}_${start_time}.sh"
							else
								qsub -sync y "${main_dir}/csstp_${sample}_${start_time}.sh"
							fi
							mv "${shareScript}/csstp_${sample}.out" ${main_dir}
							mv "${shareScript}/csstp_${sample}.err" ${main_dir}
						else
							echo "${project}/${sample} plasmid already has the newest ResGANNOT (${resGANNOT_srst2_filename})"
							echo -e "$(date)" > "${main_dir}/complete/${sample}_csstarp_complete.txt"
						fi
					fi
					break
				else
					timer=$(( timer + 5 ))
					echo "${main_dir}/complete/${waiting_sample}_csstarn_complete.txt not ready, sleeping for 5 seconds, so far slept for ${timer}"
					sleep 5
				fi
			done
		fi
	fi
	counter=$(( counter + 1 ))
done

# Check for completion of all samples
timer=0

for item in "${arr[@]}"; do
	waiting_sample=$(echo "${item}" | cut -d'/' -f2)
	if [[ -f "${main_dir}/complete/${waiting_sample}_csstarn_complete.txt" ]] || [[ ! -s "${processed}/${project}/${waiting_sample}/Assembly/${waiting_sample}_scaffolds_trimmed.fasta" ]]; then
		echo "${item} is complete"
	else
		while :
		do
				if [[ ${timer} -gt 3600 ]]; then
					echo "Timer exceeded limit of 3600 seconds = 60 minutes"
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

exit 0
