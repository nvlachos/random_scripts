#!/bin/sh -l

#$ -o ablmq-16s.out
#$ -e ablmq-16s.err
#$ -N ablmq-16s
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
. ./config.sh
#Import the module file that loads all necessary mods
. "${mod_changers}/pipeline_mods"

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
  arr+=("$line")
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
main_dir="${share}/mass_subs/blast16s_subs"
if [[ ! -d "${share}/mass_subs/blast16s_subs" ]]; then
	mkdir "${share}/mass_subs/blast16s_subs"
	mkdir "${share}/mass_subs/blast16s_subs/complete"
elif [[ ! -d "${share}/mass_subs/blast16s_subs/complete" ]]; then
	mkdir "${share}/mass_subs/blast16s_subs/complete"
fi

start_time=$(DATE)

while [ ${counter} -lt ${arr_size} ] ; do
	sample=$(echo "${arr[${counter}]}" | cut -d'/' -f2)
	project=$(echo "${arr[${counter}]}" | cut -d'/' -f1)
	#rm -r "${processed}/${project}/${sample}/16s/${sample}.16s"
	#echo ${counter}"-${processed}/${project}/${sample}/kraken/postAssembly/${sample}_kraken_summary_assembled_BP_data.txt"
	if [[ -s "${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed.fasta" ]]; then
		if [[ ${counter} -lt ${max_subs} ]]; then
			if [[ ! -f "${processed}/${project}/${sample}/16s/${sample}_16s_blast_id.txt" ]]; then
				echo  "Index is below max submissions, submitting"
				echo -e "#!/bin/bash -l\n" > "${main_dir}/blast16s_${sample}_${start_time}.sh"
				echo -e "#$ -o blast16s_${sample}.out" >> "${main_dir}/blast16s_${sample}_${start_time}.sh"
				echo -e "#$ -e blast16s_${sample}.err" >> "${main_dir}/blast16s_${sample}_${start_time}.sh"
				echo -e "#$ -N blast16s_${sample}" >> "${main_dir}/blast16s_${sample}_${start_time}.sh"
				echo -e "#$ -cwd"  >> "${main_dir}/blast16s_${sample}_${start_time}.sh"
				echo -e "#$ -q short.q\n" >> "${main_dir}/blast16s_${sample}_${start_time}.sh"
				# Defaulting to gapped/98, change if you want to include user preferences
				echo -e "\"${shareScript}/16s_blast.sh\" \"${sample}\" \"${project}\"" >> "${main_dir}/blast16s_${sample}_${start_time}.sh"
				echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_blast16s_complete.txt\"" >> "${main_dir}/blast16s_${sample}_${start_time}.sh"
				qsub "${main_dir}/blast16s_${sample}_${start_time}.sh"
			else
				echo "${project}/${sample} already has 16s summary"
				echo "$(date)" > "${main_dir}/complete/${sample}_blast16s_complete.txt"
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
				if [[ -f "${main_dir}/complete/${waiting_sample}_blast16s_complete.txt" ]]; then
					if [[ ! -f "${processed}/${project}/${sample}/16s/${sample}_16s_blast_id.txt" ]]; then
						echo  "Index is below max submissions, submitting"
						echo -e "#!/bin/bash -l\n" > "${main_dir}/blast16s_${sample}_${start_time}.sh"
						echo -e "#$ -o blast16s_${sample}.out" >> "${main_dir}/blast16s_${sample}_${start_time}.sh"
						echo -e "#$ -e blast16s_${sample}.err" >> "${main_dir}/blast16s_${sample}_${start_time}.sh"
						echo -e "#$ -N blast16s_${sample}" >> "${main_dir}/blast16s_${sample}_${start_time}.sh"
						echo -e "#$ -cwd"  >> "${main_dir}/blast16s_${sample}_${start_time}.sh"
						echo -e "#$ -q short.q\n" >> "${main_dir}/blast16s_${sample}_${start_time}.sh"
						echo -e "\"${shareScript}/16s_blast.sh\" \"${sample}\" \"${project}\"" >> "${main_dir}/blast16s_${sample}_${start_time}.sh"
						echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_blast16s_complete.txt\"" >> "${main_dir}/blast16s_${sample}_${start_time}.sh"
						qsub "${main_dir}/blast16s_${sample}_${start_time}.sh"
					else
						echo "${project}/${sample} already has 16s summary"
						echo "$(date)" > "${main_dir}/complete/${sample}_blast16s_complete.txt"
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
