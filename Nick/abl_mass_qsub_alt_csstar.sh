#!/bin/sh -l

#$ -o ablmq-csn.out
#$ -e ablmq-csn.err
#$ -N ablmq-csn
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
	echo "Usage is ./abl_mass_qsub_csstar_alt.sh path_to_list_file(single sample ID per line, e.g. B8VHY/1700128 (it must include project id also)) max_submissions path_to_alt_database"
	echo "Output location varies depending on which tasks are performed but will be found somewhere under ${processed}"
	exit 0
fi

arr=()
while IFS= read -r line || [[ "$line" ]];  do
  arr+=("$line")
done < ${1}

arr_size="${#arr[@]}"
echo "-${arr_size}:${arr[@]}-"


# Loop through and act on each sample name in the passed/provided list
counter=0
max_subs=${2}
main_dir="${share}/mass_subs/csstar_alt_subs"
if [[ ! -d "${share}/mass_subs/csstar_alt_subs" ]]; then
	mkdir "${share}/mass_subs/csstar_alt_subs"
	mkdir "${share}/mass_subs/csstar_alt_subs/complete"
elif [[ ! -d  "${share}/mass_subs/csstar_alt_subs/complete" ]]; then
	mkdir "${share}/mass_subs/csstar_alt_subs/complete"
fi


alt_database_path=$(basename -- "${5}")
alt_database=$(echo ${alt_database_path##*/} | cut -d'.' -f1)
alt_database=${alt_database//_srst2/}

exit

start_time=$(date "+%F-%T")

while [ ${counter} -lt ${arr_size} ] ; do
	sample=$(echo "${arr[${counter}]}" | cut -d'/' -f2)
	project=$(echo "${arr[${counter}]}" | cut -d'/' -f1)
	rm -r "${processed}/${project}/${sample}/c-sstar/alt_database=_gapped/"
	rm -r "${processed}/${project}/${sample}/c-sstar/${sample}.alt_database=.gapped_98_sstar_summary.txt"
	rm -r "${processed}/${project}/${sample}/c-sstar/${alt_database}_gapped/"
	rm -r "${processed}/${project}/${sample}/c-sstar/${sample}.${alt_database}.gapped_98_sstar_summary.txt"
	echo ${counter}
	if [[ -s "${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed.fasta" ]]; then
		if [[ ${counter} -lt ${max_subs} ]]; then
			echo  "Index is below max submissions, submitting"

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
				echo -e "\"${shareScript}/run_c-sstar_on_single_alternate_DB.sh\" \"${sample}\" g h \"${project}\" \"${alt_database_path}\"" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
				echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_csstarn_complete.txt\"" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
				qsub "${main_dir}/csstn_${sample}_${start_time}.sh"
			else
				echo "${project}/${sample} already has 0608"
				echo "$(date)" > "${main_dir}/complete/${sample}_csstarn_complete.txt"
			fi

			if [[ -d "${processed}/${project}/${sample}/c-sstar_plasmid" ]]; then
				if [[ ! -f "${processed}/${project}/${sample}/c-sstar_plasmid/${sample}.${alt_database}.gapped_40_sstar_summary.txt" ]]; then
					echo -e "#!/bin/bash -l\n" > "${main_dir}/csstp_${sample}_${start_time}.sh"
					echo -e "#$ -o csstp_${sample}.out" >> "${main_dir}/csstp_${sample}_${start_time}.sh"
					echo -e "#$ -e csstp_${sample}.err" >> "${main_dir}/csstp_${sample}_${start_time}.sh"
					echo -e "#$ -N csstp_${sample}"   >> "${main_dir}/csstp_${sample}_${start_time}.sh"
					echo -e "#$ -cwd"  >> "${main_dir}/csstp_${sample}_${start_time}.sh"
					echo -e "#$ -q short.q\n"  >> "${main_dir}/csstp_${sample}_${start_time}.sh"
					echo -e "module load Python/3.6.1\n" >> "${main_dir}/csstp_${sample}_${start_time}.sh"
					# Defaulting to gapped/98, change if you want to include user preferences
					echo -e "\"${shareScript}/run_c-sstar_on_single_alternate_DB.sh\" \"${sample}\" g o \"${project}\" \"${alt_database_path}\" \"--plasmid\"" >> "${main_dir}/csstp_${sample}_${start_time}.sh"
					echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_csstarp_complete.txt\"" >> "${main_dir}/csstp_${sample}_${start_time}.sh"
					qsub "${main_dir}/csstp_${sample}_${start_time}.sh"
					if [[ -f "${processed}/${project}/${sample}/c-sstar_plasmid/${sample}.${alt_database}.gapped_98_sstar_summary.txt" ]]; then
						rm "${processed}/${project}/${sample}/c-sstar_plasmid/${sample}.${alt_database}.gapped_98_sstar_summary.txt"
					fi
				else
					echo "${project}/${sample} already has 0608 PLASMID"
					echo "$(date)" > "${main_dir}/complete/${sample}_csstarp_complete.txt"
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
						echo -e "\"${shareScript}/run_c-sstar_on_single_alternate_DB.sh\" \"${sample}\" g h \"${project}\" \"${alt_database_path}\"" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
						echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_csstarn_complete.txt\"" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
						qsub "${main_dir}/csstn_${sample}_${start_time}.sh"
					else
						echo "${project}/${sample} already has 0608"
						echo "$(date)" > "${main_dir}/complete/${sample}_csstarn_complete.txt"
					fi

					if [[ -d "${processed}/${project}/${sample}/c-sstar_plasmid" ]]; then
						if [[ ! -f "${processed}/${project}/${sample}/c-sstar_plasmid/${sample}.${alt_database}.gapped_40_sstar_summary.txt" ]]; then
							echo -e "#!/bin/bash -l\n" > "${main_dir}/csstp_${sample}_${start_time}.sh"
							echo -e "#$ -o csstp_${sample}.out" >> "${main_dir}/csstp_${sample}_${start_time}.sh"
							echo -e "#$ -e csstp_${sample}.err" >> "${main_dir}/csstp_${sample}_${start_time}.sh"
							echo -e "#$ -N csstp_${sample}"   >> "${main_dir}/csstp_${sample}_${start_time}.sh"
							echo -e "#$ -cwd"  >> "${main_dir}/csstp_${sample}_${start_time}.sh"
							echo -e "#$ -q short.q\n"  >> "${main_dir}/csstp_${sample}_${start_time}.sh"
							echo -e "module load Python/3.6.1\n" >> "${main_dir}/csstp_${sample}_${start_time}.sh"
							# Defaulting to gapped/98, change if you want to include user preferences
							echo -e "\"${shareScript}/run_c-sstar_on_single_alternate_DB.sh\" \"${sample}\" g o \"${project}\" \"/scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/DBs/star/${alt_database}.fasta\" \"--plasmid\"" >> "${main_dir}/csstp_${sample}_${start_time}.sh"
							echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_csstarp_complete.txt\"" >> "${main_dir}/csstp_${sample}_${start_time}.sh"
							qsub "${main_dir}/csstp_${sample}_${start_time}.sh"
							if [[ -f "${processed}/${project}/${sample}/c-sstar_plasmid/${sample}.${alt_database}.gapped_98_sstar_summary.txt" ]]; then
								rm "${processed}/${project}/${sample}/c-sstar_plasmid/${sample}.${alt_database}.gapped_98_sstar_summary.txt"
							fi
						else
							echo "${project}/${sample} already has 0608 PLASMID"
							echo "$(date)" > "${main_dir}/complete/${sample}_csstarp_complete.txt"
						fi
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
		echo "${project}/${sample} does not have an assembly file???"
	fi
	counter=$(( counter + 1 ))
done

echo "All isolates completed"
global_end_time=$(date "+%m-%d-%Y @ %Hh_%Mm_%Ss")
#Script exited gracefully (unless something else inside failed)
#printf "%s %s" "Act_by_list.sh has completed ${2}" "${global_end_time}" | mail -s "act_by_list complete" nvx4@cdc.gov
exit 0
