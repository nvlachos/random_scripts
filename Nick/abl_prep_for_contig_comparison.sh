#!/bin/sh -l

#$ -o ablmq_contig.out
#$ -e ablmq_contig.err
#$ -N ablmq_contig
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
# Usage ./abl_mass_qsub_contig_fix.sh path_to_list max_concurrent_submission output_directory_for_scripts
#

# Number regex to test max concurrent submission parametr
number='^[0-9]+$'

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to abl_mass_sub_contig.sh, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./abl_mass_qsub_contig.sh path_to_list_file(single sample ID per line, e.g. B8VHY/1700128 (it must include project id also)) max_submissions output_directory_for_scripts"
	echo "Output location varies depending on which tasks are performed but will be found somewhere under ${processed}"
	exit 0
elif [[ ! -f "${1}" ]]; then
	echo "${1} (list) does not exist...exiting"
	exit 1
elif ! [[ ${2} =~ $number ]] || [[ -z "${2}" ]]; then
	echo "${2} is not a number or is empty. Please input max number of concurrent qsub submissions...exiting"
	exit 2
elif [[ -z "${3}" ]]; then
	echo "No script output directory given...exiting"
	exit 3
fi

# Make array to hold all sample names to be processed
arr=()
while IFS= read -r line || [[ "$line" ]];  do
  arr+=("$line")
done < ${1}

arr_size="${#arr[@]}"
last_index=$(( arr_size -1 ))
echo "-${arr_size}:${arr[@]}-"

max_subs=${2}

# Loop through and act on each sample name in the passed/provided list
counter=0
main_dir="${3}/contig_subs"
if [[ ! -d "${3}/contig_subs" ]]; then
	mkdir "${3}/contig_subs"
	mkdir "${3}/contig_subs/complete"
elif [[ ! -d  "${3}/contig_subs/complete" ]]; then
	mkdir "${3}/contig_subs/complete"
fi

# Loop through all samples in the array
while [ ${counter} -lt ${arr_size} ] ; do
	sample=$(echo "${arr[${counter}]}" | cut -d'/' -f2)
	project=$(echo "${arr[${counter}]}" | cut -d'/' -f1)
	# echo ${counter}
	if [ ${counter} -lt ${max_subs} ]; then
		if [[ -s "${processed}/${project}/${sample}/Assembly/contigs.fasta" ]]; then
			echo  "Index is below max submissions, submitting"
			echo -e "#!/bin/bash -l\n" > "${main_dir}/contig_${sample}_${start_time}.sh"
			echo -e "#$ -o contig_${sample}.out" >> "${main_dir}/contig_${sample}_${start_time}.sh"
			echo -e "#$ -e contig_${sample}.err" >> "${main_dir}/contig_${sample}_${start_time}.sh"
			echo -e "#$ -N contig_${sample}"   >> "${main_dir}/contig_${sample}_${start_time}.sh"
			echo -e "#$ -cwd"  >> "${main_dir}/contig_${sample}_${start_time}.sh"
			echo -e "#$ -q short.q\n"  >> "${main_dir}/contig_${sample}_${start_time}.sh"
			echo -e "cd ${shareScript}" >> "${main_dir}/contig_${sample}_${start_time}.sh"
			echo -e "python \"${shareScript}/removeShortContigs.py\" \"${processed}/${project}/${sample}/Assembly/contigs.fasta\" 500" >> "${main_dir}/contig_${sample}_${start_time}.sh"
			echo -e "mkdir ${processed}/${project}/${sample}/Contig_check" >> "${main_dir}/contig_${sample}_${start_time}.sh"
			echo -e "mv \"${processed}/${project}/${sample}/Assembly/contigs.fasta.TRIMMED.fasta" "${processed}/${project}/${sample}/Contig_check/${sample}_contigs_trimmed_original.fasta" >> "${main_dir}/contig_${sample}_${start_time}.sh"
			echo -e "python3 \"${shareScript}/fasta_headers.py\" \"${processed}/${project}/${sample}/Contig_check/${sample}_contigs_trimmed_original.fasta\" \"${processed}/${project}/${sample}/Contig_check/${sample}_contigs_trimmed.fasta\"" >> "${main_dir}/contig_${sample}_${start_time}.sh"
			echo -e "${shareScript}/run_c-sstar_on_contigs.sh \"${sample}\" g h \"${project}\"" >> "${main_dir}/contig_${sample}_${start_time}.sh"
			echo -e "${shareScript}/run_plasmidFinder_on_contigs.sh \"${sample}\" \"${project}\" plasmid" >> "${main_dir}/contig_${sample}_${start_time}.sh"
			echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_contig_complete.txt\"" >> "${main_dir}/contig_${sample}_${start_time}.sh"
			cd "${main_dir}"
			qsub "${main_dir}/contig_${sample}_${start_time}.sh"
		# No contig assembly is present to work on
		else
			echo "${project}/${sample} has no contig Assembly to fix"
			echo "$(date)" > "${main_dir}/complete/${sample}_contig_complete.txt"
		fi
		# Counter is above max subs, must wait until a slot opens up
	else
		waiting_for_index=$(( counter - max_subs ))
		waiting_sample=$(echo "${arr[${waiting_for_index}]}" | cut -d'/' -f2)
		timer=0
		echo "Index is above max submissions, waiting for index ${waiting_for_index}:${waiting_sample} to complete"
		# Loop to wait until slot opens up
		while :
		do
			# Check if max allowed time has elapsed
			if [[ ${timer} -gt 1800 ]]; then
				echo "Timer exceeded limit of 1800 seconds 30 minutes"
				break
			fi
			# Check if "waiting" sample has completed
			if [ -f "${main_dir}/complete/${waiting_sample}_contig_complete.txt" ]; then
				# Check if an acceptable assembly is available to use
				if [[ -s "${processed}/${project}/${sample}/Assembly/contigs.fasta" ]]; then
					# Check if contig fixer has already been run on sample
					echo  "Index is below max submissions, submitting"
					echo -e "#!/bin/bash -l\n" > "${main_dir}/contig_${sample}_${start_time}.sh"
					echo -e "#$ -o contig_${sample}.out" >> "${main_dir}/contig_${sample}_${start_time}.sh"
					echo -e "#$ -e contig_${sample}.err" >> "${main_dir}/contig_${sample}_${start_time}.sh"
					echo -e "#$ -N contig_${sample}"   >> "${main_dir}/contig_${sample}_${start_time}.sh"
					echo -e "#$ -cwd"  >> "${main_dir}/contig_${sample}_${start_time}.sh"
					echo -e "#$ -q short.q\n"  >> "${main_dir}/contig_${sample}_${start_time}.sh"
					echo -e "cd ${shareScript}" >> "${main_dir}/contig_${sample}_${start_time}.sh"
					echo -e "python \"${shareScript}/removeShortContigs.py\" \"${processed}/${project}/${sample}/Assembly/contigs.fasta\" 500" >> "${main_dir}/contig_${sample}_${start_time}.sh"
					echo -e "mkdir ${processed}/${project}/${sample}/Contig_check" >> "${main_dir}/contig_${sample}_${start_time}.sh"
					echo -e "mv \"${processed}/${project}/${sample}/Assembly/contigs.fasta.TRIMMED.fasta" "${processed}/${project}/${sample}/Contig_check/${sample}_contigs_trimmed_original.fasta" >> "${main_dir}/contig_${sample}_${start_time}.sh"
					echo -e "python3 \"${shareScript}/fasta_headers.py\" \"${processed}/${project}/${sample}/Contig_check/${sample}_contigs_trimmed_original.fasta\" \"${processed}/${project}/${sample}/Contig_check/${sample}_contigs_trimmed.fasta\"" >> "${main_dir}/contig_${sample}_${start_time}.sh"
					echo -e "${shareScript}/run_c-sstar_on_contigs.sh \"${sample}\" g h \"${project}\"" >> "${main_dir}/contig_${sample}_${start_time}.sh"
					echo -e "${shareScript}/run_plasmidFinder_on_contigs.sh \"${sample}\" \"${project}\" plasmid" >> "${main_dir}/contig_${sample}_${start_time}.sh"
					echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_contig_complete.txt\"" >> "${main_dir}/contig_${sample}_${start_time}.sh"
					cd "${main_dir}"
					qsub "${main_dir}/contig_${sample}_${start_time}.sh"
					break
				# No assembly file to fix
				else
					echo "${project}/${sample} has no Assembly to fix"
					echo "$(date)" > "${main_dir}/complete/${sample}_contig_complete.txt"
					break
				fi
				# Wait 5 seconds and then check if "waiting" sample completed yet
			else
				timer=$(( timer + 5 ))
				echo "sleeping for 5 seconds, so far slept for ${timer}"
				sleep 5
			fi
		done
	fi
	counter=$(( counter + 1 ))
done

# Check for completion of all samples
timer=0
for item in "${arr[@]}"; do
	waiting_sample=$(echo "${item}" | cut -d'/' -f2)
	if [[ -f "${main_dir}/complete/${waiting_sample}_contig_complete.txt" ]] || [[ ! -s "${processed}/${project}/${waiting_sample}/Assembly/${waiting_sample}_contigs_trimmed.fasta" ]]; then
		echo "${item} is complete"
		if [[ -f "${shareScript}/contig_${sample}.out" ]]; then
			mv "${shareScript}/contig_${sample}.out" ${main_dir}
		fi
		if [[ -f "${shareScript}/contig_${sample}.err" ]]; then
 			mv "${shareScript}/contig_${sample}.err" ${main_dir}
		fi
	else
		while :
		do
				if [[ ${timer} -gt 3600 ]]; then
					echo "Timer exceeded limit of 3600 seconds = 60 minutes"
					exit 1
				fi
				if [[ -f "${main_dir}/complete/${waiting_sample}_contig_complete.txt" ]]; then
					echo "${item} is complete"
					if [[ -f "${shareScript}/contig_${sample}.out" ]]; then
						mv "${shareScript}/contig_${sample}.out" ${main_dir}
					fi
					if [[ -f "${shareScript}/contig_${sample}.err" ]]; then
			 			mv "${shareScript}/contig_${sample}.err" ${main_dir}
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
printf "%s %s" "abl_mass_qsub_contig.sh has completed" "${global_end_time}" | mail -s "abl_mass_qsub_contig.sh complete" nvx4@cdc.gov

#Script exited gracefully (unless something else inside failed)
exit 0
