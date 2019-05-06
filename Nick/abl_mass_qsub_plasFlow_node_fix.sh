#!/bin/sh -l

#$ -o ablmq_plnode.out
#$ -e ablmq_plnode.err
#$ -N ablmq_plnode
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
# Usage ./abl_mass_qsub_plasFlow_node_fix.sh path_to_list max_concurrent_submission output_directory_for_scripts clobberness[keep|clobber]
#

# Number regex to test max concurrent submission parametr
number='^[0-9]+$'

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to abl_mass_sub_node.sh, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./abl_mass_qsub_plasFlow_node.sh path_to_list_file(single sample ID per line, e.g. B8VHY/1700128 (it must include project id also)) max_submissions output_directory_for_scripts"
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
main_dir="${3}/node_subs"
if [[ ! -d "${3}/node_subs" ]]; then
	mkdir "${3}/node_subs"
	mkdir "${3}/node_subs/complete"
elif [[ ! -d  "${3}/node_subs/complete" ]]; then
	mkdir "${3}/node_subs/complete"
fi

# Loop through all samples in the array
while [ ${counter} -lt ${arr_size} ] ; do
	sample=$(echo "${arr[${counter}]}" | cut -d'/' -f2)
	project=$(echo "${arr[${counter}]}" | cut -d'/' -f1)
	if [[ "${clobberness}" == "clobber" ]]; then
		if [[ -s "${processed}/${project}/${sample}/plasFlow/Unicycler_assemblies/${sample}_uni_assembly/${sample}_plasmid_assembly.fasta" ]]; then
			header=$(head -n1 "${processed}/${project}/${sample}/plasFlow/Unicycler_assemblies/${sample}_uni_assembly/${sample}_plasmid_assembly.fasta" | cut -d'_' -f1)
			if [[ "${header}" = ">${sample}" ]]; then
				:
			else
				"${processed}/${project}/${sample}/plasFlow/Unicycler_assemblies/${sample}_uni_assembly/${sample}_plasmid_assembly.fasta"
				mv "${processed}/${project}/${sample}/plasFlow/Unicycler_assemblies/${sample}_uni_assembly/${sample}_plasmid_assembly_original.fasta" "${processed}/${project}/${sample}/plasFlow/Unicycler_assemblies/${sample}_uni_assembly/${sample}_plasmid_assembly.fasta"
			fi
		fi
	fi
	# echo ${counter}
	if [ ${counter} -lt ${max_subs} ]; then
		if [[ -s "${processed}/${project}/${sample}/plasFlow/Unicycler_assemblies/${sample}_uni_assembly/${sample}_plasmid_assembly.fasta" ]]; then
			header=$(head -n1 "${processed}/${project}/${sample}/plasFlow/Unicycler_assemblies/${sample}_uni_assembly/${sample}_plasmid_assembly.fasta" | cut -d'_' -f1)
			# Check if Node fixer has been run already, skip if so
			if [[ "${header}" != ">${sample}" ]]; then
				echo  "Index is below max submissions, submitting"
				echo -e "#!/bin/bash -l\n" > "${main_dir}/plnode_${sample}_${start_time}.sh"
				echo -e "#$ -o plnode_${sample}.out" >> "${main_dir}/plnode_${sample}_${start_time}.sh"
				echo -e "#$ -e plnode_${sample}.err" >> "${main_dir}/plnode_${sample}_${start_time}.sh"
				echo -e "#$ -N plnode_${sample}"   >> "${main_dir}/plnode_${sample}_${start_time}.sh"
				echo -e "#$ -cwd"  >> "${main_dir}/plnode_${sample}_${start_time}.sh"
				echo -e "#$ -q short.q\n"  >> "${main_dir}/plnode_${sample}_${start_time}.sh"
				echo -e "cd ${shareScript}" >> "${main_dir}/plnode_${sample}_${start_time}.sh"
				echo -e "mv \"${processed}/${project}/${sample}/plasFlow/Unicycler_assemblies/${sample}_uni_assembly/${sample}_plasmid_assembly.fasta\" \"${processed}/${project}/${sample}/plasFlow/Unicycler_assemblies/${sample}_uni_assembly/${sample}_plasmid_assembly_original.fasta\"" >> "${main_dir}/plnode_${sample}_${start_time}.sh"
				echo -e "python3 \"${shareScript}/fasta_headers_plasFlow.py\" \"${processed}/${project}/${sample}/plasFlow/Unicycler_assemblies/${sample}_uni_assembly/${sample}_plasmid_assembly_original.fasta\" \"${processed}/${project}/${sample}/plasFlow/Unicycler_assemblies/${sample}_uni_assembly/${sample}_plasmid_assembly.fasta\"" >> "${main_dir}/plnode_${sample}_${start_time}.sh"
				echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_plnode_complete.txt\"" >> "${main_dir}/plnode_${sample}_${start_time}.sh"
				cd "${main_dir}"
				#if [[ "${counter}" -lt "${last_index}" ]]; then
					qsub "${main_dir}/plnode_${sample}_${start_time}.sh"
				#else
				#	qsub -sync y "${main_dir}/plnode_${sample}_${start_time}.sh"
				#fi
			# Nodes already have been fixed
			else
				echo "${project}/${sample} already had its nodes removed"
				echo "$(date)" > "${main_dir}/complete/${sample}_plnode_complete.txt"
			fi
		# No assembly is present to work on
		else
			echo "${project}/${sample} has no Assembly to fix"
			echo "$(date)" > "${main_dir}/complete/${sample}_plnode_complete.txt"
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
			if [ -f "${main_dir}/complete/${waiting_sample}_plnode_complete.txt" ]; then
				# Check if an acceptable assembly is available to use
				if [[ -s "${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed.fasta" ]]; then
					# Check if node fixer has already been run on sample
					header=$(head -n1 "${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed.fasta" | cut -d'_' -f1)
					if [[ "${header}" = ">NODE" ]]; then
						echo  "Index is below max submissions, submitting"
						echo -e "#!/bin/bash -l\n" > "${main_dir}/plnode_${sample}_${start_time}.sh"
						echo -e "#$ -o plnode_${sample}.out" >> "${main_dir}/plnode_${sample}_${start_time}.sh"
						echo -e "#$ -e plnode_${sample}.err" >> "${main_dir}/plnode_${sample}_${start_time}.sh"
						echo -e "#$ -N plnode_${sample}"   >> "${main_dir}/plnode_${sample}_${start_time}.sh"
						echo -e "#$ -cwd"  >> "${main_dir}/plnode_${sample}_${start_time}.sh"
						echo -e "#$ -q short.q\n"  >> "${main_dir}/plnode_${sample}_${start_time}.sh"
						echo -e "cd ${shareScript}" >> "${main_dir}/plnode_${sample}_${start_time}.sh"
						echo -e "mv \"${processed}/${project}/${sample}/plasFlow/Unicycler_assemblies/${sample}_uni_assembly/${sample}_plasmid_assembly.fasta\" \"${processed}/${project}/${sample}/plasFlow/Unicycler_assemblies/${sample}_uni_assembly/${sample}_plasmid_assembly_original.fasta\"" >> "${main_dir}/plnode_${sample}_${start_time}.sh"
						echo -e "python3 \"${shareScript}/fasta_headers_plasFlow.py\" \"${processed}/${project}/${sample}/plasFlow/Unicycler_assemblies/${sample}_uni_assembly/${sample}_plasmid_assembly_original.fasta\" \"${processed}/${project}/${sample}/plasFlow/Unicycler_assemblies/${sample}_uni_assembly/${sample}_plasmid_assembly.fasta\"" >> "${main_dir}/plnode_${sample}_${start_time}.sh"
						echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_plnode_complete.txt\"" >> "${main_dir}/plnode_${sample}_${start_time}.sh"
						cd "${main_dir}"
						#if [[ "${counter}" -lt "${last_index}" ]]; then
							qsub "${main_dir}/plnode_${sample}_${start_time}.sh"
						#else
						#	qsub -sync y "${main_dir}/plnode_${sample}_${start_time}.sh"
						#fi
						break
					# Node fixer already run on sample
					else
						echo "${project}/${sample} already had its nodes removed"
						echo "$(date)" > "${main_dir}/complete/${sample}_plnode_complete.txt"
						break
					fi
				# No assembly file to fix
				else
					echo "${project}/${sample} has no Assembly to fix"
					echo "$(date)" > "${main_dir}/complete/${sample}_plnode_complete.txt"
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
	if [[ -f "${main_dir}/complete/${waiting_sample}_plnode_complete.txt" ]] || [[ ! -s "${processed}/${project}/${sample}/plasFlow/Unicycler_assemblies/${waitin_sample}_uni_assembly/${waiting_sample}_plasmid_assembly.fasta" ]]; then
		echo "${item} is complete"
		if [[ -f "${shareScript}/plnode_${sample}.out" ]]; then
			mv "${shareScript}/plnode_${sample}.out" ${main_dir}
		fi
		if [[ -f "${shareScript}/plnode_${sample}.err" ]]; then
 			mv "${shareScript}/plnode_${sample}.err" ${main_dir}
		fi
	else
		while :
		do
				if [[ ${timer} -gt 3600 ]]; then
					echo "Timer exceeded limit of 3600 seconds = 60 minutes"
					exit 1
				fi
				if [[ -f "${main_dir}/complete/${waiting_sample}_plnode_complete.txt" ]]; then
					echo "${item} is complete"
					if [[ -f "${shareScript}/plnode_${sample}.out" ]]; then
						mv "${shareScript}/plnode_${sample}.out" ${main_dir}
					fi
					if [[ -f "${shareScript}/plnode_${sample}.err" ]]; then
			 			mv "${shareScript}/plnode_${sample}.err" ${main_dir}
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
printf "%s %s" "abl_mass_qsub_node.sh has completed" "${global_end_time}" | mail -s "abl_mass_qsub_node.sh complete" nvx4@cdc.gov

#Script exited gracefully (unless something else inside failed)
exit 0
