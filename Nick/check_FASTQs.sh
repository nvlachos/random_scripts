#!/bin/sh -l

#$ -o act_by_list_barebones_2.out
#$ -e act_by_list_barebones_2.err
#$ -N ablb2
#$ -cwd
#$ -q all.q

#Import the config file with shortcuts and settings
. /scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/scripts/config.sh
#Import the module file that loads all necessary mods
. "${mod_changers}/pipeline_mods"

#
# Usage ./act_by_list.sh list_name(currently has to be placed in /scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR folder) Description of list function
#
# script changes depending on what needs to be run through the list
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to act_by_list.sh, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./act_by_list.sh path_to_list_file(single sample ID per line, e.g. B8VHY/1700128 (it must include project id also))"
	echo "Output location varies depending on which tasks are performed but will be found somewhere under ${processed}"
	exit 0
fi

# Loop through and act on each sample name in the passed/provided list
counter=0
while IFS= read -r var; do
	FR1="-"
	FR2="-"
	FQ1="-"
	FQ2="-"	
	TR1="-"
	TR2="-"
	sample_name=$(echo "${var}" | cut -d'/' -f2 | tr -d '[:space:]')
	project=$(echo "${var}" | cut -d'/' -f1 | tr -d '[:space:]')
	#echo ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R1_001.fastq
	if [[ ! -f ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R1_001.fastq ]]; then
		if [[ -f ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R1_001.fastq.gz ]]; then
			#rm -r ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R1_001.fastq
			FQ1="G"
		elif [[ -f ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R1.fastq ]]; then
			mv ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R1.fastq ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R1_001.fastq
			gzip ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R1_001.fastq
			rm -r ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R1.fastq
			FQ1="R"
		fi
	else
		if [[ -f ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R1_001.fastq.gz ]]; then
			rm -r ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R1_001.fastq 
			FQ1="G"
		else
			gzip ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R1_001.fastq
			rm -r ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R1_001.fastq
			FQ1="Z"
		fi
	fi
	#echo "${FQ1}"
	
	#echo ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R2_001.fastq
	if [[ ! -f ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R2_001.fastq ]]; then
		if [[ -f ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R2_001.fastq.gz ]]; then
			#rm -r ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R2_001.fastq
			FQ2="G"
		elif [[ -f ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R2.fastq ]]; then
			mv ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R2.fastq ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R2_001.fastq
			gzip ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R2_001.fastq
			rm -r ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R2.fastq
			FQ2="R"
		fi
	else
		if [[ -f ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R2_001.fastq.gz ]]; then
			rm -r ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R2_001.fastq 
			FQ2="G"
		else
			gzip ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R2_001.fastq
			rm -r ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R2_001.fastq
			FQ2="Z"
		fi
	fi
	#echo "${FQ2}"
	
	#echo ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R1_001.fastq.gz
	if [[ ! -f ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R1_001.fastq.gz ]]; then
		if [[ -f ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R1.fastq.gz ]]; then
			mv ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R1.fastq.gz ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R1_001.fastq.gz
			FR1="M"
		elif [[ -f ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R1.fastq ]]; then
			gzip ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R1.fastq
			mv ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R1.fastq.gz ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R1_001.fastq.gz
			FR1="Z"
		else
			FR1="X"
		fi
	else
		FR1="G"
	fi
	#echo "${FR1}"
	
	#echo ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R2_001.fastq.gz
	if [[ ! -f ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R2_001.fastq.gz ]]; then
		if [[ -f ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R2.fastq.gz ]]; then
			mv ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R2.fastq.gz ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R2_001.fastq.gz
			FR2="M"
		elif [[ -f ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R2.fastq ]]; then
			gzip ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R2.fastq
			mv ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R2.fastq.gz ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R2_001.fastq.gz
			FR2="Z"
		else
			FR2="X"
		fi
	else
		FR2="G"
	fi
	#echo "${FR2}"
	
	#echo ${processed}/${project}/${sample_name}/trimmed/${sample_name}_R1_001.paired.fq.gz
	if [[ ! -f ${processed}/${project}/${sample_name}/trimmed/${sample_name}_R1_001.paired.fq.gz ]]; then
		if [[ -f ${processed}/${project}/${sample_name}/trimmed/${sample_name}_R1.paired.fq.gz ]]; then
			mv ${processed}/${project}/${sample_name}/trimmed/${sample_name}_R1.paired.fq.gz ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R1_001.paired.fq.gz
			TR1="R"
		elif [[ -f ${processed}/${project}/${sample_name}/trimmed/${sample_name}_R1_001.paired.fq ]]; then
			gzip ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R1_001.paired.fq
			TR1="Z"
		elif [[ -f ${processed}/${project}/${sample_name}/trimmed/${sample_name}_R1.paired.fq ]]; then
			mv ${processed}/${project}/${sample_name}/trimmed/${sample_name}_R1.paired.fq ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R1_001.paired.fq
			gzip ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R1_001.paired.fq
			TR1="Y"
		fi
	else
		TR1="G"
	fi
	#echo "${TR1}"
	
	#echo ${processed}/${project}/${sample_name}/trimmed/${sample_name}_R2.paired.fq.gz
	if [[ ! -f ${processed}/${project}/${sample_name}/trimmed/${sample_name}_R2_001.paired.fq.gz ]]; then
		if [[ -f ${processed}/${project}/${sample_name}/trimmed/${sample_name}_R2.paired.fq.gz ]]; then
			mv ${processed}/${project}/${sample_name}/trimmed/${sample_name}_R2.paired.fq.gz ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R2_001.paired.fq.gz
			TR2="R"
		elif [[ -f ${processed}/${project}/${sample_name}/trimmed/${sample_name}_R2_001.paired.fq ]]; then
			gzip ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R2_001.paired.fq
			TR2="Z"
		elif [[ -f ${processed}/${project}/${sample_name}/trimmed/${sample_name}_R2.paired.fq ]]; then
			mv ${processed}/${project}/${sample_name}/trimmed/${sample_name}_R2.paired.fq ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R2_001.paired.fq
			gzip ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R2_001.paired.fq
			TR2="Y"
		fi
	else
		TR2="G"
	fi
	#echo "${TR2}"
	
	if [[ ! -f ${processed}/${project}/${sample_name}/Assembly_Stats/${sample_name}_report.txt ]]; then
		mv ${processed}/${project}/${sample_name}/Assembly_Stats/report.txt ${processed}/${project}/${sample_name}/Assembly_Stats/${sample_name}_report.txt 
	fi
	if [[ ! -f ${processed}/${project}/${sample_name}/Assembly_Stats/${sample_name}_report.tsv ]]; then
		mv ${processed}/${project}/${sample_name}/Assembly_Stats/report.tsv ${processed}/${project}/${sample_name}/Assembly_Stats/${sample_name}_report.tsv
	fi
	if [[ ! -f ${processed}/${project}/${sample_name}/Assembly_Stats/${sample_name}_toms_assembly_report.txt ]]; then
		mv ${processed}/${project}/${sample_name}/Assembly_Stats/toms_assembly_report.txt ${processed}/${project}/${sample_name}/Assembly_Stats/${sample_name}_toms_assembly_report.txt
	fi
	echo "${project}/${sample_name} ${FQ1} ${FQ2} ${FR1} ${FR2} ${TR1} ${TR2}" >> ${share}/Rs.out
done < "${share}/${1}"
echo "All isolates completed"
global_end_time=$(date "+%m-%d-%Y @ %Hh_%Mm_%Ss")
#Script exited gracefully (unless something else inside failed)
printf "%s %s" "Act_by_list.sh has completed ${2} " "${global_end_time}" | mail -s "act_by_list complete" nvx4@cdc.gov
exit 0
