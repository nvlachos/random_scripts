#!/bin/sh -l

#$ -o check_FASTQs.out
#$ -e check_FASTQs.err
#$ -N ablb2
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
. ./config.sh
#Import the module file that loads all necessary mods
. "${mod_changers}/pipeline_mods"

#
# Usage ./check_FASTQs.sh path_to_list Description of list function
#
# script changes depending on what needs to be run through the list
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./check_FASTQs.sh path_to_list_file(single sample ID per line, e.g.B8VHY/1700128 (it must include project id also)) output_file"
	echo "Output location varies depending on which tasks are performed but will be found somewhere under ${processed}"
	exit 0
elif [[ -z "${2}" ]]	; then
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
			FQ1="R1-Good"
		elif [[ -f ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R1.fastq ]]; then
			mv ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R1.fastq ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R1_001.fastq
			gzip ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R1_001.fastq
			FQ1="R1-Renamed-and-zipped"
		fi
	else
		if [[ -f ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R1_001.fastq.gz ]]; then
			rm -r ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R1_001.fastq
			FQ1="R1-Both-removed-unzipped"
		else
			gzip ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R1_001.fastq
			FQ1="R1-zipped"
		fi
	fi
	#echo "${FQ1}"

	#echo ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R2_001.fastq
	if [[ ! -f ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R2_001.fastq ]]; then
		if [[ -f ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R2_001.fastq.gz ]]; then
			FQ2="R2-Good"
		elif [[ -f ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R2.fastq ]]; then
			mv ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R2.fastq ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R2_001.fastq
			gzip ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R2_001.fastq
			FQ2="R2-Renamed-and-zipped"
		fi
	else
		if [[ -f ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R2_001.fastq.gz ]]; then
			rm -r ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R2_001.fastq
			FQ2="R2-Both-removed-unzipped"
		else
			gzip ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R2_001.fastq
			FQ2="R2-zipped"
		fi
	fi
	#echo "${FQ2}"

	#echo ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R1_001.fastq.gz
	if [[ ! -f ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R1_001.fastq.gz ]]; then
		if [[ -f ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R1.fastq.gz ]]; then
			mv ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R1.fastq.gz ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R1_001.fastq.gz
			FR1="R1-Renamed"
		elif [[ -f ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R1.fastq ]]; then
			gzip ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R1.fastq
			mv ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R1.fastq.gz ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R1_001.fastq.gz
			FR1="R1-Renamed-and-zipped"
		else
			FR1="No-R1-found"
		fi
	else
		FR1="R1-Good"
	fi
	#echo "${FR1}"

	#echo ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R2_001.fastq.gz
	if [[ ! -f ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R2_001.fastq.gz ]]; then
		if [[ -f ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R2.fastq.gz ]]; then
			mv ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R2.fastq.gz ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R2_001.fastq.gz
			FR2="R2-Renamed"
		elif [[ -f ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R2.fastq ]]; then
			gzip ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R2.fastq
			mv ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R2.fastq.gz ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R2_001.fastq.gz
			FR2="R2-Renamed-and-zipped"
		else
			FR2="No-R2-found"
		fi
	else
		FR2="R2-Good"
	fi
	#echo "${FR2}"

	#echo ${processed}/${project}/${sample_name}/trimmed/${sample_name}_R1_001.paired.fq.gz
	if [[ ! -f ${processed}/${project}/${sample_name}/trimmed/${sample_name}_R1_001.paired.fq.gz ]]; then
		if [[ -f ${processed}/${project}/${sample_name}/trimmed/${sample_name}_R1.paired.fq.gz ]]; then
			mv ${processed}/${project}/${sample_name}/trimmed/${sample_name}_R1.paired.fq.gz ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R1_001.paired.fq.gz
			TR1="R1-Paired-Renamed"
		elif [[ -f ${processed}/${project}/${sample_name}/trimmed/${sample_name}_R1_001.paired.fq ]]; then
			gzip ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R1_001.paired.fq
			TR1="R1-PairedZipped"
		elif [[ -f ${processed}/${project}/${sample_name}/trimmed/${sample_name}_R1.paired.fq ]]; then
			mv ${processed}/${project}/${sample_name}/trimmed/${sample_name}_R1.paired.fq ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R1_001.paired.fq
			gzip ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R1_001.paired.fq
			TR1="R1-Paired-Renamed-and-zipped"
		else
			TR1="R1-Paired-Not-found"
		fi
	else
		TR1="R1-Paired-Good"
	fi
	#echo "${TR1}"

	#echo ${processed}/${project}/${sample_name}/trimmed/${sample_name}_R2.paired.fq.gz
	if [[ ! -f ${processed}/${project}/${sample_name}/trimmed/${sample_name}_R2_001.paired.fq.gz ]]; then
		if [[ -f ${processed}/${project}/${sample_name}/trimmed/${sample_name}_R2.paired.fq.gz ]]; then
			mv ${processed}/${project}/${sample_name}/trimmed/${sample_name}_R2.paired.fq.gz ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R2_001.paired.fq.gz
			TR2="R2-Paired-Renamed"
		elif [[ -f ${processed}/${project}/${sample_name}/trimmed/${sample_name}_R2_001.paired.fq ]]; then
			gzip ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R2_001.paired.fq
			TR2="R2-Paired-Zipped"
		elif [[ -f ${processed}/${project}/${sample_name}/trimmed/${sample_name}_R2.paired.fq ]]; then
			mv ${processed}/${project}/${sample_name}/trimmed/${sample_name}_R2.paired.fq ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R2_001.paired.fq
			gzip ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R2_001.paired.fq
			TR2="R2-Paired-Renamed-and-zipped"
		else
			TR1="R2-Paired-Not-found"
		fi
	else
		TR2="R2-Paired-Good"
	fi
	#echo "${TR2}"

	if [[ ! -f ${processed}/${project}/${sample_name}/Assembly_Stats/${sample_name}_report.txt ]]; then
		mv ${processed}/${project}/${sample_name}/Assembly_Stats/report.txt ${processed}/${project}/${sample_name}/Assembly_Stats/${sample_name}_report.txt
	fi
	if [[ ! -f ${processed}/${project}/${sample_name}/Assembly_Stats/${sample_name}_report.tsv ]]; then
		mv ${processed}/${project}/${sample_name}/Assembly_Stats/report.tsv ${processed}/${project}/${sample_name}/Assembly_Stats/${sample_name}_report.tsv
	fi

	echo "${project}/${sample_name} ${FQ1} ${FQ2} ${FR1} ${FR2} ${TR1} ${TR2}" >> ${2}
done < "${1}"
echo "All isolates completed"
global_end_time=$(date "+%m-%d-%Y @ %Hh_%Mm_%Ss")
#Script exited gracefully (unless something else inside failed)
printf "%s %s" "check_FASTQs.sh has completed ${2} " "${global_end_time}" | mail -s "check_FASTQs complete" nvx4@cdc.gov
exit 0
