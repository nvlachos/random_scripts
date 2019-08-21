#!/bin/sh -l

#$ -o bionumerics_get_reads_MAFF.out
#$ -e bionumerics_get_reads_MAFF.err
#$ -N bionumerics_get_reads_MAFF
#$ -cwd
#$ -q short.q

# Import the config file with shortcuts and settings
. ./config.sh

#
# Will find all fastq.gz files matching samples in provided list file
#
# Usage ./bionumerics_get_reads.sh  destination_database_folder(Genus_species) list_file (Maybe add sample name or list functionality later)
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
elif [[ -z "${1}" ]]; then
	echo "Empty project name supplied to $0, exiting"
	exit 1
elif [[ "${1}" = "-h" ]]; then
	echo "Usage is ./bionumerics_get_reads.sh  Genus_species file_list (use extra / and machine id if data has been moved instruments"
	exit 0
elif [[ -z "${2}" ]]; then
	echo "Empty file list supplied to $0, exiting"
	exit 1
fi


# Sets folder to where files will be downloaded to
OUTDATADIR="/scicomp/groups/OID/NCEZID/DHQP/CEMB/analysis/calcengine/${1}"

# Loop through and act on each sample name in the passed/provided list
while IFS= read -r var; do
	sample_name=$(echo "${var}" | awk -F/ '{ print $2}' | tr -d '[:space:]')
	project=$(echo "${var}" | awk -F/ '{ print $1}' | tr -d '[:space:]')
	if [[ "${full_name}" == "_" ]]; then
		echo "${sample_name} not found on instrument: ${full_instrument}"
		continue
	else
		echo "Attempting to copy: ${sample_name}"
		if [[ ! -s "${OUTDATADIR}/${sample_name}_R1_001.fastq" ]]; then
#			gunzip -c "${full_name1}" > "${OUTDATADIR}/${sample_name}_R1_001.fastq"
			echo "Unzipping R1-${sample_name}"
			gunzip -c "${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R1_001.fastq.gz" > "${OUTDATADIR}/${sample_name}_R1_001.fastq"
		else
			echo "R1 already exists-${sample_name}"
		fi
		if [[ ! -s "${OUTDATADIR}/${samlpe_name}_R2_001.fastq" ]]; then
#			gunzip -c "${full_name2}" > "${OUTDATADIR}/${sample_name}_R2_001.fastq"
			echo "Unzipping R2-${sample_name}"
			gunzip -c "${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R2_001.fastq.gz" > "${OUTDATADIR}/${sample_name}_R2_001.fastq"
		else
			echo "R2 already exists-${sample_name}"
		fi
	fi
done < ${2}

#Script exited gracefully (unless something else inside failed)
exit 0
