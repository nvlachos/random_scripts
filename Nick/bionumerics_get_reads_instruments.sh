#!/bin/sh -l

#$ -o bionumerics_get_reads.out
#$ -e bionumerics_get_reads.err
#$ -N bionumerics_get_reads
#$ -cwd
#$ -q all.q

# Import the config file with shortcuts and settings
. ./config.sh

#
# Will find all fastq.gz files matching samples in provided list file
#
# Usage ./bionumerics_get_reads.sh  destination_database_folder(Genus_species) list_file (Maybe add sample name or list functionality later)
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to bionumerics_get_reads.sh, exiting"
	exit 1
elif [[ -z "${1}" ]]; then
	echo "Empty project name supplied to bionumerics_get_reads.sh, exiting"
	exit 1
elif [[ "${1}" = "-h" ]]; then
	echo "Usage is ./bionumerics_get_reads.sh  Genus_species file_list (use extra / and machine id if data has been moved instruments"
	exit 0
elif [[ -z "${2}" ]]; then
	echo "Empty file list supplied to bionumerics_get_reads.sh, exiting"
	exit 1
fi


# Sets folder to where files will be downloaded to
OUTDATADIR="/scicomp/groups/OID/NCEZID/DHQP/CEMB/analysis/calcengine/${1}"

# Loop through and act on each sample name in the passed/provided list
while IFS= read -r var; do
	sample_name=$(echo "${var}" | awk -F/ '{ print $2}' | tr -d '[:space:]')
	project=$(echo "${var}" | awk -F/ '{ print $1}' | tr -d '[:space:]')
	alt_machine=$(echo "${var}" | awk -F/ '{ print $3}' | tr -d '[:space:]')
	#echo "${sample_name} in ${project}"
	if [[ ! -z "${alt_machine}" ]]; then
		instrument="${alt_machine}"
	else
		instrument=$(echo "${project}" |  cut -d'_' -f2)
	fi	
	full_instrument=$(find /scicomp/instruments/ -maxdepth 1 -name *${instrument})
	full_instrument=${full_instrument##*/}
	#echo "${full_instrument}"
#	full_name1=$(find ${processed}/${project}/Fastq/ -maxdepth 1 -name ${sample_name}*R1*gz)
#	full_name2=$(find ${processed}/${project}/Fastq/ -maxdepth 1 -name ${sample_name}*R2*gz)
#	full_name=${sample_name}
	full_name=$(find /scicomp/instruments/${full_instrument}/${project}/Data/Intensities/BaseCalls/ -maxdepth 1 -name *${sample_name}*.gz)
	#echo "${full_name}-1"
	full_name=${full_name##*/}
	sample_number=$(echo ${full_name} | cut -d'_' -f2)
	#echo "${sample_number}-2"
	shname=$(echo "${full_name}" | cut -d'_' -f1)
	#echo "${shname}-3"
	full_name="${shname}_${sample_number}"
	if [[ "${full_name}" == "_" ]]; then
		echo "${sample_name} not found on instrument: ${full_instrument}"
		continue
	else
		echo "Attempting to copy: ${shname}"
		if [[ ! -s "${OUTDATADIR}/${shname}_R1_001.fastq" ]]; then
#			gunzip -c "${full_name1}" > "${OUTDATADIR}/${sample_name}_R1_001.fastq"
			gunzip -c "/scicomp/instruments/${full_instrument}/${project}/Data/Intensities/BaseCalls/${full_name}_L001_R1_001.fastq.gz" > "${OUTDATADIR}/${full_name}_R1_001.fastq"
		fi
		if [[ ! -s "${OUTDATADIR}/${shname}_R2_001.fastq" ]]; then
#			gunzip -c "${full_name2}" > "${OUTDATADIR}/${sample_name}_R2_001.fastq"
			gunzip -c "/scicomp/instruments/${full_instrument}/${project}/Data/Intensities/BaseCalls/${full_name}_L001_R2_001.fastq.gz" > "${OUTDATADIR}/${full_name}_R2_001.fastq"
		fi
	fi	
done < /scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/${2}

#Script exited gracefully (unless something else inside failed)
exit 0
