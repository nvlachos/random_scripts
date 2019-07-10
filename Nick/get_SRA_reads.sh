#!/bin/sh -l

#$ -o getSRA.out
#$ -e getSRA.err
#$ -N getSRA
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi
. ./config.sh

#
# Script to download SRA reads
#
# Usage ./get_SRA_reads.sh SRA_Number output_directory (will mirror standard quaisar-h output format)
#
#  sratoolkit/2.9.1
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
elif [[ -z "${1}" ]]; then
	echo "Empty SRA Number supplied to $0, exiting"
	exit 1
# Gives the user a brief usage and help section if requested with the -h option argument
elif [[ "${1}" = "-h" ]]; then
	echo "Usage is ./get_SRA_reads.sh   SRA_Number   run_id"
	echo "Output is saved to ${processed}/[References|run_id]/SRA_Number/"
	exit 0
elif [[ -z "${2}" ]]; then
	echo "No output directory given, using default ${processed}..."
	OUTDATADIR="${processed}/References/${1}"
else
	OUTDATADIR="${2}/${1}"
fi

if [ ! -d "${OUTDATADIR}" ]; then
	echo "Empty output directory supplied to get_SRA_reads.sh...creating"
	mkdir -p "${OUTDATADIR}/FASTQs"
fi



module load sratoolkit/2.9.1
module load BBMap/38.26
module load trimmomatic/0.35

fasterq-dump --split-files ${1} -O ${OUTDATADIR}/FASTQs

complete="true"
if [[ -s "${OUTDATADIR}/FASTQs/${1}_1.fastq" ]]; then
	gzip -c "${OUTDATADIR}/FASTQs/${1}_1.fastq" > "${OUTDATADIR}/FASTQs/${1}_R1_001.fastq.gz"
else
	echo "R1 does not exist for ${1}"
	complete="false"
fi
if [[ -s "${OUTDATADIR}/FASTQs/${1}_2.fastq" ]]; then
	gzip -c "${OUTDATADIR}/FASTQs/${1}_2.fastq" > "${OUTDATADIR}/FASTQs/${1}_R2_001.fastq.gz"
else
	echo "R2 does not exist for ${1}"
	complete="false"
fi

if [[ "${complete}" == "true" ]]; then
	bbduk.sh -"${bbduk_mem}" threads="${procs}" in="${OUTDATADIR}/FASTQs/${1}_R1_001.fastq" in2="${OUTDATADIR}/FASTQs/${1}_R2_001.fastq" out="${OUTDATADIR}/removedAdapters/${1}-noPhiX-R1.fsq" out2="${OUTDATADIR}/removedAdapters/${1}-noPhiX-R2.fsq" ref="${phiX_location}" k="${bbduk_k}" hdist="${bbduk_hdist}"
	mkdir ${OUTDATADIR}/trimmed
	trimmomatic "${trim_endtype}" -"${trim_phred}" -threads "${procs}" "${OUTDATADIR}/removedAdapters/${1}-noPhiX-R1.fsq" "${OUTDATADIR}/removedAdapters/${1}-noPhiX-R2.fsq" "${OUTDATADIR}/trimmed/${1}_R1_001.paired.fq" "${OUTDATADIR}/trimmed/${1}_R1_001.unpaired.fq" "${OUTDATADIR}/trimmed/${1}_R2_001.paired.fq" "${OUTDATADIR}/trimmed/${1}_R2_001.unpaired.fq" ILLUMINACLIP:"${trim_adapter_location}:${trim_seed_mismatch}:${trim_palindrome_clip_threshold}:${trim_simple_clip_threshold}:${trim_min_adapt_length}:${trim_complete_palindrome}" SLIDINGWINDOW:"${trim_window_size}:${trim_window_qual}" LEADING:"${trim_leading}" TRAILING:"${trim_trailing}" MINLEN:"${trim_min_length}"
	cat "${OUTDATADIR}/trimmed/${1}_R1_001.paired.fq" "${OUTDATADIR}/trimmed/${1}_R2_001.paired.fq" > "${OUTDATADIR}/trimmed/${1}.paired.fq"
	cat "${OUTDATADIR}/trimmed/${1}_R1_001.unpaired.fq" "${OUTDATADIR}/trimmed/${1}_R2_001.unpaired.fq" > "${OUTDATADIR}/trimmed/${1}.single.fq"
fi

module unload sratoolkit/2.9.1
module unload BBMap/38.26
module unload trimmomatic/0.35

exit 0
