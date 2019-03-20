#!/bin/sh -l

#$ -o run_MLST.out
#$ -e run_MLST.err
#$ -N run_MLST
#$ -cwd
#$ -q all.q

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
	echo "No argument supplied to get_SRA_reads.sh, exiting"
	exit 1
elif [[ -z "${1}" ]]; then
	echo "Empty SRA Number supplied to get_SRA_reads.sh, exiting"
	exit 1
# Gives the user a brief usage and help section if requested with the -h option argument
elif [[ "${1}" = "-h" ]]; then
	echo "Usage is ./get_SRA_reads.sh   sample_name   run_id"
	echo "Output is saved to ${processed}/References/SRA_Number/"
	exit 0
elif [[ -z "${2}" ]]; then
	echo "No output directory given, using default ${processed}..."
	OUTDATADIR="${processed}/References/${1}/FASTQs"
else;
	OUTDATADIR="${2}/${1}/FASTQs"
fi

if [ ! -d "${2}" ]; then
	echo "Empty output directory supplied to get_SRA_reads.sh...creating"
	mkdir -p "${OUTDATADIR}"
fi



module load sratoolkit/2.9.1

cd ${OUTDATADIR}

fastq-dump --split-files --origfmt --gzip ${1}

if [[ -s "${OUTDATADIR}/${1}_1.fastq.gz" ]]; then
	mv "${OUTDATADIR}/${1}_1.fastq.gz" "${OUTDATADIR}/${1}_R1_001.fastq.gz"
else
	echo "R1 does not exist for ${1}"
fi
if [[ -s "${OUTDATADIR}/${1}_2.fastq.gz" ]]; then
	mv "${OUTDATADIR}/${1}_2.fastq.gz" "${OUTDATADIR}/${1}_R2_001.fastq.gz"
else
	echo "R2 does not exist for ${1}"
fi


module unload sratoolkit/2.9.1
exit 0
