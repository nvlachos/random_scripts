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
	echo "Usage is ./get_SRA_reads.sh   sample_name   run_ID"
	echo "Output is saved to ${processed}/References/SRA_Number/"
	exit 0
elif [[ -z "${2}" ]]; then
	echo "No output directory given, using default ${processed}..."
	OUTDATADIR="${processed}/References/${1}"
else
	OUTDATADIR="${2}/${1}"
fi

if [ ! -d "${OUTDATADIR}" ]; then
	echo "Empty output directory supplied to $0...creating"
	mkdir -p "${OUTDATADIR}/FASTQs"
fi



module load sratoolkit/2.9.1
module load BBMap/38.26
module load trimmomatic/0.35

cd ${OUTDATADIR}/FASTQs

sra_id=$(esearch -db sra -query ${1} </dev/null | efetch -format docsum | xtract -pattern Runs -element Run@acc)
fastq-dump -A $sra_id --split-files --origfmt --gzip

complete="true"
if [[ -s "${OUTDATADIR}/FASTQs/${1}_1.fastq.gz" ]]; then
	mv "${OUTDATADIR}/FASTQs/${1}_1.fastq.gz" "${OUTDATADIR}/FASTQs/${1}_R1_001.fastq.gz"
else
	echo "R1 does not exist for ${1}"
	complete="false"
fi
if [[ -s "${OUTDATADIR}/FASTQs/${1}_2.fastq.gz" ]]; then
	mv "${OUTDATADIR}/FASTQs/${1}_2.fastq.gz" "${OUTDATADIR}/FASTQs/${1}_R2_001.fastq.gz"
else
	echo "R2 does not exist for ${1}"
	complete="false"
fi

if [[ "${complete}" == "true" ]]; then
	bbduk.sh -"${bbduk_mem}" threads="${procs}" in="${OUTDATADIR}/FASTQs/${1}_R1_001.fastq" in2="${OUTDATADIR}/FASTQs/${1}_R2_001.fastq" out="${OUTDATADIR}/removedAdapters/${1}-noPhiX-R1.fsq" out2="${OUTDATADIR}/removedAdapters/${1}-noPhiX-R2.fsq" ref="${phiX_location}" k="${bbduk_k}" hdist="${bbduk_hdist}"
	trimmomatic "${trim_endtype}" -"${trim_phred}" -threads "${procs}" "${OUTDATADIR}/removedAdapters/${1}-noPhiX-R1.fsq" "${OUTDATADIR}/removedAdapters/${1}-noPhiX-R2.fsq" "${OUTDATADIR}/trimmed/${1}_R1_001.paired.fq" "${OUTDATADIR}/trimmed/${1}_R1_001.unpaired.fq" "${OUTDATADIR}/trimmed/${1}_R2_001.paired.fq" "${OUTDATADIR}/trimmed/${1}_R2_001.unpaired.fq"
	cat "${OUTDATADIR}/trimmed/${1}_R1_001.paired.fq" "${OUTDATADIR}/trimmed/${1}_R2_001.paired.fq" > "${OUTDATADIR}/trimmed/${1}.paired.fq"
	cat "${OUTDATADIR}/trimmed/${1}_R1_001.unpaired.fq" "${OUTDATADIR}/trimmed/${1}_R2_001.unpaired.fq" > "${OUTDATADIR}/trimmed/${1}.single.fq"
fi

module unload sratoolkit/2.9.1
module unload BBMap/38.26
module unload trimmomatic/0.35

exit 0
