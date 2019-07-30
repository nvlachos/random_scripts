#!/bin/sh -l

#$ -o abl-template.out
#$ -e abl-template.err
#$ -N abl-template
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
. ./config.sh
#Import the module file that loads all necessary mods
#. "${mod_changers}/pipeline_mods"

#
# Usage ./act_by_list_template.sh path_to_list path_for_output_file
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./act_by_list_template.sh path_for_list_file"
	exit 0
fi

#ml BBMap/38.26 trimmomatic/0.35

# Loop through and act on each sample name in the passed/provided list
while IFS= read -r var; do
	sample_name=$(echo "${var}" | cut -d'/' -f2 | tr -d '[:space:]')
	project=$(echo "${var}" | cut -d'/' -f1 | tr -d '[:space:]')
	#bbduk.sh -"${bbduk_mem}" threads="${procs}" in="${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R1_001.fastq" in2="${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R2_001.fastq" out="${processed}/${project}/${sample_name}/removedAdapters/${sample_name}-noPhiX-R1.fsq" out2="${processed}/${project}/${sample_name}/removedAdapters/${sample_name}-noPhiX-R2.fsq" ref="${phiX_location}" k="${bbduk_k}" hdist="${bbduk_hdist}"
	#trimmomatic "${trim_endtype}" -"${trim_phred}" -threads "${procs}" "${processed}/${project}/${sample_name}/removedAdapters/${sample_name}-noPhiX-R1.fsq" "${processed}/${project}/${sample_name}/removedAdapters/${sample_name}-noPhiX-R2.fsq" "${processed}/${project}/${sample_name}/trimmed/${sample_name}_R1_001.paired.fq" "${processed}/${project}/${sample_name}/trimmed/${sample_name}_R1_001.unpaired.fq" "${processed}/${project}/${sample_name}/trimmed/${sample_name}_R2_001.paired.fq" "${processed}/${project}/${sample_name}/trimmed/${sample_name}_R2_001.unpaired.fq" ILLUMINACLIP:"${trim_adapter_location}:${trim_seed_mismatch}:${trim_palindrome_clip_threshold}:${trim_simple_clip_threshold}:${trim_min_adapt_length}:${trim_complete_palindrome}" SLIDINGWINDOW:"${trim_window_size}:${trim_window_qual}" LEADING:"${trim_leading}" TRAILING:"${trim_trailing}" MINLEN:"${trim_min_length}"
	#rm -r "${processed}/${project}/${sample_name}/removedAdapters"
	if [[ ! -f "${processed}/${project}/${sample_name}/trimmed/${sample_name}.paired.fq" ]]; then
		cat "${processed}/${project}/${sample_name}/trimmed/${sample_name}_R1_001.paired.fq" "${processed}/${project}/${sample_name}/trimmed/${sample_name}_R2_001.paired.fq" > "${processed}/${project}/${sample_name}/trimmed/${sample_name}.paired.fq"
	fi
	echo "About to zip - ${processed}/${project}/${sample_name}/trimmed/${sample_name}.paired.fq"
	gzip "${processed}/${project}/${sample_name}/trimmed/${sample_name}.paired.fq"
	echo "About to zip - ${processed}/${project}/${sample_name}/trimmed/${sample_name}_R1_001.paired.fq"
	gzip "${processed}/${project}/${sample_name}/trimmed/${sample_name}_R1_001.paired.fq"
	echo "About to zip - ${processed}/${project}/${sample_name}/trimmed/${sample_name}_R2_001.paired.fq"
	gzip "${processed}/${project}/${sample_name}/trimmed/${sample_name}_R2_001.paired.fq"
done < "${1}"

#ml -BBMap/38.26 -trimmomatic/0.35

echo "All isolates completed"
global_end_time=$(date "+%m-%d-%Y @ %Hh_%Mm_%Ss")
#Script exited gracefully (unless something else inside failed)
printf "%s %s" "Act_by_list_template.sh has completed check of snd MLSTs" "${global_end_time}" | mail -s "act_by_list complete" nvx4@cdc.gov
exit 0
