#!/bin/sh -l

#$ -o rename_project.out
#$ -e rename_project.err
#$ -N rename_project
#$ -cwd
#$ -q all.q

#Import the config file with shortcuts and settings
. ./config.sh
if [[ "$1" = "-h" ]]; then
	echo "Usage is ./rename_sample.sh old_project_name new_project_name"
fi

new_project_name=${2}
old_project_name=${1}


# Finding all internal instances of old project name and changing them to new preoject name
echo "Testing new internal finder"
#find . -not -name "*.fq" -not -name "*.fastq" -not -name "*.fsq" -not -name "*.fasta" -type f -print0 | xargs -0 sed -i 's/${old_project_name}/${new_project_name}/g'

grep -lr --exclude-dir="removedAdapters" --exclude-dir="trimmed" --exclude-dir="FASTQs" --exclude-dir="ANI" -e "${old_project_name}" ${processed}/${old_project_name} | xargs sed -i '' -e 's/${old_project_name}/${new_project_name}/g'
#grep -lr --exclude-dir="removedAdapters" --exclude-dir="trimmed" --exclude-dir="FASTQs" --exclude-dir="ANI" -e "${old_project_name}" ${processed}/${old_project_name} > ${processed}/${old_project_name}/groutput.txt

# for thing in ${processed}/${old_project_name}/*; do
# 	if [[ -d ${thing} ]]; then
# 		if [[ "${thing}" = *"FASTQs" ]] || [[ "${thing}" = *"removedAdapters" ]] || [[ "${thing}" = *"trimmed" ]]; then
# 			continue
# 		else
# 			echo "Trying ${thing}"
# 			find ${thing} -type f -exec sed -i 's/${old_project_name}/${new_sample_name}/g' {} +
# 		fi
# 	fi
# done

#echo "Test 1, Searching contents of files"
#for thing in /scicomp/groups/OID/NCEZID/DHQP/CEMB/MiSeqAnalysisFiles/${old_project_name}/*; do
#	if [[ -f ${thing} ]]; then
#		echo "Doing normal - $thing"
#	elif [[ -d ${thing} ]]; then
#		echo "doing directory - $thing"
#		if [[ "${thing}" = *"FASTQs" ]]; then
#			echo "Skipping FASTQs deep dive"
#		fi
#		find ${thing} -type f -exec sed -i "s/${old_project_name}/${new_project_name}/g" {} +
#	else
#		echo "Thing (${thing}) is not file or directory"
#	fi
#done


# Changing all files names containing old project name
echo "Testing new filename changer"


for i in ${processed}/${old_project_name}/* ;do mv -- "$i" "${i//${old_project_name}/${new_project_name}}";done

#for thing in ${processed}/${old_project_name}/*; do
#	if [[ "${thing}" = *"${old_project_name}"* ]]; then
#		rename ${old_project_name} ${new_project_name} ${thing}
#	fi
#done



mv "${processed}/${old_project_name}" "/${processed}/${new_project_name}"
