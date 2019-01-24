#!/bin/sh -l

#$ -o rename_sample.out
#$ -e rename_sample.err
#$ -N rename_sample
#$ -cwd
#$ -q all.q

#Import the config file with shortcuts and settings
. ./config.sh
if [[ "$1" = "-h" ]]; then
	echo "Usage is ./rename_sample.sh old_name new_name project_id"
fi

sample_name=${1}
if [[ ! -z ${4} ]]; then
	old_name=${4}
fi
new_project_name=${3}
old_project_name=${2}

echo "Test 1, Searching contents of files"
for thing in /scicomp/groups/OID/NCEZID/DHQP/CEMB/MiSeqAnalysisFiles/${old_project_name}/${sample_name}/*; do
	if [[ -f ${thing} ]]; then
		echo "Doing normal - $thing"
		sed -i "s/${old_project_name}/${new_project_name}/g" ${thing}
	elif [[ -d ${thing} ]]; then
		echo "doing directory - $thing"
		if [[ "${thing}" = *"FASTQs" ]]; then
			echo "Skipping FASTQs deep dive"
		fi
		find ${thing} -type f -exec sed -i "s/${old_project_name}/${new_project_name}/g" {} +
	else
		echo "Thing (${thing}) is not file or directory"
	fi
done

# echo "Test 2, ${old_project_name}-${new_project_name}, Searching filenames"
# find /scicomp/groups/OID/NCEZID/DHQP/CEMB/MiSeqAnalysisFiles/${project}/${old_name} -type f -name "*${old_name}*" | while read FILE ; do
#   dirname=$(dirname $FILE)
# 	filename=$(basename $FILE)
# 	#echo "Found-${FILE}"
# 	#echo "$dirname"
# 	#echo "$filename"
# 	newfile="$(echo ${filename/${old_name}/${new_name}})"
# 	echo "${newfile}"
#     mv "${FILE}" "${dirname}/${newfile}"
# done

echo "Checking summary and list files"
for thing in /scicomp/groups/OID/NCEZID/DHQP/CEMB/MiSeqAnalysisFiles/${old_project_name}/*; do
	if [[ -f ${thing} ]]; then
		echo "Doing normal - $thing"
		sed -i "s/${old_project_name}/${new_project_name}/g" ${thing}
	elif [[ -d ${thing} ]]; then
		echo "doing nothing on directory - $thing"
	else
		echo "Thing (${thing}) is not file or directory"
	fi
done

mv "/scicomp/groups/OID/NCEZID/DHQP/CEMB/MiSeqAnalysisFiles/${old_project_name}" "/scicomp/groups/OID/NCEZID/DHQP/CEMB/MiSeqAnalysisFiles/${new_project_name}"
