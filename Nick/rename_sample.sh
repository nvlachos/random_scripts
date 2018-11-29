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

old_name=${1}
if [[ ! -z ${4} ]]; then
	old_name=${4}
fi
new_name=${2}
project=${3}
echo "Test 1-${old_name}-${new_name}, filename changing"
find /scicomp/groups/OID/NCEZID/DHQP/CEMB/MiSeqAnalysisFiles/${project}/${old_name} -type f -name "*${old_name}*" | while read FILE ; do
  dirname=$(dirname $FILE)
	filename=$(basename $FILE)
	#echo "Found-${FILE}"
	#echo "$dirname"
	#echo "$filename"
	newfile="$(echo ${filename/${old_name}/${new_name}})"
	echo "${newfile}"
    mv "${FILE}" "${dirname}/${newfile}"
done
echo "Test 2, deep dive"
for thing in /scicomp/groups/OID/NCEZID/DHQP/CEMB/MiSeqAnalysisFiles/${project}/${old_name}/*; do
	if [[ -f ${thing} ]]; then
		sed -i -- 's/${old_name}/${new_name}/g' ${thing}
	elif [[ -d ${thing} ]]; then
		if [[ "${thing}" = *"FASTQs" ]]; then
			echo "Skipping FASTQs deep dive"
		fi
		find ${thing} -type f -exec sed -i 's/${old_name}/${new_name}/g' {} +
	else
		echo "Thing (${thing}) is not file or directory"
	fi
done
mv "/scicomp/groups/OID/NCEZID/DHQP/CEMB/MiSeqAnalysisFiles/${project}/${old_name}" "/scicomp/groups/OID/NCEZID/DHQP/CEMB/MiSeqAnalysisFiles/${project}/${new_name}"
rm -r "/scicomp/groups/OID/NCEZID/DHQP/CEMB/MiSeqAnalysisFiles/${project}/temp_folders.txt"
