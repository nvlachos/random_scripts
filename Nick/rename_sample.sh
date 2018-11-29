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
echo "Test 1-${old_name}-${new_name}"
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

folders=()
index=0
find /scicomp/groups/OID/NCEZID/DHQP/CEMB/MiSeqAnalysisFiles/${project}/${old_name} -type d > "/scicomp/groups/OID/NCEZID/DHQP/CEMB/MiSeqAnalysisFiles/${project}/temp_folders.txt"
echo "Test 2"
find /scicomp/groups/OID/NCEZID/DHQP/CEMB/MiSeqAnalysisFiles/${project}/${old_name} -type f -exec sed -i 's/{old_name}/${new_name}/g' {} +
echo "Test 3"
for (( idx=${index}-1 ; idx>=0 ; idx-- )) ; do
	subpath=$(echo "${folders[${idx}]}" | cut -d'/' -f1-10)
	sub_sample=$(echo "${folders[${idx}]}" | cut -d'/' -f11-)
	#echo "sp-${subpath}"
	#echo "ss-$sub_sample"
    new_dir="${subpath}/${sub_sample/${old_name}/${new_name}}"
	echo "NT-${new_dir}"
	if [[ "${new_dir}" != "${folders[${idx}]}" ]]; then
		mv "${folders[${idx}]}" "${new_dir}"
	else
		echo "Destination is the same as the source"
	fi
done
mv "/scicomp/groups/OID/NCEZID/DHQP/CEMB/MiSeqAnalysisFiles/${project}/${old_name}" "/scicomp/groups/OID/NCEZID/DHQP/CEMB/MiSeqAnalysisFiles/${project}/${new_name}"
rm -r "/scicomp/groups/OID/NCEZID/DHQP/CEMB/MiSeqAnalysisFiles/${project}/temp_folders.txt"
