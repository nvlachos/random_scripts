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

echo "Test 1, Searching contents of files"
for thing in /scicomp/groups/OID/NCEZID/DHQP/CEMB/MiSeqAnalysisFiles/${old_project_name}/*; do
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

mv "/scicomp/groups/OID/NCEZID/DHQP/CEMB/MiSeqAnalysisFiles/${old_project_name}" "/scicomp/groups/OID/NCEZID/DHQP/CEMB/MiSeqAnalysisFiles/${new_project_name}"
