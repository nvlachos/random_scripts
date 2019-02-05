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
directory_to_scan=${3}
 python3 "${shareScript}/rename_project.py" ${1} ${2} ${3}
