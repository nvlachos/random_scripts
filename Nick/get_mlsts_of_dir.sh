#!/bin/sh -l

#$ -o abl-gmlsts.out
#$ -e abl-gmlsts.err
#$ -N abl-gmlsts
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
. ./config.sh
#Import the module file that loads all necessary mods
. "${mod_changers}/pipeline_mods"

#
# Usage ./get_mlsts_of_dir.sh path_to_list path_for_output_file
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./get_mlsts_of_dir.sh path_to_list_file ResGANNOT_Identifier(YYYYMMDD) path_for_output_file"
	exit 0
elif [[ ! -d "${1}" ]]; then
	echo  "No folder, exiting..."
	exit 1
elif [[ ! -f "${2}" ]]; then
	>${2}
	echo "MLSTs from list"
fi


# Loop through and act on each sample name in the passed/provided list
counter=0
for assembly in "${1}/*"; do
	if [[ ]]
	mlst "${assembly}" > "${OUTDATADIR}/MLST/${1}.mlst"
done

echo "All isolates completed"
global_end_time=$(date "+%m-%d-%Y @ %Hh_%Mm_%Ss")
#Script exited gracefully (unless something else inside failed)
printf "%s %s" "get_mlsts_of_dir.sh has completed check of snd MLSTs" "${global_end_time}" | mail -s "act_by_list complete" nvx4@cdc.gov
exit 0
