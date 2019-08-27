#!/bin/sh -l

#$ -o AR_updater.out
#$ -e AR_updater.err
#$ -N AR_updater
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
. ./config.sh
#Import the module file that loads all necessary mods
. "${mod_changers}/pipeline_mods"

#
# Usage ./AR_updater.sh
#
#
#

# Checks for proper argumentation
if [[ "$1" = "-h" ]]; then
	echo "Usage is ./AR-updater.sh"
	exit 0
fi

# Makes list of ALL samples in MiSeqAnalysisFiles
today=$(date "+%Y%m%d")
${shareScript}/make_sample_list_from_run_list.sh ${shareScript}/sample_list_${today}.txt ${shareScript}/directory_list_${today}.txt

# Exits if a sample list from TODAY was not creaetd
if [[ ! -f ${shareScript}/sample_list_${today}.txt ]]; then
	echo "No sample list file was created (${shareScript}/sample_list_${today}.txt), must exit..."
	exit 245
fi

# Checks the last time isolates were updated for AR. If databases havent changed since then, then no need to do any updating
last_run=$(tail -n1 ${local_DBs}/ar_updater.txt)
echo "Last time run=${last_run}
ResGANNCBI created=${ResGANNCBI_srst2_filename}"
if [[ "${last_run}" == "${ResGANNCBI_srst2_filename}" ]]; then
	echo "The last run was with the newest database... ${ResGANNCBI_srst2_filename}, exiting"
	exit 0
fi

# Gets the date that the ResGANNCBI database was made
ResGANNCBI_date=$(echo ${ResGANNCBI_srst2_filename} | cut -d'_' -f2)

# Submit the list of samples for csstar and srst2
qsub "${shareScript}/abl_mass_qsub_csstar.sh" "${shareScript}/sample_list_${today}.txt" 100 "/scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/mass_subs" "keep"
qsub -sync y "${shareScript}/abl_mass_qsub_srst2.sh" "${shareScript}/sample_list_${today}.txt" 100 "/scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/mass_subs" "keep"
qsub "${shareScript}/abl_mass_qsub_runsum.sh" "${shareScript}/directory_list_${today}.txt" 100 "/scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/mass_subs"
qsub -sync y "${shareScript}/abl_AR_completion_check.sh" "${shareScript}/sample_list_${today}.txt" "${ResGANNCBI_date}" "/scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/AR_check_${today}.txt"



echo "All isolates are UTD on AR files, check /scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/AR_check_${today} to see if anything failed or is missing"
echo "${ResGANNCBI_srst2_filename}" >> ${local_DBs}/ar_updater.txt

global_end_time=$(date "+%m-%d-%Y @ %Hh_%Mm_%Ss")
rm ${shareScript}/sample_list_${today}.txt
rm ${shareScript}/directory_list_${today}.txt
#Script exited gracefully (unless something else inside failed)
printf "%s %s" "AR_updater.sh has completed ${1}" "${global_end_time}" | mail -s "AR_updater complete" nvx4@cdc.gov
exit 0
