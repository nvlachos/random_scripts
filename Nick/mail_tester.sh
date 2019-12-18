#!/bin/sh -l

#$ -o clean_list.out
#$ -e clean_list.err
#$ -N clean_list
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi
. ./config.sh

#
# Script to clean any list file of extra newlines and space
#
# Usage ./clean_list.sh path_to_list_file
#

runsum=$(${shareScript}/view_sum.sh ${1})
outarray+="${runsum}"

# Add print time the run completed in the text that will be emailed

global_end_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")
run_start_time=${global_end_time}
echo "Finished with run ${1} at ${global_end_time}"
outarray+=("
${1} finished at ${global_end_time}")


echo "${outarray}"

echo "Sending summary email to nvx4@cdc.gov"
printf "%s\\n" "${outarray}" | mail -s "Run Status for ${1}_on_${run_start_time}_run.log" "nvx4@cdc.gov"


end_date=$(date "+%m_%d_%Y_at_%Hh_%Mm")
echo "Run ended at ${end_date}"
#Script exited gracefully (unless something else inside failed)
exit 0
