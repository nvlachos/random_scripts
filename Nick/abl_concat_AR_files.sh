#!/bin/sh -l

#$ -o abl-ARC.out
#$ -e abl-ARC.err
#$ -N abl-plasCATTER
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
. ./config.sh
#Import the module file that loads all necessary mods
. "${mod_changers}/pipeline_mods"

#
# Usage ./act_by_list_concat_csstar_plasFlow_files.sh path_to_list path_for_output_file %ID_to_find_correct_csstar_file[80|95|98|99|100]
#

number='^[0-9]+$'

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./act_by_list_concat_csstar_plasFlow_files.sh path_to_list_file DB_Identifier(ResGANNCBI_date, for example) path_for_output_file %ID_to_find_correct_csstar_file[80|95|98|99|100]"
	exit 0
elif [[ -z "${2}" ]]; then     #! [[ ${2} =~ $number ]] ||
	echo "${2} is not a number or is empty. Please input correct date for ResGANNCBI DB...exiting"
	exit 2
elif [[ ! -f "${3}" ]]; then
	echo "${3} does not exist"
	outdir=$(dirname "${3}")
	if [[ ! -d "${outdir}" ]]; then
		echo "${outdir} does not exist"
		mkdir -p "${outdir}"
	fi
	#>${3}_plasFLOW
	#>${3}_csstar
	#>${3}_srst2
fi

# Checks that value given for % Identity is one of the presets for csstar
if [[ "${4}" != 80 ]] && [[ "${4}" != 95 ]] && [[ "${4}" != 98 ]] && [[ "${4}" != 99 ]] && [[ "${4}" != 100 ]]; then
	echo "Identity is not one of the presets for csstar and therefore will fail, defaulting to 98..."
	sim="h"
	simnum=98
else
	if [ "${4}" == 98 ]; then
		sim="h"
	elif [ "${4}" == 80 ]; then
		sim="l"
	elif [ "${4}" == 99 ]; then
		sim="u"
	elif [ "${4}" == 95 ]; then
		sim="m"
	elif [ "${4}" == 100 ]; then
		sim="p"
	elif [ "${4}" == 40 ]; then
		sim="o"
	fi
	simnum=${4}
fi

# Loop through and act on each sample name in the passed/provided list
counter=0
while IFS= read -r var || [ -n "$var" ]; do
	sample_name=$(echo "${var}" | cut -d'/' -f2 | tr -d '[:space:]')
	project=$(echo "${var}" | cut -d'/' -f1 | tr -d '[:space:]')
	#echo "${counter} - ${processed}/${project}/${sample_name}/MLST/"
	#if [[ "${counter}" -gt 25 ]]; then
	#	break
	#fi
	if [[ -f "${processed}/${project}/${sample_name}/c-sstar_plasFlow/${sample_name}.${2}.gapped_40_sstar_summary.txt" ]]; then
		:
		cat "${processed}/${project}/${sample_name}/c-sstar_plasFlow/${sample_name}.${2}.gapped_40_sstar_summary.txt" >> ${3}_plasFlow.tsv
	else
		echo "${project}/${sample_name} does not have c-sstar_plasFlow 40 file"
	fi
	if [[ -f "${processed}/${project}/${sample_name}/c-sstar/${sample_name}.${2}.gapped_${simnum}_sstar_summary.txt" ]]; then
		if [[ "${counter}" -eq 0 ]]; then
			echo "Source_DB	AR_Family	sstar_on_single_alternate_DB	Gene_Status	Gene_Allele	Contig#	%ID	Source_Length	Target_length	%_Length" >> ${3}_csstar.tsv
		fi
		tail -n +1 "${processed}/${project}/${sample_name}/c-sstar/${sample_name}.${2}.gapped_${simnum}_sstar_summary.txt" >> ${3}_csstar.tsv
	else
		echo "${project}/${sample_name} does not have c-sstar ${simnum} file"
	fi
	if [[ -f "${processed}/${project}/${sample_name}/srst2/${sample_name}__fullgenes__${2}_srst2__results.txt" ]]; then
		if [[ "${counter}" -eq 0 ]]; then
			tail -n +1 "${processed}/${project}/${sample_name}/srst2/${sample_name}__fullgenes__${2}_srst2__results.txt" >> ${3}_srst2.tsv
		else
			tail -n +2 "${processed}/${project}/${sample_name}/srst2/${sample_name}__fullgenes__${2}_srst2__results.txt" >> ${3}_srst2.tsv
		fi
	else
		echo "${project}/${sample_name} does not have srst2 file"
	fi
	counter=$(( counter + 1 ))
done < "${1}"

echo "All isolates completed"
global_end_time=$(date "+%m-%d-%Y @ %Hh_%Mm_%Ss")
#Script exited gracefully (unless something else inside failed)
printf "%s %s" "act_by_list_concat_csstar_plasFlow_files.sh has completed check of snd MLSTs" "${global_end_time}" | mail -s "act_by_list complete" nvx4@cdc.gov
exit 0
