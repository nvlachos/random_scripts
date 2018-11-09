#!/bin/sh -l

#$ -o act_by_list_barebones_4.out
#$ -e act_by_list_barebones_4.err
#$ -N abl-prokka
#$ -cwd
#$ -q all.q

#Import the config file with shortcuts and settings
. ./config.sh
#Import the module file that loads all necessary mods
. "${mod_changers}/pipeline_mods"

#
# Usage ./act_by_list.sh list_name(currently has to be placed in /scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR folder) Description of list function
#
# script changes depending on what needs to be run through the list
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to act_by_list.sh, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./act_by_list.sh path_to_list_file(single sample ID per line, e.g. B8VHY/1700128 (it must include project id also))"
	echo "Output location varies depending on which tasks are performed but will be found somewhere under ${processed}"
	exit 0
fi

# Loop through and act on each sample name in the passed/provided list
counter=0
echo "c-sstar:c-sstar_plasmid:srst2"
echo "ID:0608-c-sstar:0608-c-sstar_plasmid:0608-srst2:1003-c-sstar:1003-c-sstar_plasmid:1003-srst2" > "${share}/current_DBS_in_samples.txt"
while IFS= read -r var; do
	sample_name=$(echo "${var}" | cut -d'/' -f2 | tr -d '[:space:]')
	project=$(echo "${var}" | cut -d'/' -f1 | tr -d '[:space:]')
	echo "${counter}"
	echo "Looking for ${processed}/${project}/${sample_name}/c-sstar/${sample_name}.ResGANNOT_20180608_srst2.gapped_98_sstar_summary.txt"
	if [[ -f "${processed}/${project}/${sample_name}/c-sstar/${sample_name}.ResGANNOT_20180608_srst2.gapped_98_sstar_summary.txt" ]]; then
		echo "Found-renaming"
		mv "${processed}/${project}/${sample_name}/c-sstar/${sample_name}.ResGANNOT_20180608_srst2.gapped_98_sstar_summary.txt" "${processed}/${project}/${sample_name}/c-sstar/${sample_name}.ResGANNOT_20180608.gapped_98_sstar_summary.txt"
	else
		echo "Not Found"
	fi
	echo "Looking for ${processed}/${project}/${sample_name}/c-sstar/ResGANNOT_20180608_srst2"
	if [[ -d "${processed}/${project}/${sample_name}/c-sstar/ResGANNOT_20180608_srst2" ]]; then
		echo "Found - renaming"
		echo "Looking for ${processed}/${project}/${sample_name}/c-sstar/ResGANNOT_20180608_srst2/${sample}.ResGANNOT_20180608_srst2.gapped_98.sstar_grouped"
		if [[ -f "${processed}/${project}/${sample_name}/c-sstar/ResGANNOT_20180608_srst2/${sample}.ResGANNOT_20180608_srst2.gapped_98.sstar_grouped" ]]; then
			echo "Found - renaming"
			mv "${processed}/${project}/${sample_name}/c-sstar/ResGANNOT_20180608_srst2/${sample}.ResGANNOT_20180608_srst2.gapped_98.sstar_grouped" "${processed}/${project}/${sample_name}/c-sstar/ResGANNOT_20180608_srst2ResGANNOT_20180608.gapped_98.sstar_grouped"
		else
			echo "Not Found"
		fi
		echo "Looking for ${processed}/${project}/${sample_name}/c-sstar/ResGANNOT_20180608_srst2/${sample}.ResGANNOT_20180608_srst2.gapped_98.sstar"
		if [[ -f "${processed}/${project}/${sample_name}/c-sstar/ResGANNOT_20180608_srst2/${sample}.ResGANNOT_20180608_srst2.gapped_98.sstar" ]]; then
			echo "Found - renaming"
			mv "${processed}/${project}/${sample_name}/c-sstar/ResGANNOT_20180608_srst2/${sample}.ResGANNOT_20180608_srst2.gapped_98.sstar" "${processed}/${project}/${sample_name}/c-sstar/ResGANNOT_20180608_srst2ResGANNOT_20180608.gapped_98.sstar"
		else
				echo "Not Found"
		fi
		mv -r ${processed}/${project}/${sample_name}/c-sstar/ResGANNOT_20180608_srst2 ${processed}/${project}/${sample_name}/c-sstar/ResGANNOT_20180608
	else
		echo "Not Found"
	fi
	if [[ -f "${processed}/${project}/${sample_name}/c-sstar_plasmid/${sample_name}.ResGANNOT_20180608_srst2.gapped_40_sstar_summary.txt" ]]; then
		mv "${processed}/${project}/${sample_name}/c-sstar_plasmid/${sample_name}.ResGANNOT_20180608_srst2.gapped_40_sstar_summary.txt" "${processed}/${project}/${sample_name}/c-sstar/${sample_name}.ResGANNOT_20180608.gapped_40_sstar_summary.txt"
	fi
	if [[ -d "${processed}/${project}/${sample_name}/c-sstar_plasmid/ResGANNOT_20180608_srst2" ]]; then
		if [[ -f "${processed}/${project}/${sample_name}/c-sstar_plasmid/ResGANNOT_20180608_srst2/${sample}.ResGANNOT_20180608_srst2.gapped_40.sstar_grouped" ]]; then
			mv "${processed}/${project}/${sample_name}/c-sstar_plasmid/ResGANNOT_20180608_srst2/${sample}.ResGANNOT_20180608_srst2.gapped_40.sstar_grouped" "${processed}/${project}/${sample_name}/c-sstar/ResGANNOT_20180608_srst2ResGANNOT_20180608.gapped_40.sstar_grouped"
		fi
		if [[ -f "${processed}/${project}/${sample_name}/c-sstar_plasmid/ResGANNOT_20180608_srst2/${sample}.ResGANNOT_20180608_srst2.gapped_40.sstar" ]]; then
			mv "${processed}/${project}/${sample_name}/c-sstar_plasmid/ResGANNOT_20180608_srst2/${sample}.ResGANNOT_20180608_srst2.gapped_40.sstar" "${processed}/${project}/${sample_name}/c-sstar/ResGANNOT_20180608_srst2ResGANNOT_20180608.gapped_40.sstar"
		fi
		mv -r ${processed}/${project}/${sample_name}/c-sstar_plasmid/ResGANNOT_20180608_srst2 ${processed}/${project}/${sample_name}/c-sstar/ResGANNOT_20180608
	fi
	if [[ -f "${processed}/${project}/${sample_name}/srst2/${sample_name}_ResGANNOT_20180608_srst2.fasta__fullgenes__ResGANNOT_20180608_srst2__results.txt" ]]; then
		mv "${processed}/${project}/${sample_name}/srst2/${sample_name}_ResGANNOT_20180608_srst2.fasta__fullgenes__ResGANNOT_20180608_srst2__results.txt" "${processed}/${project}/${sample_name}/srst2/${sample_name}__fullgenes__ResGANNOT_20180608_srst2__results.txt"
	fi
	if [[ -f "${processed}/${project}/${sample_name}/srst2/${sample_name}_ResGANNOT_20180608_srst2.fasta__genes__ResGANNOT_20180608_srst2__results.txt" ]]; then
		mv "${processed}/${project}/${sample_name}/srst2/${sample_name}_ResGANNOT_20180608_srst2.fasta__genes__ResGANNOT_20180608_srst2__results.txt" "${processed}/${project}/${sample_name}/srst2/${sample_name}__genes__ResGANNOT_20180608_srst2__results.txt"
	fi
	if [[ -f "${processed}/${project}/${sample_name}/c-sstar/${sample_name}.ResGANNOT_20180608.gapped_98_sstar_summary.txt" ]]; then
		ohsixoheight="X"
	elif [[ ! -s "${processed}/${project}/${sample_name}/Assembly/${sample_name}_scaffolds_trimmed.fasta" ]]; then
		ohsixoheight="A"
	else
		ohsixoheight="M"
	fi
	if [[ -d "${processed}/${project}/${sample_name}/c-sstar_plasmid" ]]; then
		if [[ -f "${processed}/${project}/${sample_name}/c-sstar_plasmid/${sample_name}.ResGANNOT_20180608.gapped_40_sstar_summary.txt" ]]; then
			ohsixoheightp="X"
		elif [[ ! -s "${processed}/${project}/${sample_name}/Assembly/${sample_name}_scaffolds_trimmed.fasta" ]]; then
			ohsixoheightp="A"
		else
			ohsixoheightp="M"
		fi
	else
		ohsixoheightp="-"
	fi
	if [[ -f "${processed}/${project}/${sample_name}/c-sstar/${sample_name}.ResGANNOT_20181003.gapped_98_sstar_summary.txt" ]]; then
		oheightseventeen="X"
	elif [[ ! -s "${processed}/${project}/${sample_name}/Assembly/${sample_name}_scaffolds_trimmed.fasta" ]]; then
		oheightseventeen="A"
	else
		oheightseventeen="M"
	fi
	if [[ -d "${processed}/${project}/${sample_name}/c-sstar_plasmid" ]]; then
		if [[ -f "${processed}/${project}/${sample_name}/c-sstar_plasmid/${sample_name}.ResGANNOT_20181003.gapped_40_sstar_summary.txt" ]]; then
			oheightseventeenp="X"
		elif [[ ! -s "${processed}/${project}/${sample_name}/Assembly/${sample_name}_scaffolds_trimmed.fasta" ]]; then
			oheightseventeenp="A"
		else
			oheightseventeenp="M"
		fi
	else
		oheightseventeenp="-"
	fi
	if [[ -f "${processed}/${project}/${sample_name}/srst2/${sample_name}__genes__ResGANNOT_20180608_srst2__results.txt" ]] || [[ -f "${processed}/${project}/${sample_name}/srst2/${sample_name}__fullgenes__ResGANNOT_20180608_srst2__results.txt" ]]; then
		ohsixoheights="X"
	else
		ohsixoheights="M"
	fi
	if [[ -f "${processed}/${project}/${sample_name}/srst2/${sample_name}__genes__ResGANNOT_20181003_srst2__results.txt" ]] || [[ -f "${processed}/${project}/${sample_name}/srst2/${sample_name}__fullgenes__ResGANNOT_20181003_srst2__results.txt" ]]; then
		oheightseventeens="X"
	else
		oheightseventeens="M"
	fi
	continue
	echo "${counter}:${project}/${sample_name}:	0608	:${ohsixoheight}:${ohsixoheightp}:${ohsixoheights}:	1003	:${oheightseventeen}:${oheightseventeenp}:${oheightseventeens}:"
	echo "${project}/${sample_name}:${ohsixoheight}:${ohsixoheightp}:${ohsixoheights}:${oheightseventeen}:${oheightseventeenp}:${oheightseventeens}:" >> "${share}/current_DBS_in_samples.txt"
	counter=$(( counter + 1 ))
done < "${1}"
echo "All isolates completed"
global_end_time=$(date "+%m-%d-%Y @ %Hh_%Mm_%Ss")
#Script exited gracefully (unless something else inside failed)
printf "%s %s" "Act_by_list.sh has completed ${2} " "${global_end_time}" | mail -s "act_by_list complete" nvx4@cdc.gov
exit 0
