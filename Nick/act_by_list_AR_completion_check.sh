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
# Usage ./act_by_list_AR_completion_check.sh path_to_list ResGANNOT_identifier(YYYYMMDD)
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to act_by_list.sh, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./act_by_list_AR_completion_check.sh path_to_list_file ResGANNOT_Identifier(YYYYMMDD)"
	exit 0
elif [[ -z "${2}" ]]; then
	echo "ResGANNOT identifier is empty, please input the DB identifier using YYYYMMDD that it was created, exiting..."
	exit 1
fi

# Loop through and act on each sample name in the passed/provided list
counter=0
echo "c-sstar:c-sstar_plasmid:srst2"
echo "Identification:0608-c-sstar:0608-c-sstar_plasmid:0608-srst2:1204-c-sstar:1204-c-sstar_plasmid:1204-srst2" > "${share}/current_DBS_in_samples.txt"
while IFS= read -r var; do
	sample_name=$(echo "${var}" | cut -d'/' -f2 | tr -d '[:space:]')
	project=$(echo "${var}" | cut -d'/' -f1 | tr -d '[:space:]')
#	echo "${counter}"

#	Checking for misnamed files and folders from a short span when the script wasnt not correctly parsing names
	if [[ -d "${processed}/${project}/${sample_name}/c-sstar/ResGANNOT_20180608_srst2_gapped" ]]; then
		echo "Found - renaming"
		mv "${processed}/${project}/${sample_name}/c-sstar/ResGANNOT_20180608_srst2_gapped" "${processed}/${project}/${sample_name}/c-sstar/ResGANNOT_20180608_gapped"
	fi
	if [[ -d "${processed}/${project}/${sample_name}/c-sstar_plasmid/ResGANNOT_20180608_srst2_gapped" ]]; then
		echo "Found-renaming"
		mv "${processed}/${project}/${sample_name}/c-sstar_plasmid/ResGANNOT_20180608_srst2_gapped" "${processed}/${project}/${sample_name}/c-sstar/ResGANNOT_20180608_gapped"
	fi
	if [[ -f "${processed}/${project}/${sample_name}/c-sstar/${sample_name}.ResGANNOT_20180608_srst2.gapped_98_sstar_summary.txt" ]]; then
		echo "Found-renaming"
		mv "${processed}/${project}/${sample_name}/c-sstar/${sample_name}.ResGANNOT_20180608_srst2.gapped_98_sstar_summary.txt" "${processed}/${project}/${sample_name}/c-sstar/${sample_name}.ResGANNOT_20180608.gapped_98_sstar_summary.txt"
	fi
	if [[ -f "${processed}/${project}/${sample_name}/c-sstar_plasmid/${sample_name}.ResGANNOT_20180608_srst2.gapped_40_sstar_summary.txt" ]]; then
		echo "Found-renaming"
		mv "${processed}/${project}/${sample_name}/c-sstar_plasmid/${sample_name}.ResGANNOT_20180608_srst2.gapped_40_sstar_summary.txt" "${processed}/${project}/${sample_name}/c-sstar/${sample_name}.ResGANNOT_20180608.gapped_40_sstar_summary.txt"
	fi
	if [[ -f "${processed}/${project}/${sample_name}/srst2/${sample_name}_ResGANNOT_20180608_srst2.fasta__fullgenes__ResGANNOT_20180608_srst2__results.txt" ]]; then
		echo "Found-renaming"
		rm "${processed}/${project}/${sample_name}/srst2/${sample_name}_ResGANNOT_20180608_srst2.fasta__fullgenes__ResGANNOT_20180608_srst2__results.txt" "${processed}/${project}/${sample_name}/srst2/${sample_name}__fullgenes__ResGANNOT_20180608_srst2__results.txt"
	fi
	if [[ -f "${processed}/${project}/${sample_name}/srst2/${sample_name}_ResGANNOT_20180608_srst2.fasta__genes__ResGANNOT_20180608_srst2__results.txt" ]]; then
		echo "Found-renaming"
		mv "${processed}/${project}/${sample_name}/srst2/${sample_name}_ResGANNOT_20180608_srst2.fasta__genes__ResGANNOT_20180608_srst2__results.txt" "${processed}/${project}/${sample_name}/srst2/${sample_name}__genes__ResGANNOT_20180608_srst2__results.txt"
	fi
	if [[ -f "${processed}/${project}/${sample_name}/c-sstar/ResGANNOT_20180608_gapped/${sample_name}.ResGANNOT_20180608_srst2.gapped_98.sstar_grouped" ]]; then
		echo "Found - renaming"
		mv "${processed}/${project}/${sample_name}/c-sstar/ResGANNOT_20180608_gapped/${sample_name}.ResGANNOT_20180608_srst2.gapped_98.sstar_grouped" "${processed}/${project}/${sample_name}/c-sstar/ResGANNOT_20180608_gapped/${sample_name}.ResGANNOT_20180608.gapped_98.sstar_grouped"
	fi
	if [[ -f "${processed}/${project}/${sample_name}/c-sstar/ResGANNOT_20180608_gapped/${sample_name}.ResGANNOT_20180608_srst2.gapped_98.sstar" ]]; then
		echo "Found - renaming"
		mv "${processed}/${project}/${sample_name}/c-sstar/ResGANNOT_20180608_gapped/${sample_name}.ResGANNOT_20180608_srst2.gapped_98.sstar" "${processed}/${project}/${sample_name}/c-sstar/ResGANNOT_20180608_gapped/${sample_name}.ResGANNOT_20180608.gapped_98.sstar"
	fi
	if [[ -f "${processed}/${project}/${sample_name}/c-sstar_plasmid/ResGANNOT_20180608_gapped/${sample_name}.ResGANNOT_20180608_srst2.gapped_40.sstar_grouped" ]]; then
		echo "Found - renaming"
		mv "${processed}/${project}/${sample_name}/c-sstar_plasmid/ResGANNOT_20180608_gapped/${sample_name}.ResGANNOT_20180608_srst2.gapped_40.sstar_grouped" "${processed}/${project}/${sample_name}/c-sstar/ResGANNOT_20180608_gapped/${sample_name}.ResGANNOT_20180608.gapped_40.sstar_grouped"
	fi
	if [[ -f "${processed}/${project}/${sample_name}/c-sstar_plasmid/ResGANNOT_20180608_gapped/${sample_name}.ResGANNOT_20180608_srst2.gapped_40.sstar" ]]; then
		echo "Found - renaming"
		mv "${processed}/${project}/${sample_name}/c-sstar_plasmid/ResGANNOT_20180608_gapped/${sample_name}.ResGANNOT_20180608_srst2.gapped_40.sstar" "${processed}/${project}/${sample_name}/c-sstar/ResGANNOT_20180608_gapped/${sample_name}.ResGANNOT_20180608.gapped_40.sstar"
	fi
	if [[ -f "${processed}/${project}/${sample_name}/srst2/${sample_name}_ResGANNOT_20180608_srst2.fasta__fullgenes__ResGANNOT_20180608_srst2__results.txt" ]]; then
		echo "Found - renaming"
		mv "${processed}/${project}/${sample_name}/srst2/${sample_name}_ResGANNOT_20180608_srst2.fasta__fullgenes__ResGANNOT_20180608_srst2__results.txt" "${processed}/${project}/${sample_name}/srst2/${sample_name}__fullgenes__ResGANNOT_20180608_srst2__results.txt"
	fi
	if [[ -f "${processed}/${project}/${sample_name}/srst2/${sample_name}_ResGANNOT_20180608_srst2.fasta__genes__ResGANNOT_20180608_srst2__results.txt" ]]; then
		echo "Found - renaming"
		mv "${processed}/${project}/${sample_name}/srst2/${sample_name}_ResGANNOT_20180608_srst2.fasta__genes__ResGANNOT_20180608_srst2__results.txt" "${processed}/${project}/${sample_name}/srst2/${sample_name}__genes__ResGANNOT_20180608_srst2__results.txt"
	fi


# Check to see if all samples have the normal 0608 AR csstar output
	if [[ -f "${processed}/${project}/${sample_name}/c-sstar/${sample_name}.ResGANNOT_20180608.gapped_98_sstar_summary.txt" ]]; then
		header=$(head -n1 "${processed}/${project}/${sample_name}/c-sstar/${sample_name}.ResGANNOT_20180608.gapped_98_sstar_summary.txt")
		if [[ "${header}" = "No anti-microbial genes were found"* ]]; then
			ohsixoheight="None"
		else
			ohsixoheight="AR_Found"
		fi
	elif [[ ! -s "${processed}/${project}/${sample_name}/Assembly/${sample_name}_scaffolds_trimmed.fasta" ]]; then
		ohsixoheight="0_size"
	else
		ohsixoheight="File_Missing"
	fi
	if [[ -d "${processed}/${project}/${sample_name}/c-sstar_plasmid" ]]; then
		if [[ -f "${processed}/${project}/${sample_name}/c-sstar_plasmid/${sample_name}.ResGANNOT_20180608.gapped_40_sstar_summary.txt" ]]; then
			header=$(head -n1 "${processed}/${project}/${sample_name}/c-sstar_plasmid/${sample_name}.ResGANNOT_20180608.gapped_40_sstar_summary.txt")
			if [[ "${header}" = "No anti-microbial genes were found"* ]]; then
				ohsixoheightp="None"
			else
				ohsixoheightp="AR_Found"
			fi
		elif [[ ! -s "${processed}/${project}/${sample_name}/Assembly/${sample_name}_scaffolds_trimmed.fasta" ]]; then
			ohsixoheightp="0_size"
		else
			ohsixoheightp="File_Missing"
		fi
	else
		ohsixoheightp="No_plasmid_Assembly"
	fi
	if [[ -f "${processed}/${project}/${sample_name}/c-sstar/${sample_name}.ResGANNOT_${2}.gapped_98_sstar_summary.txt" ]]; then
		header=$(head -n1 "${processed}/${project}/${sample_name}/c-sstar/${sample_name}.ResGANNOT_${2}.gapped_98_sstar_summary.txt")
		if [[ "${header}" = "No anti-microbial genes were found"* ]]; then
			input_DB_csstar="None"
		else
			input_DB_csstar="AR_Found"
		fi
	elif [[ ! -s "${processed}/${project}/${sample_name}/Assembly/${sample_name}_scaffolds_trimmed.fasta" ]]; then
		input_DB_csstar="0_size"
	else
		input_DB_csstar="File_Missing"
	fi
	if [[ -d "${processed}/${project}/${sample_name}/c-sstar_plasmid" ]]; then
		if [[ -f "${processed}/${project}/${sample_name}/c-sstar_plasmid/${sample_name}.ResGANNOT_${2}.gapped_40_sstar_summary.txt" ]]; then
			header=$(head -n1 "${processed}/${project}/${sample_name}/c-sstar_plasmid/${sample_name}.ResGANNOT_${2}.gapped_40_sstar_summary.txt")
			if [[ "${header}" = "No anti-microbial genes were found"* ]]; then
				input_DB_csstar_plasmid="None"
			else
				input_DB_csstar_plasmid="AR_Found"
			fi
		elif [[ ! -s "${processed}/${project}/${sample_name}/Assembly/${sample_name}_scaffolds_trimmed.fasta" ]]; then
			input_DB_csstar_plasmid="0_size"
		else
			input_DB_csstar_plasmid="File_Missing"
		fi
	else
		input_DB_csstar_plasmid="No_plasmid_Assembly"
	fi
	# Brief check if srst2 files exist, dont really have time to check for content at the moment
	if [[ -f "${processed}/${project}/${sample_name}/srst2/${sample_name}__genes__ResGANNOT_20180608_srst2__results.txt" ]] || [[ -f "${processed}/${project}/${sample_name}/srst2/${sample_name}__fullgenes__ResGANNOT_20180608_srst2__results.txt" ]]; then
		input_DB_srst2="X"
	else
		input_DB_srst2="M"
	fi
	if [[ -f "${processed}/${project}/${sample_name}/srst2/${sample_name}__genes__ResGANNOT_${2}_srst2__results.txt" ]] || [[ -f "${processed}/${project}/${sample_name}/srst2/${sample_name}__fullgenes__ResGANNOT_${2}_srst2__results.txt" ]]; then
		input_DB_srst2="X"
	else
		input_DB_srst2="M"
	fi
	echo "${counter}:${project}/${sample_name}:	20180608	:${ohsixoheight}:${ohsixoheightp}:${ohsixoheights}:	${2}	:${input_DB_csstar}:${input_DB_csstar_plasmid}:${input_DB_srst2}:"
	echo "${project}/${sample_name}:${ohsixoheight}:${ohsixoheightp}:${ohsixoheights}:${input_DB_csstar}:${input_DB_csstar_plasmid}:${input_DB_srst2}:" >> "${share}/current_DBS_in_samples.txt"
	counter=$(( counter + 1 ))
done < "${1}"
echo "All isolates completed"
global_end_time=$(date "+%m-%d-%Y @ %Hh_%Mm_%Ss")
#Script exited gracefully (unless something else inside failed)
printf "%s %s" "Act_by_list.sh has completed check for AR completion " "${global_end_time}" | mail -s "act_by_list complete" nvx4@cdc.gov
exit 0
