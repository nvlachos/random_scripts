#!/bin/sh -l

#$ -o abl-ARC.out
#$ -e abl-ARC.err
#$ -N abl-AR_checker
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
. ./config.sh
#Import the module file that loads all necessary mods
. "${mod_changers}/pipeline_mods"

#
# Usage ./act_by_list_AR_completion_check.sh path_to_list ResGANNCBI_identifier(YYYYMMDD) path_for_output_file
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./act_by_list_AR_completion_check.sh path_to_list_file ResGANNCBI_Identifier(YYYYMMDD) path_for_output_file"
	exit 0
elif [[ -z "${2}" ]]; then
	echo "ResGANNCBI identifier is empty, please input the DB identifier using YYYYMMDD that it was created, exiting..."
	exit 1
elif [[ -z "${3}" ]]; then
	echo  "No output file input, exiting..."
	exit 1
fi

# Loop through and act on each sample name in the passed/provided list
counter=0
echo "c-sstar:c-sstar_plasmid:srst2"
echo "Identification	20180608-c-sstar	20180608-srst2	${2}-c-sstar	${2}-srst2 plaFlow ${2}-c-sstar-plasFlow	plasmidFinder	plasmidFinder_on_plasFlow" > "${3}"
while IFS= read -r var || [ -n "$var" ]; do
	sample_name=$(echo "${var}" | cut -d'/' -f2 | tr -d '[:space:]')
	project=$(echo "${var}" | cut -d'/' -f1 | tr -d '[:space:]')
#	echo "${counter}"

#	Checking for misnamed files and folders from a short span when the script wasnt not correctly parsing names
	if [[ -d "${processed}/${project}/${sample_name}/c-sstar/ResGANNOT_20180608_srst2_gapped" ]]; then
		echo "Found - renaming"
		mv "${processed}/${project}/${sample_name}/c-sstar/ResGANNOT_20180608_srst2_gapped" "${processed}/${project}/${sample_name}/c-sstar/ResGANNOT_20180608_gapped"
	fi
	if [[ -d "${processed}/${project}/${sample_name}/c-sstar_plasmid/ResGANNOT_20180608_srst2_gapped" ]]; then
		echo "Found - renaming"
		mv "${processed}/${project}/${sample_name}/c-sstar_plasmid/ResGANNOT_20180608_srst2_gapped" "${processed}/${project}/${sample_name}/c-sstar/ResGANNOT_20180608_gapped"
	fi
	if [[ -f "${processed}/${project}/${sample_name}/c-sstar/${sample_name}.ResGANNOT_20180608_srst2.gapped_98_sstar_summary.txt" ]]; then
		echo "Found - renaming"
		mv "${processed}/${project}/${sample_name}/c-sstar/${sample_name}.ResGANNOT_20180608_srst2.gapped_98_sstar_summary.txt" "${processed}/${project}/${sample_name}/c-sstar/${sample_name}.ResGANNOT_20180608.gapped_98_sstar_summary.txt"
	fi
	if [[ -f "${processed}/${project}/${sample_name}/c-sstar_plasmid/${sample_name}.ResGANNOT_20180608_srst2.gapped_40_sstar_summary.txt" ]]; then
		echo "Found - renaming"
		mv "${processed}/${project}/${sample_name}/c-sstar_plasmid/${sample_name}.ResGANNOT_20180608_srst2.gapped_40_sstar_summary.txt" "${processed}/${project}/${sample_name}/c-sstar/${sample_name}.ResGANNOT_20180608.gapped_40_sstar_summary.txt"
	fi
	if [[ -f "${processed}/${project}/${sample_name}/srst2/${sample_name}_ResGANNOT_20180608_srst2.fasta__fullgenes__ResGANNOT_20180608_srst2__results.txt" ]]; then
		echo "Found - renaming"
		rm "${processed}/${project}/${sample_name}/srst2/${sample_name}_ResGANNOT_20180608_srst2.fasta__fullgenes__ResGANNOT_20180608_srst2__results.txt" "${processed}/${project}/${sample_name}/srst2/${sample_name}__fullgenes__ResGANNOT_20180608_srst2__results.txt"
	fi
	if [[ -f "${processed}/${project}/${sample_name}/srst2/${sample_name}_ResGANNOT_20180608_srst2.fasta__genes__ResGANNOT_20180608_srst2__results.txt" ]]; then
		echo "Found - renaming"
		mv "${processed}/${project}/${sample_name}/srst2/${sample_name}_ResGANNOT_20180608_srst2.fasta__genes__ResGANNOT_20180608_srst2__results.txt" "${processed}/${project}/${sample_name}/srst2/${sample_name}__genes__ResGANNOT_20180608_srst2__results.txt"
	fi
	if [[ -f "${processed}/${project}/${sample_name}/srst2/${sample_name}_ResGANNOT_20180608__fullgenes__ResGANNOT_20180608_srst2__results.txt" ]]; then
		echo "Found - renaming"
		rm "${processed}/${project}/${sample_name}/srst2/${sample_name}_ResGANNOT_20180608__fullgenes__ResGANNOT_20180608_srst2__results.txt" "${processed}/${project}/${sample_name}/srst2/${sample_name}__fullgenes__ResGANNOT_20180608_srst2__results.txt"
	fi
	if [[ -f "${processed}/${project}/${sample_name}/srst2/${sample_name}_ResGANNOT_20180608__genes__ResGANNOT_20180608_srst2__results.txt" ]]; then
		echo "Found - renaming"
		mv "${processed}/${project}/${sample_name}/srst2/${sample_name}_ResGANNOT_20180608__genes__ResGANNOT_20180608_srst2__results.txt" "${processed}/${project}/${sample_name}/srst2/${sample_name}__genes__ResGANNOT_20180608_srst2__results.txt"
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
	if [[ -s "${processed}/${project}/${sample_name}/Assembly/${sample_name}_scaffolds_trimmed.fasta" ]]; then
		if [[ -f "${processed}/${project}/${sample_name}/c-sstar/${sample_name}.ResGANNOT_20180608.gapped_98_sstar_summary.txt" ]]; then
			header=$(head -n1 "${processed}/${project}/${sample_name}/c-sstar/${sample_name}.ResGANNOT_20180608.gapped_98_sstar_summary.txt")
			if [[ "${header}" = "No anti-microbial genes were found"* ]]; then
				ohsixoheight="No_chromo_AR"
			else
				ohsixoheight="AR_Found"
			fi
		else
			ohsixoheight="NO_CSSTAR_file(HAS_ASSEMBLY)"
		fi
	else
		ohsixoheight="NO_ASSEMBLY_TO_RUN_CSSTAR_ON"
	fi

	if [[ -s "${processed}/${project}/${sample_name}/Assembly/${sample_name}_scaffolds_trimmed.fasta" ]]; then
		if [[ -f "${processed}/${project}/${sample_name}/c-sstar/${sample_name}.ResGANNCBI_${2}.gapped_98_sstar_summary.txt" ]]; then
			header=$(head -n1 "${processed}/${project}/${sample_name}/c-sstar/${sample_name}.ResGANNCBI_${2}.gapped_98_sstar_summary.txt")
			if [[ "${header}" = "No anti-microbial genes were found"* ]]; then
				input_DB_csstar="No_chromo_AR"
			else
				input_DB_csstar="AR_Found"
			fi
		else
			input_DB_csstar="NO_CSSTAR_file(HAS_ASSEMBLY)"
		fi
	else
		input_DB_csstar="No_ASSEMBLY_TO_RUN_CSSTAR_ON"
	fi

	# Brief check if srst2 files exist, dont really have time to check for content at the moment
	if [[ -f "${processed}/${project}/${sample_name}/srst2/${sample_name}__genes__ResGANNOT_20180608_srst2__results.txt" ]] || [[ -f "${processed}/${project}/${sample_name}/srst2/${sample_name}__fullgenes__ResGANNOT_20180608_srst2__results.txt" ]]; then
		ohsixoheights="File Found"
	else
		ohsixoheights="File Missing"
	fi
	if [[ -f "${processed}/${project}/${sample_name}/srst2/${sample_name}__genes__ResGANNCBI_${2}_srst2__results.txt" ]] || [[ -f "${processed}/${project}/${sample_name}/srst2/${sample_name}__fullgenes__ResGANNCBI_${2}_srst2__results.txt" ]]; then
		input_DB_srst2="File Found"
	else
		input_DB_srst2="File Missing"
	fi

	if [[ -f "${processed}/${project}/${sample_name}/${sample_name}.tax" ]]; then
		# Parse tax file
		family=""
		while IFS= read -r line  || [ -n "$line" ]; do
			# Grab first letter of line (indicating taxonomic level)
			first=${line:0:1}
			# Assign taxonomic level value from 4th value in line (1st-classification level, 2nd-% by kraken, 3rd-true % of total reads, 4th-identifier)
			if [ "${first}" = "F" ];
			then
				family=$(echo "${line}" | awk -F ' ' '{print $2}')
			fi
		done < "${processed}/${project}/${sample_name}/${sample_name}.tax"
		if [[ "${family}" == "Enterobacteriaceae" ]]; then
			if [[ -s "${processed}/${project}/${sample_name}/plasFlow/Unicycler_assemblies/${sample_name}_uni_assembly/${sample_name}_plasmid_assembly_trimmed.fasta" ]]; then
				plasFlow="Completed"
				if [[ -f "${processed}/${project}/${sample_name}/c-sstar_plasFlow/${sample_name}.ResGANNCBI_${2}.gapped_40_sstar_summary.txt" ]]; then
					header=$(head -n1 "${processed}/${project}/${sample_name}/c-sstar_plasFlow/${sample_name}.ResGANNCBI_${2}.gapped_40_sstar_summary.txt")
					if [[ "${header}" = "No anti-microbial genes were found"* ]]; then
						cplas="No_chromo_AR"
					else
						csplas="AR_Found"
					fi
					## Check run_plasmidFinder
				if [[ -f "${processed}/${project}/${sample_name}/plasmid_on_plasFlow/${sample_name}_results_table_summary.txt" ]]; then
					pfin_plas="Found"
				else
					pfin_plas="NOT_found"
				fi
				else
					cplas="NO_plasFlow_CSSTAR_file(HAS_plasFlow_ASSEMBLY)"
				fi
			fi
		elif [[ "${family}" == "" ]]; then
			plasFlow="No_family_in_TAX_file"
			csplas="No_plasFlow"
			pfin_plas="No_plasFlow"
		else
			plasFlow="Not_ENTEROBACTERIACEAE_family"
			cplas="No_plasFlow"
			pfin_plas="No_plasFlow"
		fi
	else
		plasFlow="No_TAX_file"
		csplas="No_plaslow"
		pfin_plas="No_plasFlow"
	fi

	## Check run_plasmidFinder
if [[ -f "${processed}/${project}/${sample_name}/plasmid/${sample_name}_results_table_summary.txt" ]]; then
	pfin="Found"
else
	pfin="NOT_found"
fi


	echo "${counter}:${project}/${sample_name}:	20180608	:${ohsixoheight}:${ohsixoheightp}:${ohsixoheights}:	${2}	:${input_DB_csstar}:${input_DB_srst2}:${plasFlow}:${cplas}:${pfin}:${pfin_plas}"
	echo "${project}/${sample_name}	${ohsixoheight}	${ohsixoheights}	${input_DB_csstar}	${input_DB_srst2}	${plasFlow}	${cplas}	${pfin}	${pfin_plas}" >> "${3}"
	counter=$(( counter + 1 ))
done < "${1}"
echo "All isolates completed"
global_end_time=$(date "+%m-%d-%Y @ %Hh_%Mm_%Ss")
#Script exited gracefully (unless something else inside failed)
printf "%s %s" "act_by_list_AR_completion_check.sh has completed check for AR completion " "${global_end_time}" | mail -s "act_by_list complete" nvx4@cdc.gov
exit 0
