#!/bin/sh -l

#$ -o abl-ARC.out
#$ -e abl-ARC.err
#$ -N abl-AR_checker
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
. ./config.sh
#Import the module file that loads all necessary mods
#. "${mod_changers}/pipeline_mods"

#
# Usage ./act_by_list_AR_completion_check.sh path_to_list ResGANNOT_identifier(YYYYMMDD) path_for_output_file
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./act_by_list_AR_completion_check.sh path_to_list_file ResGANNOT_Identifier(YYYYMMDD) path_for_output_file"
	exit 0
elif [[ -z "${2}" ]]; then
	echo "ResGANNOT identifier is empty, please input the DB identifier using YYYYMMDD that it was created, exiting..."
	exit 1
elif [[ -z "${3}" ]]; then
	echo  "No output file input, exiting..."
	exit 1
fi

# Loop through and act on each sample name in the passed/provided list
counter=0
echo "c-sstar:c-sstar_plasmid:srst2"
echo "Identification	Raw_Read_1	Zipped_Read_1	Raw_Read_2	Zipped_Read_2	Trimmed_R1	Zipped_Trimmed_R1	Trimmed_R2	Zipped_Trimmed_R2	Raw_QC_Counts	Trimmed_QC_Counts	Kraker_reads	Gottcha	Assembly	Assembly_Stats	Kraker_Assembly	BUSCO	PROKKA	family	Genus	Species	16s	mash_ANI	ANI_genus	ANI_All	MLST	MLST_SRST2	plasmidFinder	csstar-20180608	srst2AR-20180608	csstar-ResGANNOT_${2}	srst2-ResGANNOT_${2}	plasFlow	plasFlow_Assembly_Stats	csstar_plasFlow-ResGANNOT_${2}	plasmidFinder_on_plasFlow" > "${3}"
#echo "Identification	20180608-c-sstar	20180608-srst2	${2}-c-sstar	${2}-srst2 plaFlow ${2}-c-sstar-plasFlow	plasmidFinder	plasmidFinder_on_plasFlow" > "${3}"
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
		if [[ -f "${processed}/${project}/${sample_name}/c-sstar/${sample_name}.ResGANNOT_${2}.gapped_98_sstar_summary.txt" ]]; then
			header=$(head -n1 "${processed}/${project}/${sample_name}/c-sstar/${sample_name}.ResGANNOT_${2}.gapped_98_sstar_summary.txt")
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
	if [[ -f "${processed}/${project}/${sample_name}/srst2/${sample_name}__genes__ResGANNOT_${2}_srst2__results.txt" ]] || [[ -f "${processed}/${project}/${sample_name}/srst2/${sample_name}__fullgenes__ResGANNOT_${2}_srst2__results.txt" ]]; then
		input_DB_srst2="File Found"
	else
		input_DB_srst2="File Missing"
	fi

	## Check run_plasmidFinder
	if [[ -f "${processed}/${project}/${sample_name}/plasmid/${sample_name}_results_table_summary.txt" ]]; then
		pFin="Found"
	else
		pFin="NOT_found"
	fi

	if [[ -s "${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R1_001.fastq" ]]; then
		FQR1="Found"
	else
		FQR1="NOT FOUND"
	fi

	if [[ -s "${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R1_001.fastq.gz" ]]; then
		FQZR1="Found"
	else
		FQZR1="NOT FOUND"
	fi

	if [[ -s "${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R2_001.fastq" ]]; then
		FQR2="Found"
	else
		FQR2="NOT FOUND"
	fi

	if [[ -s "${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R2_001.fastq.gz" ]]; then
		FQZR2="Found"
	else
		FQZR2="NOT FOUND"
	fi

	if [[ -s "${processed}/${project}/${sample_name}/trimmed/${sample_name}_R1_001.paired.fq" ]]; then
		FQTR1="Found"
	else
		FQTR1="NOT FOUND"
	fi

	if [[ -s "${processed}/${project}/${sample_name}/trimmed/${sample_name}_R1_001.paired.fq.gz" ]]; then
		FQTZR1="Found"
	else
		FQTZR1="NOT FOUND"
	fi

	if [[ -s "${processed}/${project}/${sample_name}/trimmed/${sample_name}_R2_001.paired.fq" ]]; then
		FQTR2="Found"
	else
		FQTR2="NOT FOUND"
	fi

	if [[ -s "${processed}/${project}/${sample_name}/trimmed/${sample_name}_R2_001.paired.fq.gz" ]]; then
		FQTZR2="Found"
	else
		FQTZR2="NOT FOUND"
	fi

	if [[ -s "${processed}/${project}/${sample_name}/preQCcounts/${sample_name}_counts.txt" ]]; then
		preQCr="Found"
	else
		preQCr="NOT FOUND"
	fi

	if [[ -s "${processed}/${project}/${sample_name}/preQCcounts/${sample_name}_trimmed_counts.txt" ]]; then
		preQCt="Found"
	else
		preQCt="NOT FOUND"
	fi

	if [[ -s "${processed}/${project}/${sample_name}/kraken/preAssembly/${sample_name}_kraken_summary_paired.txt" ]]; then
		krakr="Found"
	else
		krakr="NOT_FOUND"
	fi

	if [[ -s "${processed}/${project}/${sample_name}/kraken/postAssembly/${sample_name}_kraken_summary_assembled_BP_data.txt" ]]; then
		kraka="Found"
	else
		kraka="NOT_FOUND"
	fi

	if [[ -s "${processed}/${project}/${sample_name}/gottcha/${sample_name}_gottcha_species_summary.txt" ]]; then
		gott="Found"
	else
		gott="NOT_FOUND"
	fi

	if [[ -s "${processed}/${project}/${sample_name}/Assembly/${sample_name}_scaffolds_trimmed.fasta" ]]; then
		Assembly="Found"
	else
		Assembly="NOT_FOUND"
	fi

	if [[ -s "${processed}/${project}/${sample_name}/Assembly_Stats/${sample_name}_report.tsv" ]]; then
		Assembly_stats="Found"
	else
		Assembly_stats="NOT_FOUND"
	fi

	if [[ -s "${processed}/${project}/${sample_name}/BUSCO/short_summary_${sample_name}.txt" ]]; then
		busco="Found"
	else
		busco="NOT_FOUND"
	fi

	if [[ -s "${processed}/${project}/${sample_name}/${sample_name}.tax" ]]; then
		tax="Found"
	else
		tax="NOT_FOUND"
	fi

	if [[ -s "${processed}/${project}/${sample_name}/16s/${sample_name}_16s_blast_id.txt" ]]; then
		best16s_result=$(head -n1 "${processed}/${project}/${sample_name}/16s/${sample_name}_16s_blast_id.txt" | cut -d'	' -f3)
		largest16s_result=$(tail -n1 "${processed}/${project}/${sample_name}/16s/${sample_name}_16s_blast_id.txt" | cut -d'	' -f3)
		if [[ "${best16s_result}" != "" ]]; then
			best16s="Found"
		else
			best16s="NOT_FOUND"
		fi
		if [[ "${largest16s_result}" != "" ]]; then
			largest16s="Found"
		else
			largest16s="NOT_FOUND"
		fi
		sixteenS="B:${best16s}|L:${largest16s}"
	else
		sixteenS="NOT_FOUND"
	fi

	if [[ -s "${processed}/${project}/${sample_name}/prokka/${sample_name}_PROKKA.gbf" ]] || [[ -s "${processed}/${project}/${sample_name}/prokka/${sample_name}_PROKKA.gbk" ]]; then
		prokka="Found"
	else
		prokka="NOT_FOUND"
	fi

	if [[ -s "${processed}/${project}/${sample_name}/ANI/best_ANI_hits_ordered(${sample_name}_vs_All).txt" ]]; then
		aniAll="Found"
	else
		aniAll="NOT_FOUND"
	fi

	if [[ -s "${processed}/${project}/${sample_name}/MLST/${sample_name}.mlst" ]]; then
		mlst=$(head -n1 "${processed}/${project}/${sample_name}/MLST/${sample_name}.mlst" | cut -d'	' -f3)
	else
		mlst="NOT_FOUND"
	fi




	if [[ "${tax}" == "Found" ]]; then
		# Parse tax file
		family="Unknown"
		genus="Unknown"
		species="Unknown"
		while IFS= read -r line  || [ -n "$line" ]; do
			# Grab first letter of line (indicating taxonomic level)
			first=${line:0:1}
			# Assign taxonomic level value from 4th value in line (1st-classification level, 2nd-% by kraken, 3rd-true % of total reads, 4th-identifier)
			if [ "${first}" = "F" ];
			then
				family=$(echo "${line}" | awk -F ' ' '{print $2}')
			fi
			if [ "${first}" = "G" ];
			then
				genus=$(echo "${line}" | awk -F ' ' '{print $2}')
			fi
			if [ "${first}" = "s" ];
			then
				species=$(echo "${line}" | awk -F ' ' '{print $2}')
			fi
		done < "${processed}/${project}/${sample_name}/${sample_name}.tax"
		if [[ "${family}" == "Enterobacteriaceae" ]]; then
			if [[ -s "${processed}/${project}/${sample_name}/plasFlow/Unicycler_assemblies/${sample_name}_uni_assembly/${sample_name}_plasmid_assembly_trimmed.fasta" ]]; then
				plasFlow="Completed"
				if [[ -f "${processed}/${project}/${sample_name}/c-sstar_plasFlow/${sample_name}.ResGANNOT_${2}.gapped_40_sstar_summary.txt" ]]; then
					header=$(head -n1 "${processed}/${project}/${sample_name}/c-sstar_plasFlow/${sample_name}.ResGANNOT_${2}.gapped_40_sstar_summary.txt")
					if [[ "${header}" = "No anti-microbial genes were found"* ]]; then
						cplas="No_chromo_AR"
					else
						cplas="AR_Found"
					fi
					## Check run_plasmidFinder
				if [[ -f "${processed}/${project}/${sample_name}/plasmid_on_plasFlow/${sample_name}_results_table_summary.txt" ]]; then
					pFin_plas="Found"
				else
					pFin_plas="NOT_FOUND"
				fi
				else
					cplas="NO_plasFlow_CSSTAR_file(HAS_plasFlow_ASSEMBLY)"
				fi
				if [[ -s "${processed}/${project}/${sample_name}/Assembly_Stats_plasFlow/${sample_name}_report.tsv" ]]; then
					plasFlow_Stats="Found"
				else
					plasFlow_Stats="NOT_FOUND"
				fi
			fi
		elif [[ "${family}" == "Unknown" ]]; then
			plasFlow="No_family_in_TAX_file"
			cplas="No_plasFlow"
			pFin_plas="No_plasFlow"
			plasFlow_Stats="No_plasFlow"
		else
			plasFlow="Not_ENTEROBACTERIACEAE_family"
			cplas="No_plasFlow"
			pFin_plas="No_plasFlow"
			plasFlow_Stats="No_plasFlow"
		fi

		if [[ -s "${processed}/${project}/${sample_name}/ANI/${genus}_and_${sample_name}_mashtree.dnd" ]]; then
			animash="Found"
		else
			animash="NOT_FOUND"
		fi
		if [[ -s "${processed}/${project}/${sample_name}/ANI/best_ANI_hits_ordered(${sample_name}_vs_${genus}).txt" ]]; then
			anigenus="Found"
		else
			anigenus="NOT_FOUND"
		fi

		if [[ "${genus}" == "Acinetobacter" ]] && [[ "${species}" == "baumannii" ]]; then
			if [[ -s "${processed}/${project}/${sample_name}/MLST/${sample_name}.mlst" ]]; then
				mlst_result1=$(head -n1 "${processed}/${project}/${sample_name}/MLST/${sample_name}.mlst" | cut -d'	' -f3)
			else
				mlst_result1="NOT_FOUND"
			fi
			if [[ -s "${processed}/${project}/${sample_name}/MLST/${sample_name}_abaumannii.mlst" ]]; then
				mlst_result2=$(head -n1 "${processed}/${project}/${sample_name}/MLST/${sample_name}_abaumannii.mlst" | cut -d'	' -f3)
			else
				mlst_result2="NOT_FOUND"
			fi
			mlst="${mlst_result1}|${mlst_result2}"
			if [[ -s "${processed}/${project}/${sample_name}/MLST/${sample_name}_srst2_Acinetobacter_baumannii#1-Oxford.mlst" ]]; then
				srst2_mlst_result1=$(tail -n1 "${processed}/${project}/${sample_name}/MLST/${sample_name}_srst2_Acinetobacter_baumannii#1-Oxford.mlst" | cut -d'	' -f2)
			else
				srst2_mlst_result1="NOT_FOUND"
			fi
			if [[ -s "${processed}/${project}/${sample_name}/MLST/${sample_name}_srst2_Acinetobacter_baumannii#2-Pasteur.mlst" ]]; then
				srst2_mlst_result2=$(tail -n1 "${processed}/${project}/${sample_name}/MLST/${sample_name}_srst2_Acinetobacter_baumannii#2-Pasteur.mlst" | cut -d'	' -f2)
			else
				srst2_mlst_result2="NOT_FOUND"
			fi
			srst2_mlst="${mlst_result1}|${mlst_result2}"
		elif [[ "${genus}" == "Escherichia" ]] && [[ "${species}" == "coli" ]]; then
			if [[ -s "${processed}/${project}/${sample_name}/MLST/${sample_name}.mlst" ]]; then
				mlst_result1=$(head -n1 "${processed}/${project}/${sample_name}/MLST/${sample_name}.mlst" | cut -d'	' -f3)
			else
				mlst_result1="NOT_FOUND"
			fi
			if [[ -s "${processed}/${project}/${sample_name}/MLST/${sample_name}_ecoli_2.mlst" ]]; then
				mlst_result2=$(head -n1 "${processed}/${project}/${sample_name}/MLST/${sample_name}_ecoli_2.mlst" | cut -d'	' -f3)
			else
				mlst_result2="NOT_FOUND"
			fi
			mlst="${mlst_result1}|${mlst_result2}"
			if [[ -s "${processed}/${project}/${sample_name}/MLST/${sample_name}_srst2_Escherichia_coli#1-Achtman.mlst" ]]; then
				srst2_mlst_result1=$(tail -n1 "${processed}/${project}/${sample_name}/MLST/${sample_name}_srst2_Escherichia_coli#1-Achtman.mlst" | cut -d'	' -f2)
			else
				srst2_mlst_result1="NOT_FOUND"
			fi
			if [[ -s "${processed}/${project}/${sample_name}/MLST/${sample_name}_srst2_Escherichia_coli#2-Pasteur.mlst" ]]; then
				srst2_mlst_result2=$(tail -n1 "${processed}/${project}/${sample_name}/MLST/${sample_name}_srst2_Escherichia_coli#2-Pasteur.mlst" | cut -d'	' -f2)
			else
				srst2_mlst_result2="NOT_FOUND"
			fi
			srst2_mlst="${srst2_mlst_result1}|${srst2_mlst_result2}"
		else
			if [[ -s "${processed}/${project}/${sample_name}/MLST/${sample_name}_${genus}_${species}.mlst" ]]; then
				srst2_mlst=$(tail -n1 "${processed}/${project}/${sample_name}/MLST/${sample_name}_${genus}_${species}.mlst" | cut -d'	' -f2)
			else
				srst2_mlst="NOT_FOUND"
			fi
		fi
	else
		plasFlow="No_TAX_file"
		cplas="No_plasFlow"
		pFin_plas="No_plasFlow"
		srst2_mlst="NO_TAX_file"
		animash="NO_TAX_file"
		anigenus="NO_TAX_FILE"
	fi

	echo -e "${counter}:${project}/${sample_name}:	${FQR1}	${FQZR1}	${FQR2}	${FQZR2}	${FQTR1}	${FQTZR1}	${FQTR2}	${FQTZR2}	${preQCr}	${preQCt}	${krakr}	${gott}	${Assembly}	${Assembly_stats}	${kraka}	${busco}	${prokka}	${tax}	${family}	${genus}	${species}	${sixteenS}	${animash}	${anigenus}	${aniAll}	${mlst}	${srst2_mlst}	${pFin}	${ohsixoheight}	${ohsixoheights}	${input_DB_csstar}	${input_DB_srst2}	${plasFlow}	${plasFlow_Stats}	${cplas}	${pFin_plas}"
	echo -e "${project}/${sample_name}:	${FQR1}	${FQZR1}	${FQR2}	${FQZR2}	${FQTR1}	${FQTZR1}	${FQTR2}	${FQTZR2}	${preQCr}	${preQCt}	${krakr}	${gott}	${Assembly}	${Assembly_stats}	${kraka}	${busco}	${prokka}	${tax}	${family}	${genus}	${species}	${sixteenS}	${animash}	${anigenus}	${aniAll}	${mlst}	${srst2_mlst}	${pFin}	${ohsixoheight}	${ohsixoheights}	${input_DB_csstar}	${input_DB_srst2}	${plasFlow}	${plasFlow_Stats}	${cplas}	${pFin_plas}" >> "${3}"
	counter=$(( counter + 1 ))
done < "${1}"
echo "All isolates completed"
global_end_time=$(date "+%m-%d-%Y @ %Hh_%Mm_%Ss")
#Script exited gracefully (unless something else inside failed)
printf "%s %s" "act_by_list_AR_completion_check.sh has completed check for AR completion " "${global_end_time}" | mail -s "act_by_list complete" nvx4@cdc.gov
exit 0
