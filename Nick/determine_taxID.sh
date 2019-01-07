#!/bin/sh -l

#$ -o getTax.out
#$ -e getTax.err
#$ -N getTax
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi
. ./config.sh

#
# Creates a single file that attempts to pull the best taxonomic information from the isolate. Currently, it operates in a linera fashion, e.g. 1.ANI, 2.kraken, 3.Gottcha, 4.16s
# The taxon is chosen based on the highest ranked classifier first
#
# Usage ./determine_texID.sh sample_name project_ID
#
# No modules required
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied todetermine_taxID.sh, exiting"
	exit 1
elif [[ -z "${1}" ]]; then
	echo "Empty sample_id supplied to determine_taxID.sh, exiting"
	exit 1
elif [[ "${1}" = "-h" ]]; then
	echo "Usage is ./determine_taxID.sh sample_ID run_ID"
	echo "Output is saved to ${processed}/run_ID/sample_ID/taxonomy.csv"
	exit 0
elif [[ -z "${2}" ]]; then
	echo "Empty run_id supplied to determine_taxID.sh, exiting"
	exit 1
fi

# Set default values after setting variables
sample=${1}
project=${2}

Domain="Not_assigned"
Phylum="Not_assigned"
Class="Not_assigned"
Order="Not_assigned"
Family="Not_assigned"
Genus="Not_assigned"
species="Not_assigned"
ani_files=0

if [[ $(ls -t "${processed}/${project}/${sample}/ANI/best_ANI_hits_ordered"* | wc -l) -ne 0 ]]; then
	ani_files=$(ls -t "${processed}/${project}/${sample}/ANI/best_ANI_hits_ordered"* | wc -l)
else
	ani_files=0
fi

if [[ ${ani_files} -gt 0 ]]; then
	source="ANI"
	echo "${source}"
	if [[ -f "${processed}/${project}/${sample}/ANI/best_ANI_hits_ordered(${sample}_vs_All).txt" ]]; then
		source_file="${processed}/${project}/${sample}/ANI/best_ANI_hits_ordered(${sample}_vs_All).txt"
	else
		source_file=$(ls -t "${processed}/${project}/${sample}/ANI/best_ANI_hits_ordered"* | head -n 1)
	fi
	header=$(head -n 1 "${source_file}")
	echo "${header}"
	Genus=$(echo "${header}" | cut -d' ' -f1 | cut -d'-' -f2)
	species=$(echo "${header}" | cut -d' ' -f2 | cut -d'(' -f1)
	echo "${Genus}-${species}"
elif [[ -s "${processed}/${project}/${sample}/16s/${sample}_16s_blast_id.txt" ]]; then
	source="16s"
	#echo "${source}"
		# Lookup Taxonomy
		line=$(tail -n 1 "${processed}/${project}/${sample}/16s/${sample}_16s_blast_id.txt")
		Genus=$(echo "${line}" | cut -d"	" -f3 | cut -d" " -f1)
		species=$(echo "${line}" | cut -d"	" -f3 | cut -d" " -f2)
elif [[ -s "${processed}/${project}/${sample}/gottcha/${sample}_gottcha_species_summary.txt" ]]; then
	source="GOTTCHA"
	#echo "${source}"
	while IFS= read -r line;
	do
		# Grab first letter of line (indicating taxonomic level)
		first=${line::1}
		# Assign taxonomic level value from 4th value in line (1st-classification level,2nd-% by kraken, 3rd-true % of total reads, 4th-identifier)
		if [ "${first}" = "S" ]
		then
			species=$(echo "${line}" | awk -F ' ' '{print $5}')
		elif [ "${first}" = "G" ]
		then
			Genus=$(echo "${line}" | awk -F ' ' '{print $4}')
		fi
	done < "${processed}/${project}/${sample}/gottcha/${sample}_gottcha_species_summary.txt"
elif [[ -s "${processed}/${project}/${sample}/kraken/postAssembly/${sample}_kraken_summary_assembled_BP_data.txt" ]]; then
	source="Kraken"
	#echo "${source}"
	while IFS= read -r line;
	do
		# Grab first letter of line (indicating taxonomic level)
		first=${line::1}
		# Assign taxonomic level value from 4th value in line (1st-classification level,2nd-% by kraken, 3rd-true % of total reads, 4th-identifier)
		if [ "${first}" = "S" ]
		then
			species=$(echo "${line}" | awk -F ' ' '{print $4}')
		elif [ "${first}" = "G" ]
		then
			Genus=$(echo "${line}" | awk -F ' ' '{print $4}')
		fi
	done < "${processed}/${project}/${sample}/kraken/postAssembly/${sample}_kraken_summary_assembled_BP_data.txt"
else
	echo "Exiting, no reliable taxonomy files exist"
	exit 1
fi

if [[ ${Genus} == "Peptoclostridium" ]]; then
	Genus="Clostridium"
fi

while IFS= read -r line;
do
	DB_genus=$(echo ${line} | cut -d"," -f1)
	#echo ":${Genus}:${DB_genus}:"
	if [[ "${Genus,}" = "${DB_genus}" ]]; then
			tax_DB="${local_DBs}/taxes.csv"
			Domain=$(echo "${line}" | cut -d"," -f2)
			Phylum=$(echo "${line}" | cut -d"," -f3)
			Class=$(echo "${line}" | cut -d"," -f4)
			Order=$(echo "${line}" | cut -d"," -f5)
			Family=$(echo "${line}" | cut -d"," -f6 | tr -d '\r' )
			#echo ":${Family}:"
			break
	fi
done < "${local_DBs}/taxes.csv"
printf "(${source}) \nD:	${Domain}\nP:	${Phylum}\nC:	${Class}\nO:	${Order}\nF:	${Family}\nG:	${Genus}\ns:	${species}\n" > "${processed}/${project}/${sample}/${sample}.tax"
