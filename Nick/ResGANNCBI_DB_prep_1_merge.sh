#!/bin/sh -l

#$ -o prepCSSTARDB_1.out
#$ -e prepCSSTARDB_1.err
#$ -N prepCSSTARDB_1
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
. ./config.sh
#Import list of modds used during pipeline analysis (or downstream)
module load Python/2.7.15
. "${mod_changers}/list_modules.sh"

#
# Consolidates resFinders multi fasta to one
#
# Usage ./ResGANNCBI_DB_prep_1.sh path_to_resFinder_zip	path_to_ARGANNOT_fasta	path_to_output_dir
#

$(python2 -V)

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, using DEFAULT DATADIR=${local_DBs}/star/db_prep"
	DATADIR="${local_DBs}/star/db_prep"
	if [ -d "${DATADIR}" ]; then
		#echo "Directory already exists, try another or delete ${DATADIR}, exiting"
		#exit
		:
	else
		mkdir "${DATADIR}"
	fi
else
	if [[ "${1}" == "-f" ]]; then
		if [[ -z "${2}" ]]; then
			echo "Empty filename supplied to $0, exiting"
			exit 1
		else
			non_duplicated="true"
			ResGANNCBI_source="${DATADIR}/${2}"
		fi
	# Gives the user a brief usage and help section if requested with the -h option argument
	elif [[ "${1}" = "-h" ]]; then
		echo "Usage is ./ResGANNCBI_DB_prep_1.sh path_to_resFinder_zip	path_to_ARGANNOT_fasta	(path_to_output_dir)"
		echo "Converts the ARGANNOT fasta and resFinder ZIP files into a single fasta file to a single srst2 db fasta for use with csstar and srst2"
		echo "Output is ResGANNCBI_date_srst2.fasta"
		exit 0
	elif [[ ! -z ${2} ]]; then
		if [[ -f "${1}" ]] && [[ -f "${2}" ]]; then
			NCBI_source=${3}
			ARGANNOT_source=${2}
			resFinder_source=${1}
			if [[ ! -z "${3}" ]]; then
				DATADIR="${3}"
				if [[ ! -d "${DATADIR}" ]] ;then
					mkdir -p ${DATADIR}
				fi
			fi
		else
			echo "One of the source files does not exist"
		fi
	fi
fi
date
today=$(date '+%Y%m%d')

# Skips the automatic consolidation dual source DBs, in favor of an already created DB (does not check for duplicates if skipped)
if [[ "${non_duplicated}" != "true" ]]; then
	# #FASTA files containing most up to date databases from resFinder and ARGANNOT
	#echo ":${ARGANNOT_source}:${resFinder_source}:"
	if [[ "${ARGANNOT_source}" == "" ]]; then
		#echo "looking ar:${DATADIR}"
		ARGANNOT_source=$(find ${DATADIR} -name 'arg-annot-nt*')
	fi
	echo "::${ARGANNOT_source}:${resFinder_source}:${NCBI_source}::"

	cp "${ARGANNOT_source}" "${DATADIR}/argannot_${today}.fasta"
	ARGANNOT_source="${DATADIR}/argannot_${today}.fasta"
	if [[ ! -f "${ARGANNOT_source}" ]]; then
		echo "No argannot fasta found, exiting"
	fi
	sed -i 's/\//_/g' "${ARGANNOT_source}"

	if [[ "${resFinder_source}" == "" ]]; then
	#echo "Looking rf:${DATADIR}"
		resFinder_zip=$(find ${DATADIR} -name 'genomicepidemiology-resfinder_db-*')
	fi
	echo "::${ARGANNOT_source}:${resFinder_zip}:${NCBI_source}::"

	temp_dir=$(basename "${resFinder_zip}" ".zip")
	#echo "${resFinder_zip}, ${temp_dir}"
	unzip "${resFinder_zip}" "-d" "${DATADIR}"

	for files in ${DATADIR}/${temp_dir}/*
	do
		echo "${files}"
		mv "${files}" "${DATADIR}"
	done

	if [[ "${NCBI_source}" == "" ]]; then
		#echo "looking ar:${DATADIR}"
		NCBI_source=$(find ${DATADIR} -name 'NCBI*')
		cp "${NCBI_source}" "${DATADIR}/NCBI_${today}.fasta"
		NCBI_source="${DATADIR}/NCBI_${today}.fasta"
		if [[ ! -f  "${NCBI_source}" ]]; then
			echo "No NCBI fasta found, exiting"
		fi
	fi

	echo "::${ARGANNOT_source}:${resFinder_zip}:${NCBI_source}::"


	resFinder_source="${DATADIR}/resFinder_${today}.fasta"
	ResGANNCBI_source="${DATADIR}/ResGANNCBI_${today}.fasta"

	echo "arg-source=${ARGANNOT_source}"
	echo "res-source=${resFinder_source}"
	echo "NCBI-source=${NCBI_source}"
	echo "resGANNCBI-source=${ResGANNCBI_source}"
	#Consolidate all resFinder excel files into a single large fasta file
	#for file in ${DATADIR}/*.fsa
	#do
	#	while IFS= read -r line || [ -n "$line" ]; do
	#		echo ${line} >> ${resFinder_source}
	#	done < ${file}
	#done
	cat ${DATADIR}/*.fsa > ${resFinder_source}
	rm -r ${DATADIR}/*.fsa
	rm -r ${DATADIR}/${temp_dir}

	python "${shareScript}/ResGANNOT_Combiner_Non_Duplicate_Exe.py" "${resFinder_source}" "${ARGANNOT_source}" "${NCBI_source}" "${ResGANNCBI_source}" "${DATADIR}/Copies.fasta" "RES" "ARG" "NCBI"
fi

#Creates an associative array for looking up the genes to what they confer
declare -A groups
echo "Creating reference array"
counter=0
while IFS= read -r line || [ -n "$line" ]; do
	line=${line,,}
	gene=$(echo "${line}" | cut -d ':' -f1)
	first=${gene:0:1}
	if [ "$first" == "#" ]; then
		continue
	fi
	confers=$(echo "${line}" | cut -d ':' -f2)
	groups[${gene}]="${confers}"
	#echo "${counter}:${gene}:${confers}"
	counter=$(( counter + 1))
done < "${local_DBs}/star/group_defs.txt"
#echo "${#groups[@]}"

for gene in "${!groups[@]}"; do
	echo "${gene}:${groups[${gene}]}:" >> "${local_DBs}/star/test_groups.txt"
done

# unique_groups=($(echo "${groups[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))

# for gene in ${unique_groups[@]}; do
	# echo ${gene} >> "${local_DBs}/star/test_groups_after.txt"
# done


echo "checking for new groups in the source DB fasta file"
new_groups=()
while IFS= read -r line || [ -n "$line" ]; do
	#echo "=${line}="
	if [[ "${line:0:1}" == ">" ]]; then
		source=$(echo "${line:2}" | cut -d']' -f1)
		header=$(echo "${line}" | cut -d']' -f2)
		if [[ "${source}" == "ARG" ]]; then
			#echo "|${header}|"
			gene_type=$(echo ${header:1} | cut -d')' -f1)
			#echo "${gene_type,,}"
			confers="${groups[${gene_type,,}]}"
			#echo "confers-${confers}-"
			if [[ "${confers}" ]];then
				continue
			else
				echo "-${line}-"
				if [[ ! -z "${line}" ]] && [[ "${header:0:1}" == "(" ]]; then
					new_groups+=(${gene_type,,})
					#read -p "${gene_type,,} is not in the group definitions...what resistance does it confer?" new_conf
					#echo "OK, ${gene_type,,} will confer ${new_conf} resistance"
					#groups[${gene_type,,}]="${new_conf}"
					#echo -e  "${gene_type,,}:${new_conf}_resistance:" >> "${local_DBs}/star/group_defs.txt"
				else
					echo "${line} is malformed, please investigate and rerun"
					exit
				fi
			fi
		elif [[ "${source}" == "RES" ]]; then
			# This is the one case...so far....where resFinder has an issue with the first 3 letters being able to determine resistance
			if [[ "${header:0:4}" = "cmrA" ]]; then
				gene_type=${header:0:4}
				gene_type=${gene_type,,}
			else
				gene_type=${header:0:3}
				gene_type=${gene_type,,}
			fi
			#echo "${gene_type}"
			confers="${groups[${gene_type}]}"
			#echo "confers-${confers}-"
			if [[ "${confers}" ]]; then
				continue
			else
				echo "RES -${line}-"
				new_groups+=(${gene_type})
			fi
		elif [[ "${source}" == "NCBI" ]]; then
			gene_type=$(echo ${header} | cut -d'_' -f1)
			gene_type=${gene_type:0:3}
			gene_type=${gene_type,,}
			#echo "${gene_type}"
			confers="${groups[${gene_type}]}"
			#echo "confers-${confers}-"
			if [[ "${confers}" ]]; then
				continue
			else
				echo "NCBI -${line}-"
				new_groups+=(${gene_type})
			fi
		fi
	fi
done < ${ResGANNCBI_source}

echo "Assigning definitions to new group members"
for i in ${new_groups[@]}
do
	if [[ ! -z "${i}" ]]; then
		read -p "${i} is not in the group definitions...what resistance does it confer?" new_conf
		echo "OK, ${i} will confer ${new_conf} resistance"
		groups[${i}]="${new_conf}"
		echo -e  "${i}:${new_conf}_resistance:" >> "${local_DBs}/star/group_defs.txt"
	fi
done

# Create the gene lookup file to match gene and conferred resistance downstream
echo "There are ${#groups[@]} different groups"
echo -e "Please be sure to check ${local_DBs}/star/group_defs.txt for proper assignment of resistance to new groups\nand ${DATADIR}/ResGANNCBI_${today}.bad for any actual sequences that might be able to be added"

exit
