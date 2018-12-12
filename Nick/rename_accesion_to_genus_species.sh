#!/bin/sh -l

#$ -o mmshfolder.out
#$ -e mmshfolder.err
#$ -N mmshfolder
#$ -cwd
#$ -q all.q

#echo ":${1}:"
if [[ ! -d ${1} ]]; then
	echo "${1} does not exist"
fi
#if [[ ! -d ${2} ]]; then
#	mkdir -p ${2}
#fi

zipped_num=$(ls ${1}/*.gz | wc -l)
if [[ "${zipped_num}" -gt 0 ]]; then
	gunzip ${1}/*.gz
else
	echo "No zipped files found to unzip"
fi

counter=0
for j in ${1}/*.${2}; do
	filename=$(basename ${j} | cut -d'_' -f1,2)
	prefix=${filename:0:2}
	echo ":${prefix}:"
	if [[ "{$prefix}" = *"GC"* ]] || [[ "{$prefix}" = *"NZ"* ]]; then
		#echo "good-${j}"
		accession=$(head -n1 ${j} | cut -d' ' -f1 | cut -d'>' -f2)
		#echo "${j}" >> "${1}/good.list"
	elif [[ "${prefix}" = "_" ]]; then
		accession=$(echo ${j} | rev | cut -d'_' -f1,2 | rev | cut -d'.' -f1,2)
		echo "${accession}"
	else
		echo "Already completed or not a normal accession identifier-${j}"
		continue
		#echo "${j}" >> "${1}/bad.list"
	fi
	#echo ${accession}
	genus=""
	species=""
	sub_counter=0
	echo "Into loop with ${sub_counter}:${genus}:${species}"
	while [[ ${sub_counter} -lt 5 ]]; do
		echo "in"
		genus_species_info=$(python ./entrez_get_taxon_from_accession.py ${accession} nvx4@cdc.gov)
		genus=$(echo "${genus_species_info}" | cut -d' ' -f1)
		species=$(echo "${genus_species_info}" | cut -d' ' -f2)
		if [[ ${genus} ]] && [[ ${species} ]]; then
			break
		else
			sub_counter=$(( sub_counter + 1 ))
		fi
		sleep 2
	done
	echo "out"
	#echo ${genus_species}
	#echo "${genus_species}" >> "${1}/species.list"
	#echo "g-${genus};s-${species};a-${accession}"
	#exit
	echo "saving ${j} to ${1}/${genus}_${species}_${accession}.${3}"
	mv "${j}" "${1}/${genus}_${species}_${accession}.${2}"
	echo "${counter}"
	counter=$(( counter + 1))
done

echo ${command}
