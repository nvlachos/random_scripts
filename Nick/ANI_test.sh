#!/bin/sh -l

#$ -o test_ANI.out
#$ -e test_ANI.err
#$ -N test_ANI
#$ -cwd
#$ -q all.q

pwd
#Import the config file with shortcuts and settings
if [[ ! -d ./config.sh ]]; then
	cp config_template.sh config.sh
fi
. ./config.sh
#   ${shareScript}/module_changers/list_modules.sh

#
# Script to calculate the average nucleotide identity of a sample to numerous other samples from the same genus (genus dependent)
# The most similar match is identified and provided for confirmation
#
# Usage ./run_ANI.sh sample_name	run_id
#
# Python/3.5.2 (pyani is located in Nick_DIR/script folder, not run from scicomp module)
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to run_ANI.sh, exiting"
	exit 1
elif [[ -z "${1}" ]]; then
	echo "Empty sample name supplied to run_ANI.sh, exiting"
	exit 1
# Gives the user a brief usage and help section if requested with the -h option argument
elif [[ "${1}" = "-h" ]]; then
	echo "Usage is ./run_ANI.sh sample_name run_id"
	echo "Output is saved to in ${processed}/sample_name/ANI"
	exit 0
#elif [ -z "$2" ]; then
#	echo "Empty database name supplied to run_ANI.sh. Second argument should be a genus found in ${share}/DBs/ANI/  ...Exiting"
#	exit 1
fi

start_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")
echo "Started ANI at ${start_time}"

# Sets the output folder to the sample_name folder in processed samples
OUTDATADIR="${local_DBs}/aniDB"
if [[ ! -d ${OUTDATADIR}/all_named_test/dists ]]; then
	mkdir ${OUTDATADIR}/all_named_test/dists
fi

# Gets persons name to use as email during entrez request to identify best matching sample
me=$(whoami)
#echo ${me}"___"${1}___${2}___${3}___${4}

#Creates a local copy of the database folder
# cp ${share}/DBs/aniDB/all/compound_sketch_all.msh "${OUTDATADIR}/ANI/"

#mash dist "${local_DBs}/aniDB/refseq.genomes.k21s1000.msh" "${OUTDATADIR}/Assembly/${1}_scaffolds_trimmed.fasta" > "${OUTDATADIR}/ANI/${1}_all_refSeq.dists"
counter=0
for ref in ${local_DBs}/aniDB/all_named_test/*.fna; do
	echo ${ref}
	counter=$(( counter + 1 ))
	filename=$(basename ${ref})
	mash dist "${local_DBs}/aniDB/all_named_test/all_named.msh" "${ref}" > "${OUTDATADIR}/all_named_test/dists/${filename}_unsorted.dists"
	sort -k3 -n -o "${OUTDATADIR}/all_named_test/dists/${filename}.dists" "${OUTDATADIR}/all_named_test/dists/${filename}_unsorted.dists"
	rm -r "${OUTDATADIR}/all_named_test/dists/${filename}_unsorted.dists"
done

for distfile in ${local_DBs}/aniDB/all_named_test/*.dists; do
	taxa=$(basename ${distfile} | cut -d'_' -f1,2)
	if [[ ! -d ${local_DBs}/aniDB/all_named_test/${taxa} ]]; then
		mkdir -p ${local_DBs}/aniDB/all_named_test/${taxa}/localANIDB
	fi
	counter=0
	max_ani_samples=30
	> "${local_DBs}/aniDB/all_named_test/${taxa}/thirty_closest_dists.txt"
	while IFS= read -r line || [[ "$line" ]];  do
		if [[ ! -d ${local_DBs}/aniDB/all_named_test/${taxa}/localANIDB ]]; then
			if [[ ${counter} -eq 0 ]]; then
				ref_path=$(echo "${line}" | cut -d'	' -f2)
				"${ref_path}" >> "${local_DBs}/aniDB/all_named_test/${taxa}/thirty_closest_dists.txt"
				cp ${ref_path} ${local_DBs}/aniDB/all_named_test/${taxa}/localANIDB
			fi
			if [[ ${counter} -gt ${max_ani_samples} ]]; then
				break
			else
				source_path=$(echo "${line}" | cut -d'	' -f1)
				cp ${source_path} ${local_DBs}/aniDB/all_named_test/${taxa}/localANIDB
		fi
	done < ${distfile}
done
echo ${counter}
exit

#counter=0
#threshold=0
#cp "${OUTDATADIR}/Assembly/${1}_scaffolds_trimmed.fasta" "${OUTDATADIR}/ANI/localANIDB/sample_${1}.fasta"
#> "${OUTDATADIR}/ANI/twenty_closest_mash.list"
#while IFS='' read -r line;
#do
#	source_path=$(echo "${line}" | cut -d'	' -f1)
#	source=$(echo "${source_path}" | rev | cut -d'/' -f1 | rev)
#	source=${source:0:-4}
#	echo ${source}
#	distance=$(echo "${line}" | cut -d'	' -f3)
#	if [[ ${counter} -gt 50 ]] && [[ ${distance} != ${threshold} ]]; then
#		break
#	else
#		echo "Going to copy ${source} to localANIDB"
#		cp "${source_path}" "${OUTDATADIR}/ANI/localANIDB/${source}.fasta"
#		echo "${OUTDATADIR}/ANI/localANIDB/${source}.fasta" >> "${OUTDATADIR}/ANI/twenty_closest_mash.list"
#	fi
#	counter=$(( counter + 1))
#	threshold="${distance}"
#done < "${OUTDATADIR}/ANI/${1}_all_sorted.dists"

# Checks for a previous copy of the aniM folder, removes it if found
if [ -d "${OUTDATADIR}/ANI/aniM" ]; then  #checks for and removes old results folder for ANIm
	echo "Removing old ANIm results in ${OUTDATADIR}/ANI/aniM"
	rm -rf "${OUTDATADIR}/ANI/aniM"
fi

python "/apps/x86_64/pyani/pyani/bin/average_nucleotide_identity.py" -i "${local_DBs}/aniDB/all_named_test/${taxa}/localANIDB" -o "${local_DBs}/aniDB/all_named_test/${taxa}/aniM" --write_excel
python "/apps/x86_64/pyani/pyani/bin/average_nucleotide_identity.py" -i "${local_DBs}/aniDB/all_named_test/${taxa}/localANIDB" -o "${local_DBs}/aniDB/all_named_test/${taxa}/aniB" --write_excel






#Extracts the query sample info line for percentage identity from the percent identity file
while IFS='' read -r line;
do
#	echo "!-${line}"
	if [[ ${line:0:7} = "sample_" ]]; then
		sampleline=${line}
#		echo "found it-"$sampleline
		break
	fi
done < "${local_DBs}/aniDB/all_named_test/${taxa}/aniM/ANIm_percentage_identity.tab"

#Extracts the top line from the %id file to get all the sample names used in analysis (they are tab separated along the top row)
if [[ -s "${local_DBs}/aniDB/all_named_test/${taxa}/aniM/ANIm_percentage_identity.tab" ]]; then
	firstline=$(head -n 1 "${local_DBs}/aniDB/all_named_test/${taxa}/aniM/ANIm_percentage_identity.tab")
else
	echo "No "${local_DBs}/aniDB/all_named_test/${taxa}/aniM/ANIm_percentage_identity.tab" file, exiting"
	exit 1
fi

#Arrays to read sample names and the %ids for the query sample against those other samples
IFS="	" read -r -a samples <<< "${firstline}"
IFS="	" read -r -a percents <<< "${sampleline}"

#How many samples were compared
n=${#samples[@]}

#Extracts all %id against the query sample (excluding itself) and writes them to file
for (( i=0; i<n; i++ ));
do
#	echo ${i}-${samples[i]}
	if [[ ${samples[i]:0:7} = "sample_" ]];
	then
#		echo "Skipping ${i}"
		continue
	fi
	definition=$(head -1 "${local_DBs}/aniDB/all_named_test/${taxa}/localANIDB/${samples[i]}.fna")
	# Prints all matching samples to file (Except the self comparison) by line as percent_match  sample_name  fasta_header
	echo "${percents[i+1]} ${samples[i]} ${definition}" >> "${local_DBs}/aniDB/all_named_test/${taxa}/best_hits_aniM.txt"
done

#Sorts the list in the file based on %id (best to worst)
sort -nr -t' ' -k1 -o "${local_DBs}/aniDB/all_named_test/${taxa}/best_hits__aniM_ordered.txt" "${local_DBs}/aniDB/all_named_test/${taxa}/best_hits_aniM.txt"

#Extracts the query sample info line for percentage identity from the percent identity file
while IFS='' read -r line;
do
#	echo "!-${line}"
	if [[ ${line:0:7} = "sample_" ]]; then
		sampleline=${line}
#		echo "found it-"$sampleline
		break
	fi
done < "${local_DBs}/aniDB/all_named_test/${taxa}/aniB/ANIb_percentage_identity.tab"

#Extracts the top line from the %id file to get all the sample names used in analysis (they are tab separated along the top row)
if [[ -s "${local_DBs}/aniDB/all_named_test/${taxa}/aniB/ANIb_percentage_identity.tab" ]]; then
	firstline=$(head -n 1 "${local_DBs}/aniDB/all_named_test/${taxa}/aniB/ANIb_percentage_identity.tab")
else
	echo "No "${local_DBs}/aniDB/all_named_test/${taxa}/aniB/ANIb_percentage_identity.tab" file, exiting"
	exit 1
fi

#Arrays to read sample names and the %ids for the query sample against those other samples
IFS="	" read -r -a samples <<< "${firstline}"
IFS="	" read -r -a percents <<< "${sampleline}"

#How many samples were compared
n=${#samples[@]}

#Extracts all %id against the query sample (excluding itself) and writes them to file
for (( i=0; i<n; i++ ));
do
#	echo ${i}-${samples[i]}
	if [[ ${samples[i]:0:7} = "sample_" ]];
	then
#		echo "Skipping ${i}"
		continue
	fi
	definition=$(head -1 "${local_DBs}/aniDB/all_named_test/${taxa}/localANIDB/${samples[i]}.fna")
	# Prints all matching samples to file (Except the self comparison) by line as percent_match  sample_name  fasta_header
	echo "${percents[i+1]} ${samples[i]} ${definition}" >> "${local_DBs}/aniDB/all_named_test/${taxa}/best_hits_aniB.txt"
done

#Sorts the list in the file based on %id (best to worst)
sort -nr -t' ' -k1 -o "${local_DBs}/aniDB/all_named_test/${taxa}/best_hits__aniB_ordered.txt" "${local_DBs}/aniDB/all_named_test/${taxa}/best_hits_aniB.txt"

end_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")
echo "ENDed ANI at ${end_time}"

#Script exited gracefully (unless something else inside failed)
exit 0
