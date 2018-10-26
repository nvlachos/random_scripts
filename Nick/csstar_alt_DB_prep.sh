#!/bin/sh -l

#$ -o prepCSSTARDB.out
#$ -e prepCSSTARDB.err
#$ -N prepCSSTARDB
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
# Usage ./run_resPrep.sh path_to_dir
#

$(python2 -V)

DATADIR="/scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/DBs/star/db_prep"
DB_source=${1}

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to csstar_DB_prep.sh"
	exit
else
	if [[ -z "${2}" ]]; then
		echo "Empty path supplied to csstar_DB_prep.sh, exiting"
		exit 1
	# Gives the user a brief usage and help section if requested with the -h option argument
	elif [[ "${1}" = "-h" ]]; then
		echo "Usage is ./DB_prep.sh input_fasta output_directory"
		echo "Output is ResGANNOT_date_srst2.fasta"
		exit 0
	else
		DATADIR=${2}
		if [[ ! -d "${DATADIR}" ]]; then
			mkdir -p "${DATADIR}"
		fi
		DB_source="${1}"
		DB_short_name=$(basename ${1})
	fi
fi
date

today=$(date '+%Y%m%d')

 module load CD-HIT/4.6
# #switch to python 2.7 temporarily for cdhit to csv scripts/config
 cp "${share}/DBs/star/cdhit_to_csv.py" "${DATADIR}"
 cp "${share}/DBs/star/csv_to_gene_db.py" "${DATADIR}"
 
 cd "${DATADIR}"

 #module unload Python/3.5.2
 module load Python/2.7.13
 python2 "-V"

echo "--- About to CD-HIT-EST on ResGANNOT DB ---"
cd-hit-est "-i" "${DB_source}" "-o" "${DATADIR}/${DB_short_name}_${today}_cdhit90" "-d" "0" > "${DATADIR}/${DB_short_name}_${today}_cdhit90.stdout"

echo "--- About to CD-HIT-to-CSV on ResGANNOT DB ---"
python "${DATADIR}/cdhit_to_csv.py" "--cluster_file" "${DATADIR}/${DB_short_name}_${today}_cdhit90.clstr" "--infasta" "${DB_source}" "--outfile" "${DATADIR}/${DB_short_name}_${today}_clustered.csv"
# python /scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/DBs/star/db_prep/cdhit_to_csv.py --cluster_file /scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/DBs/star/db_prep/ResGANNOT_040518_cdhit90.clstr --infasta /scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/DBs/star/db_prep/Combined.fasta --outfile /scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/DBs/star/db_prep/ResGANNOT_040518_clustered.csv



#Creates an associative array for matching accession to dna sequence
declare -A seqarr
counter=0
count=0
seq=""
echo "Creating seq reference array"
while IFS= read -r line;
do
	line=$(echo ${line} | tr -d '\040\011\012\015')
	if [[ "${line:0:1}" == ">" ]]; then
		if [[ "${counter}" -gt 0 ]]; then
			#echo "${header}:${seq}" >> ${DATADIR}/hits.txt
			seqarr[${header}]="${seq}"
			seq=""
		fi
		header="${line:1}"
		count=$(( count + 1 ))
		#echo "Adding ${count} to list"
	else
		seq="${seq}${line}"
	fi
	counter=$(( counter + 1 ))
	echo "${counter}"
done < "${DB_source}"

#add the last sequence to the array, since there is no > to trigger its addition
#echo "${header}:${seq}" >> ${DATADIR}/hits.txt
echo "Adding final item, number ${count}"
seqarr[${header}]="${seq}"
echo "There are ${#seqarr[@]} different entries, and counter equals ${count}"

counter=0
echo "Before"
while IFS= read -r gene_line
do
	echo "It: ${counter}"
	if [[ ${counter} -eq 0 ]]; then
		seqID="seqID"
		clusterID="clusterID"
		gene="gene"
		allele="allele"
		dnaseq="DNA"
		accession="accession"
	else
		echo "Test-${gene_line}"
		IFS=',' read -ra line_items <<< "${gene_line}"
		seqID=${line_items[0]}
		clusterID=${line_items[1]}
		gene=$(echo "${line_items[3]}" | cut -d':' -f2 | cut -d'_' -f2- | cut -d'_' -f1)
		#allele=$(echo ${line_items[3]} | cut -d':' -f2 | cut -d'_' -f2- | cut -d'_' -f2)
		accession=$(echo "${line_items[3]}" | cut -d':' -f1)
		dnaseq="${seqarr[${line_items[3]}]}"
		#echo "result:${dnaseq}:"
	fi
	echo "${seqID}/${clusterID}/${gene}/${allele}/${accession}/${group}"
#	echo "${dnaseq}"
	echo -e "${seqID},${clusterID},${gene},${gene},${dnaseq},${accession},,," >> "${DATADIR}/${DB_short_name}_${today}_cluster90_edited.csv"
	counter=$(( counter + 1 ))
done < ${DATADIR}/${DB_short_name}_${today}_clustered.csv
echo "After"

echo "--- About to csv to gene ---"
python "${DATADIR}/csv_to_gene_db.py" "-t" "${DATADIR}/${DB_short_name}_${today}_cluster90_edited.csv" "-o" "${DATADIR}/${DB_short_name}_${today}_srst2.fasta" "-f" "${DB_source}" "-c" "6" "-s" "5"

exit
