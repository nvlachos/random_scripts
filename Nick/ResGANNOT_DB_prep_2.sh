#!/bin/sh -l

#$ -o prepCSSTARDB_2.out
#$ -e prepCSSTARDB_2.err
#$ -N prepCSSTARDB_2
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
. ./config.sh
#Import list of modds used during pipeline analysis (or downstream)
module load Python/2.7.13
. "${mod_changers}/list_modules.sh"

#
# Consolidates resFinders multi fasta to one
#
# Usage ./run_resPrep.sh path_to_dir
#

$(python2 -V)


today=$(date '+%Y%m%d')

DATADIR="/scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/DBs/star/db_prep"
ResGANNOT_source="${DATADIR}/ResGANNOT_${today}.fasta"

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to prep_RESGANNOT_DB.sh, using default DATADIR=/scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/DBs/star/db_prep"
	#if [ -d "${DATADIR}" ]; then
	#	rm -r "${DATADIR}"
	#fi
	#mkdir "${DATADIR}"
else
	if [[ -z "${1}" ]]; then
		echo "Empty path supplied to ResGANNOT_DB_prep.sh, exiting"
		exit 1
	# Gives the user a brief usage and help section if requested with the -h option argument
	elif [[ "${1}" = "-h" ]]; then
		echo "Usage is ./ResGANNOT_DB_prep.sh "
		echo "Converts the ARGANNOT fasta and resFinder ZIP files into a single fasta file to a single srst2 db fasta for use with csstar and srst2"
		echo "Output is ResGANNOT_date_srst2.fasta"
		exit 0
	elif [[ "${1}" = "-f" ]]; then
		if [[ -z "${2}" ]]; then
			echo "Empty filename supplied to ResGANNOT_DB_prep.sh, exiting"
			exit 1
		else
			sourced="true"
			ResGANNOT_source="${DATADIR}/${2}"
		fi
	else
		DATADIR="${1}"
	fi
fi
date


#Creates an associative array for looking up the genes to what they confer
declare -A groups
echo "Creating reference array"
while IFS= read -r line;
do
	line=${line,,}
	gene=$(echo "${line}" | cut -d ':' -f1)
	first=${gene:0:1}
	if [ "$first" == "#" ]; then
		continue
	fi
	confers=$(echo "${line}" | cut -d ':' -f2)
	groups[${gene}]="${confers}"
done < "${share}/DBs/star/group_defs.txt"



module load CD-HIT/4.6
# #switch to python 2.7 temporarily for cdhit to csv scripts/config
cp "${local_DBs}/star/cdhit_to_csv.py" "${DATADIR}"
cp "${local_DBs}/star/csv_to_gene_db.py" "${DATADIR}"

cd "${DATADIR}"

echo "--- About to CD-HIT-EST on ResGANNOT DB-${ResGANNOT_source} ---"
cd-hit-est "-i" "${ResGANNOT_source}" "-o" "${DATADIR}/ResGANNOT_${today}_cdhit90" "-d" "0" > "${DATADIR}/ResGANNOT_${today}_cdhit90.stdout"
# cd-hit-est -i /scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/DBs/star/db_prep/Combined.fasta -o /scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/DBs/star/db_prep/ResGANNOT_040518_cdhit90 -d 0 > /scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/DBs/star/db_prep/ResGANNOT_040518_cdhit90.stdout

echo "--- About to CD-HIT-to-CSV on ResGANNOT DB-${ResGANNOT_source} ---"

python "${DATADIR}/cdhit_to_csv.py" "--cluster_file" "${DATADIR}/ResGANNOT_${today}_cdhit90.clstr" "--infasta" "${ResGANNOT_source}" "--outfile" "${DATADIR}/ResGANNOT_${today}_clustered.csv"
# python /scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/DBs/star/db_prep/cdhit_to_csv.py --cluster_file /scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/DBs/star/db_prep/ResGANNOT_040518_cdhit90.clstr --infasta /scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/DBs/star/db_prep/Combined.fasta --outfile /scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/DBs/star/db_prep/ResGANNOT_040518_clustered.csv

rm -r "${DATADIR}/"*".fsa"

#edit the clustered.csv file to creipo
#	is seqID, clusterid, gene, allele, cluster contains multiple genes, gene found in multiple clusters
#	should be seqID, clusterid, gene, allele, DNA seqeunce, resFinder DB name, accession, anitbiotic class

if [[ -s "${DATADIR}/hits.txt" ]]; then
	rm -r "${DATADIR}/hits.txt"
fi

#Creates an associative array for matching accession to dna sequence
declare -A seqarr
counter=0
count=0
seq=""
echo "Creating ResGANNOT seq reference array"
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
done < "${ResGANNOT_source}"
#add the last sequence to the array, since there is no > to trigger its addition
#echo "${header}:${seq}" >> ${DATADIR}/hits.txt
echo "Adding final item, number ${count}"
seqarr[${header}]="${seq}"
echo "There are ${#seqarr[@]} different entries, and counter equals ${count}"

if [ -f ${DATADIR}/ResGANNOT_${today}_clustered.csv ]; then
	echo "found ${DATADIR}/ResGANNOT_${today}_clustered.csv"
fi




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
		DB_ID="Database_ID"
		accession="accession"
		group="antibiotic class"
		from_source="Database"
	else
		echo "Test-${gene_line}"
		IFS=',' read -ra line_items <<< "${gene_line}"
		source=$(echo "${line_items[2]}" | cut -d']' -f1 | cut -d'[' -f2)
		if [[ "${source}" == "ARG" ]]; then
			seqID=${line_items[0]}
			clusterID=${line_items[1]}
			gene=${line_items[2]}
			from_source="ARGANNOT"
			if [[ "${gene}" == *":"* ]]; then
				gene=$(echo ${gene} | cut -d':' -f1)
			fi
			for (( i=0; i<${#gene}; i++ ));
			do
				if [[ "${gene:${i}:1}" == ")" ]]; then
					i=$(( i + 1 ))
					gene=${gene:${i}}
					break
				fi
			done
			allele=$(echo "${line_items[3]}" | cut -d':' -f1 | cut -d')' -f2-)
			group_raw=$(echo "${line_items[2]}" | cut -d')' -f1 | cut -d'(' -f2)
			accession=$(echo "${line_items[3]}" | cut -d':' -f2)
			allele_location=$(echo "${line_items[3]}" | cut -d':' -f3)
			allele_length=$(echo "${line_items[3]}" | cut -d':' -f4)
			group="${groups[${group_raw,,}]}"
			DB_ID="${line_items[3]}"
			#echo "looking up sequence:${line_items[3]}:"
			dnaseq="${seqarr[${line_items[3]}]}"
			#echo "result:${dnaseq}:"
		# Create a new way to identify other DB entries
		else
			seqID=${line_items[0]}
			clusterID=${line_items[1]}
			gene=$(echo "${line_items[2]}" | cut -d']' -f2 | cut -d'_' -f1)
			if [[ "${source}" == "RES" ]]; then
				from_source="ResFinder"
			else
				from_source="${source}"
			fi
			allele=$(echo ${line_items[3]} | cut -d'_' -f1 | cut -d']' -f2)
			DB_ID=${line_items[3]}
			group_lookup=${gene:0:3}
			#echo "Looking up group ${group_lookup,,} from ${line_items[3]}"
			group_lookup=${group_lookup,,}
			accession=$(echo "${line_items[3]}" | rev | cut -d'_' -f1 | rev)
			group=$(echo "${groups[${group_lookup}]}")
			#echo "looking up sequence:${line_items[3]}:"
			dnaseq="${seqarr[${line_items[3]}]}"
			#echo "result:${dnaseq}:"
		fi
	fi
	echo "${seqID}/${clusterID}/${gene}/${allele}/${accession}/${group}"
#	echo "${dnaseq}"
	echo -e "${seqID},${clusterID},${gene},${allele},${dnaseq},${DB_ID},${accession},${group},${from_source}" >> "${DATADIR}/ResGANNOT_${today}_cluster90_edited.csv"
	counter=$(( counter + 1 ))
done < ${DATADIR}/ResGANNOT_${today}_clustered.csv
echo "After"

cp "${DATADIR}/ResGANNOT_${today}_cluster90_edited.csv" "/scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/DBs/star/"

echo "--- About to csv to gene ---"
python "${DATADIR}/csv_to_gene_db.py" "-t" "${DATADIR}/ResGANNOT_${today}_cluster90_edited.csv" "-o" "${DATADIR}/ResGANNOT_${today}_srst2.fasta" "-f" "${ResGANNOT_source}" "-c" "6" "-s" "5"

#Additional Enumerator, specfically to address fos gene ambiguity
# declare -A allele_tally
# while IFS= read -r line;
# do
	# srst2=${echo ${line} | cut -d' ' -f1)
	# info1=$(echo ${srst2} | awk -v delimeter="__" '{split($0,a,delimeter)} END{print a[1]}')
	# supp_info=$(echo ${line} | cut -d' ' -2)
	# db=$(echo ${supp_info} | cut -d']' -f1 | cut -d'[' -f2)
	# if [[ "${db}" == "ARG" ]]; then
		# acc_info=$(echo ${line}
	# elif [[ "${db}" == "RES" ]]; then

	# else
		# echo "Unknown data, not editing line"
	# fi
# done < "${DATADIR}/ResGANNOT_${today}_srst2.fasta"






exit
