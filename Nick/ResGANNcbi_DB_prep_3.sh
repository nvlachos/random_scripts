#!/bin/sh -l

#$ -o prepCSSTARDB_3.out
#$ -e prepCSSTARDB_3.err
#$ -N prepCSSTARDB_3
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
# Usage ./ResGANNCBI_DB_prep_3.sh path_to_dir
#

$(python2 -V)


today=$(date '+%Y%m%d')

DATADIR="${local_DBs}/star/db_prep"
ResGANNCBI_source="${DATADIR}/ResGANNCBI_${today}.fasta"

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, using default DATADIR=${local_DBs}/star/db_prep"
	#if [ -d "${DATADIR}" ]; then
	#	rm -r "${DATADIR}"
	#fi
	#mkdir "${DATADIR}"
else
	if [[ -z "${1}" ]]; then
		echo "Empty path supplied to $0, exiting"
		exit 1
	# Gives the user a brief usage and help section if requested with the -h option argument
	elif [[ "${1}" = "-h" ]]; then
		echo "Usage is ./ResGANNCBI_DB_prep_3.sh "
		echo "Converts the ARGANNOT fasta and resFinder ZIP files into a single fasta file to a single srst2 db fasta for use with csstar and srst2"
		echo "Output is ResGANNCBI_date_srst2.fasta"
		exit 0
	elif [[ "${1}" = "-f" ]]; then
		if [[ -z "${2}" ]]; then
			echo "Empty filename supplied to $0, exiting"
			exit 1
		else
			sourced="true"
			ResGANNCBI_source="${DATADIR}/${2}"
		fi
	else
		DATADIR="${1}"
	fi
fi
date


#Creates an associative array for looking up the genes to what they confer
declare -A groups
echo "Creating reference array"
while IFS= read -r line  || [ -n "$line" ]; do
	line=${line,,}
	gene=$(echo "${line}" | cut -d ':' -f1)
	first=${gene:0:1}
	if [ "$first" == "#" ]; then
		continue
	fi
	confers=$(echo "${line}" | cut -d ':' -f2)
	groups[${gene}]="${confers}"
done < "${local_DBs}/star/group_defs.txt"



ml CD-HIT/4.6
# #switch to python 2.7 temporarily for cdhit to csv scripts/config
cp "${local_DBs}/star/cdhit_to_csv.py" "${DATADIR}"
cp "${local_DBs}/star/csv_to_gene_db.py" "${DATADIR}"

cd "${DATADIR}"

echo "--- About to CD-HIT-EST on ResGANNCBI DB-${ResGANNCBI_source} ---"
cd-hit-est "-i" "${ResGANNCBI_source}" "-o" "${DATADIR}/ResGANNCBI_${today}_cdhit90" "-d" "0" > "${DATADIR}/ResGANNCBI_${today}_cdhit90.stdout"
# cd-hit-est -i /scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/DBs/star/db_prep/Combined.fasta -o /scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/DBs/star/db_prep/ResGANNCBI_040518_cdhit90 -d 0 > /scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/DBs/star/db_prep/ResGANNCBI_040518_cdhit90.stdout

echo "--- About to CD-HIT-to-CSV on ResGANNCBI DB-${ResGANNCBI_source} ---"

python "${DATADIR}/cdhit_to_csv.py" "--cluster_file" "${DATADIR}/ResGANNCBI_${today}_cdhit90.clstr" "--infasta" "${ResGANNCBI_source}" "--outfile" "${DATADIR}/ResGANNCBI_${today}_clustered.csv"
# python /scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/DBs/star/db_prep/cdhit_to_csv.py --cluster_file /scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/DBs/star/db_prep/ResGANNCBI_040518_cdhit90.clstr --infasta /scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/DBs/star/db_prep/Combined.fasta --outfile /scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/DBs/star/db_prep/ResGANNCBI_040518_clustered.csv

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
echo "Creating ResGANNCBI seq reference array"
while IFS= read -r line  || [ -n "$line" ]; do
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
done < "${ResGANNCBI_source}"
#add the last sequence to the array, since there is no > to trigger its addition
#echo "${header}:${seq}" >> ${DATADIR}/hits.txt
echo "Adding final item, number ${count}"
seqarr[${header}]="${seq}"
echo "There are ${#seqarr[@]} different entries, and counter equals ${count}"

if [ -f ${DATADIR}/ResGANNCBI_${today}_clustered.csv ]; then
	echo "found ${DATADIR}/ResGANNCBI_${today}_clustered.csv"
fi




counter=0
echo "Before"
while IFS= read -r gene_line || [ -n "$gene_line" ]; do
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
			echo "allele-${allele}:accession-${accession}"
			allele+="_${accession}"
			allele_location=$(echo "${line_items[3]}" | cut -d':' -f3)
			allele_length=$(echo "${line_items[3]}" | cut -d':' -f4)
			group="${groups[${group_raw,,}]}"
			DB_ID="${line_items[3]}"
			#echo "looking up sequence:${line_items[3]}:"
			dnaseq="${seqarr[${line_items[3]}]}"
			#echo "result:${dnaseq}:"
		elif [[ "${source}" == "NCBI" ]]; then
		# Create a new way to identify other DB entries
			seqID=${line_items[0]}
			clusterID=${line_items[1]}
			gene=${line_items[2]}
			from_source="NCBI"
			if [[ "${gene}" == *":"* ]]; then
				gene=$(echo ${gene} | cut -d':' -f1)
			fi
			for (( i=0; i<${#gene}; i++ ));
			#do
			#	if [[ "${gene:${i}:1}" == ")" ]]; then
			#		i=$(( i + 1 ))
			#		gene=${gene:${i}}
			#		break
			#	fi
			#done
			allele=$(echo "${line_items[3]}" | cut -d':' -f1 | cut -d']' -f2-)
			group_raw=$(echo "${line_items[2]}" | cut -d')' -f1 | cut -d'(' -f2)
			accession=$(echo "${line_items[3]}" | cut -d':' -f3)
			echo "allele-${allele}:accession-${accession}"
			allele+="_${accession}"
			allele_location=$(echo "${line_items[3]}" | cut -d':' -f4)
			allele_start=$(echo "${allele_location}" | cut -d'-' -f1)
			allele_end=$(echo "${allele_location}" | cut -d'-' -f2)
			allele_length=$(( allele end - allele_start ))
			group="${groups[${group_raw,,}]}"
			DB_ID="${line_items[3]}"
			#echo "looking up sequence:${line_items[3]}:"
			dnaseq="${seqarr[${line_items[3]}]}"
			#echo "result:${dnaseq}:"
		else
			seqID=${line_items[0]}
			clusterID=${line_items[1]}
			gene=$(echo "${line_items[2]}" | cut -d']' -f2 | cut -d'_' -f1)
			if [[ "${source}" == "RES" ]]; then
				from_source="ResFinder"
			else
				from_source="${source}"
			fi
			allele=$(echo ${line_items[3]} | cut -d'_' -f1,2 | cut -d']' -f2)
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
	echo -e "${seqID},${clusterID},${gene},${allele},${dnaseq},${DB_ID},${accession},${group},${from_source}" >> "${DATADIR}/ResGANNCBI_${today}_cluster90_edited.csv"
	counter=$(( counter + 1 ))
done < ${DATADIR}/ResGANNCBI_${today}_clustered.csv
echo "After"

cp "${DATADIR}/ResGANNCBI_${today}_cluster90_edited.csv" "/scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/DBs/star/"

echo "--- About to csv to gene ---"
python "${DATADIR}/csv_to_gene_db.py" "-t" "${DATADIR}/ResGANNCBI_${today}_cluster90_edited.csv" "-o" "${DATADIR}/ResGANNCBI_${today}_srst2.fasta" "-f" "${ResGANNCBI_source}" "-c" "6" "-s" "5"

#Additional Enumerator, specfically to address fos gene ambiguity
# declare -A allele_tally
# while IFS= read -r line  || [ -n "$line" ]; do
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
# done < "${DATADIR}/ResGANNCBI_${today}_srst2.fasta"






exit
