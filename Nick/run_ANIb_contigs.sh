#!/bin/sh -l

#$ -o run_ANIbc.out
#$ -e run_ANIbc.err
#$ -N run_ANIbc
#$ -cwd
#$ -q all.q

#Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp config_template.sh config.sh
fi
. ./config.sh
module load pyani/0.2.7
module load mashtree/0.29
module load ncbi-blast+/2.6.0
# ${mod_changers}/list_modules.sh

#
# Script to calculate the average nucleotide identity of a sample to numerous other samples from the same genus (genus dependent)
# The most similar match is identified and provided for confirmation
#
# Usage ./run_ANIb_contigs.sh sample_name   DB(for looking up reference, just relative path, also is genus)   Species   run_id  list_samples_to_include(optional)
#
# Python/3.5.2
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
elif [[ -z "${1}" ]]; then
	echo "Empty sample name supplied to $0, exiting"
	exit 1
# Gives the user a brief usage and help section if requested with the -h option argument
elif [[ "${1}" = "-h" ]]; then
	echo "Usage is ./run_ANI.sh sample_name ani_database(which is also genus) species run_id list_of_samples_to_include(optional)"
	echo "Output is saved to in ${processed}/sample_name/ANI"
	exit 0
elif [ -z "$2" ]; then
	echo "Empty database name supplied to $0. Second argument should be a genus found in ${local_DBs}/ANI/  ...Exiting"
	exit 1
elif [ ! -s "${local_DBs}/aniDB/${2,}" ]; then
	echo "The genus does not exist in the ANI database. This will be noted and the curator of the database will be notified. However, since nothing can be done at the moment....exiting"
	# Create a dummy folder to put non-results into (if it doesnt exist
	if [ ! -d "${processed}/${4}/${1}/ANI" ]; then  #create outdir if absent
		echo "${processed}/${4}/${1}/ANI"
		mkdir -p "${processed}/${4}/${1}/ANI"
	fi
	# Write non-results to a file in ANI folder
	echo "No matching ANI database found for ${2}(genus)" >> "${processed}/${4}/${1}/ANI/best_ANIb_contigs_hits_ordered(${1}_vs_${2}).txt"
	# Add genus to list to download and to database
	global_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")
	echo "ANI: ${2} - Found as ${1} on ${global_time}" >> "${shareScript}/maintenance_To_Do.txt"
	exit 1
elif [ -z "$3" ]; then
	echo "Empty species name supplied to $0. Third argument should be the suspected species of the sample. Exiting"
	exit 1
elif [ -z "$4" ]; then
	echo "Empty miseq_run_id name supplied to $0. Fourth argument should be the run id. Exiting"
	exit 1
elif [ ! -z "$5" ]; then
	others="true"
fi

start_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")
echo "Started ANI at ${start_time}"

# Sets the output folder to the sample_name folder in processed samples
OUTDATADIR="${processed}/${4}/${1}"

# Checks to see if an ANI folder already exists and creates it if not
if [ ! -d "$OUTDATADIR/ANI" ]; then
	echo "Creating $OUTDATADIR/ANI"
	mkdir -p "$OUTDATADIR/ANI"
fi

# Checks to see if the local DB ANI folder already exists and creates it if not. This is used to store a local copy of all samples in DB to be compared to (as to not disturb the source DB)
if [ ! -d "$OUTDATADIR/Contig_check/ANI/localANIDB_aniB_contigs" ]; then  #create outdir if absent
	echo "Creating $OUTDATADIR/Contig_check/ANI/localANIDB_aniB_contigs"
	mkdir -p "$OUTDATADIR/Contig_check/ANI/localANIDB_aniB_contigs"
else
	rm -r "$OUTDATADIR/Contig_check/ANI/localANIDB_aniB_contigs"
	mkdir -p "$OUTDATADIR/Contig_check/ANI/localANIDB_aniB_contigs"
fi

# Checks to see if the local DB ANI folder already exists and creates it if not. This is used to store a local copy of all samples in DB to be compared to (as to not disturb the source DB)
if [ ! -d "$OUTDATADIR/Contig_check/ANI/temp" ]; then  #create outdir if absent
	echo "Creating $OUTDATADIR/Contig_check/ANI/temp"
	mkdir -p "$OUTDATADIR/Contig_check/ANI/temp"
else
	rm -r "$OUTDATADIR/Contig_check/ANI/temp"
	mkdir -p "$OUTDATADIR/Contig_check/ANI/temp"
fi

# Gets persons name to use as email during entrez request to identify best matching sample
me=$(whoami)
#echo ${me}"___"${1}___${2}___${3}___${4}

# Sets the genus as the database that was passed in (The $2 seemed to be getting changed somewhere, so I just set it as a local variable)
genus_in=${2}

#Creates a local copy of the database folder
echo "trying to copy ${local_DBs}/aniDB/${genus_in,}/"
cp "${local_DBs}/aniDB/${genus_in,}/"*".fna" "${OUTDATADIR}/Contig_check/ANI/localANIDB_aniB_contigs/"
gunzip ${OUTDATADIR}/Contig_check/ANI/localANIDB_aniB_contigs/*.gz

#Copies the samples assembly contigs to the local ANI db folder
cp "${OUTDATADIR}/Contig_check/${1}_contigs_trimmed.fasta" "${OUTDATADIR}/Contig_check/ANI/localANIDB_aniB_contigs/sample_${2}_${3}.fasta"


# Add in all other assemblies to compare using list provided as argument
if [[ "${others}" = "true" ]]; then
	if [[ -f "${5}" ]]; then
		while IFS= read -r var; do
			sample_name=$(echo "${var}" | cut -d'/' -f2 | tr -d '[:space:]')
			project=$(echo "${var}" | cut -d'/' -f1 | tr -d '[:space:]')
			echo "${sample_name}, ${project}"
			if [[ "${project}" == "${4}" ]] && [[ "${sample_name}" == "${1}" ]]; then
				echo "Already in there as ref sample"
			else
				cp "${processed}/${project}/${sample_name}/Assembly/${sample_name}_contigs_trimmed.fasta" "${OUTDATADIR}/Contig_check/ANI/localANIDB_aniB_contigs/temp"
			fi
		done < ${5}
	else
		echo "List file: ${5}, does not exist. Analysis will continue with only database samples"
	fi
else
	echo "Analysis will be completed using only database isolates"
fi

#Renames all files in the localANIDB_aniB_contigs folder by changing extension from fna to fasta (which pyani needs)
for file in ${OUTDATADIR}/Contig_check/ANI/localANIDB_aniB_contigs/*.fna;
do
	fasta_name=$(basename "${file}" .fna)".fasta"
	mv "${file}" "${OUTDATADIR}/Contig_check/ANI/localANIDB_aniB_contigs/${fasta_name}"
done

# Mashtree trimming to reduce run time for ANI
owd=$(pwd)
cd ${OUTDATADIR}/Contig_check/ANI/localANIDB_aniB_contigs/
. "${mod_changers}/perl_5221_to_5161mt.sh"
mashtree.pl --numcpus ${procs} *.fasta --tempdir ${OUTDATADIR}/Contig_check/ANI/temp > ${OUTDATADIR}/Contig_check/ANI/"${genus_in}_and_${1}_mashtree_ANIb_contigs.dnd";
. "${mod_changers}/perl_5161mt_to_5221.sh"

# Get total number of isolates compared in tree
sample_count=$(find ${OUTDATADIR}/Contig_check/ANI/localANIDB_aniB_contigs/ -type f | wc -l)
# Must remove sample of interest
sample_count=$(( sample_count - 1 ))
# Check if sample count is greater than the max samples for tree size, if so then reduce tree size to max closest samples balanced around submitted isolate
if [[ ${sample_count} -gt ${max_ani_samples} ]]; then
	sleep 2
	tree=$(head -n 1 "${OUTDATADIR}/Contig_check/ANI/${genus_in}_and_${1}_mashtree_ANIb_contigs.dnd")
	echo $tree
	tree=$(echo "${tree}" | tr -d '()')
	echo $tree
	IFS=',' read -r -a samples <<< "${tree}"
	counter=0
	half_max=$(( (max_ani_samples+1) / 2 ))
	echo "Halfsies = ${half_max}"
	for sample in ${samples[@]};
	do
		counter=$(( counter + 1 ))
		filename=$(echo ${sample} | cut -d':' -f1)
		echo "${filename}"
		filename="${filename}.fasta"
		if [[ "${filename}" == "sample_${2}_${3}.fasta" ]]; then
			match=${counter}
			#echo "Match @ ${counter} and half=${half_max}"
			if [[ ${match} -le ${half_max} ]]; then
				#echo "LE"
				samples_trimmed=${samples[@]:0:$(( max_ani_samples + 1 ))}
			elif [[ ${match} -ge $(( sample_count - half_max)) ]]; then
				#echo "GE"
				samples_trimmed=${samples[@]:$(( max_ani_samples * -1 - 1 ))}
			else
				#echo "MID - $(( match - half_max )) to $(( counter + half_max + 1))"
				samples_trimmed=${samples[@]:$(( match - half_max )):${max_ani_samples}}
			fi
				#echo "${#samples_trimmed[@]}-${samples_trimmed[@]}"
				break
		fi
				#echo ${filename}
	done
	mkdir "${OUTDATADIR}/Contig_check/ANI/localANIDB_aniB_contigs_trimmed"
	for sample in ${samples_trimmed[@]};
	do
		filename=$(echo ${sample} | cut -d':' -f1)
		filename="${filename}.fasta"
		echo "Moving ${filename}"
		cp ${OUTDATADIR}/Contig_check/ANI/localANIDB_aniB_contigs/${filename} ${OUTDATADIR}/Contig_check/ANI/localANIDB_aniB_contigs_trimmed/
	done
	rm -r "${OUTDATADIR}/Contig_check/ANI/localANIDB_aniB_contigs_full"
	mv "${OUTDATADIR}/Contig_check/ANI/localANIDB_aniB_contigs" "${OUTDATADIR}/Contig_check/ANI/localANIDB_aniB_contigs_full"
	mv "${OUTDATADIR}/Contig_check/ANI/localANIDB_aniB_contigs_trimmed" "${OUTDATADIR}/Contig_check/ANI/localANIDB_aniB_contigs"
# Continue without reducing the tree, as there are not enough samples to require reduction
else
	echo "Sample count below limit, not trimming ANI database"
fi

cd ${owd}
# Resume normal ANI analysis after mashtree reduction

# Checks for a previous copy of the aniM folder, removes it if found
if [ -d "${OUTDATADIR}/Contig_check/ANI/aniB" ]; then
	echo "Removing old aniB results in ${OUTDATADIR}/Contig_check/ANI/aniB"
	rm -r "${OUTDATADIR}/Contig_check/ANI/aniB"
fi

#Calls pyani on local db folder
python -V
#python "${shareScript}/pyani/average_nucleotide_identity.py" -i "${OUTDATADIR}/Contig_check/ANI/localANIDB_aniB_contigs" -o "${OUTDATADIR}/Contig_check/ANI/aniB_contigs" --write_excel
average_nucleotide_identity.py -i "${OUTDATADIR}/Contig_check/ANI/localANIDB_aniB_contigs" -o "${OUTDATADIR}/Contig_check/ANI/aniB_contigs" -m ANIb --write_excel

#Extracts the query sample info line for percentage identity from the percent identity file
while IFS='' read -r line;
do
#	echo "!-${line}"
	if [[ ${line:0:7} = "sample_" ]]; then
		sampleline=${line}
#		echo "found it-"$sampleline
		break
	fi
done < "${OUTDATADIR}/Contig_check/ANI/aniB_contigs/ANIb_percentage_identity.tab"

#Extracts the top line from the %id file to get all the sample names used in analysis (they are tab separated along the top row)
if [[ -s "${OUTDATADIR}/Contig_check/ANI/aniB_contigs/ANIb_percentage_identity.tab" ]]; then
	firstline=$(head -n 1 "${OUTDATADIR}/Contig_check/ANI/aniB_contigs/ANIb_percentage_identity.tab")
else
	echo "No "${OUTDATADIR}/Contig_check/ANI/aniB_contigs/ANIb_percentage_identity.tab" file, exiting"
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
	definition=$(head -1 "${OUTDATADIR}/Contig_check/ANI/localANIDB_aniB_contigs/${samples[i]}.fasta")
	# Prints all matching samples to file (Except the self comparison) by line as percent_match  sample_name  fasta_header
	echo "${percents[i+1]} ${samples[i]} ${definition}" >> "${OUTDATADIR}/Contig_check/ANI/best_hits.txt"
done

#Sorts the list in the file based on %id (best to worst)
sort -nr -t' ' -k1 -o "${OUTDATADIR}/Contig_check/ANI/best_hits_ordered.txt" "${OUTDATADIR}/Contig_check/ANI/best_hits.txt"
#Extracts the first line of the file (best hit)
best=$(head -n 1 "${OUTDATADIR}/Contig_check/ANI/best_hits_ordered.txt")
#Creates an array from the best hit
IFS=' ' read -r -a def_array <<< "${best}"
#echo -${def_array[@]}+
#Captures the assembly file name that the best hit came from
best_file=${def_array[1]}
#Formats the %id to standard percentage (xx.xx%)
best_percent=$(awk -v per="${def_array[0]}" 'BEGIN{printf "%.2f", per * 100}')
#echo "${best_file}"
# If the best match comes from the additional file, extract the taxonomy from that file
if [[ "${best_file}" = *"_contigs_trimmed" ]]; then
	best_outbreak_match=$(echo "${best_file}" | rev | cut -d'_' -f3- | rev)
	while IFS= read -r var
	do
		sample_name=$(echo "${var}" | cut -d'/' -f2 | tr -d '[:space:]')
		if [[ "${sample_name}" = "${best_outbreak_match}" ]]; then
			project=$(echo "${var}" | cut -d'/' -f1 | tr -d '[:space:]')
			while IFS= read -r pstats_line
				do
					tool=$(echo "${pstats_line}" | cut -d':' -f1 | tr -s " ")
					#echo ":${tool}:"
					if [[ "${tool}" = "weighted Classify " ]]; then
						best_organism_guess=$(echo "${pstats_line}" | cut -d':' -f3 | cut -d' ' -f3,4)
						break 2
					fi
				done < ${processed}/${project}/${sample_name}/${sample_name}_pipeline_stats.txt
		fi
	done < ${5}
# if the best hit cmoes from the aniDB then lookup the taxonomy on ncbi
else
	#Extracts the accession number from the definition line
	accession=$(echo "${def_array[2]}" | cut -d' ' -f1  | cut -d'>' -f2)
	#Looks up the NCBI genus species from the accession number
	if [[ "${accession}" == "No_Accession_Number" ]]; then
		best_organism_guess="${def_array[3]} ${def_array[4]}"
	else
		best_organism_guess=$(python "${shareScript}/entrez_get_taxon_from_accession.py" "${accession}" "${me}")
	fi
fi
# Uncomment this if you want to restrict ID to only genus species, without more resolute definition
#best_organism_guess_arr=($best_organism_guess})
#best_organism_guess="${best_organism_guess_arr[@]:0:2}"

#Creates a line at the top of the file to show the best match in an easily readable format that matches the style on the MMB_Seq log
echo -e "${best_percent}%-${best_organism_guess}(${best_file}.fna)\\n$(cat "${OUTDATADIR}/Contig_check/ANI/best_hits_ordered.txt")" > "${OUTDATADIR}/Contig_check/ANI/best_ANIb_contigs_hits_ordered(${1}_vs_${genus_in}).txt"

#Removes the transient hit files
if [ -s "${OUTDATADIR}/Contig_check/ANI/best_hits.txt" ]; then
	rm "${OUTDATADIR}/Contig_check/ANI/best_hits.txt"
#	echo "1"
fi
if [ -s "${OUTDATADIR}/Contig_check/ANI/best_hits_ordered.txt" ]; then
	rm "${OUTDATADIR}/Contig_check/ANI/best_hits_ordered.txt"
#	echo "2"
fi

end_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")
echo "ENDed ANI at ${end_time}"

#Script exited gracefully (unless something else inside failed)
exit 0
