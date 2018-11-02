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

local_DBs="/scicomp/groups/OID/NCEZID/DHQP/CEMB/databases"

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
if [[ ! -d ${OUTDATADIR}/all_test/dists ]]; then
	mkdir -p ${OUTDATADIR}/all_test/dists
fi

# Gets persons name to use as email during entrez request to identify best matching sample
me=$(whoami)
#echo ${me}"___"${1}___${2}___${3}___${4}

#Creates a local copy of the database folder
# cp ${share}/DBs/aniDB/all/compound_sketch_all.msh "${OUTDATADIR}/ANI/"

#mash dist "${local_DBs}/aniDB/refseq.genomes.k21s1000.msh" "${OUTDATADIR}/Assembly/${1}_scaffolds_trimmed.fasta" > "${OUTDATADIR}/ANI/${1}_all_refSeq.dists"
counter=0
#for ref in ${local_DBs}/aniDB/all/*.fna; do
#	echo ${ref}
#	counter=$(( counter + 1 ))
#	filename=$(basename ${ref} | cut -d'_' -f1,2)
#	mash dist "${local_DBs}/aniDB/all/all_sketch.msh" "${ref}" > "${OUTDATADIR}/all_test/dists/${filename}_unsorted.dists"
#	sort -k3 -n -o "${OUTDATADIR}/all_test/dists/${filename}.dists" "${OUTDATADIR}/all_test/dists/${filename}_unsorted.dists"
#	rm -r "${OUTDATADIR}/all_test/dists/${filename}_unsorted.dists"
#done

for distfile in ${local_DBs}/aniDB/all_test/dists/*.dists; do
	[ -f "$distfile" ] || break
	taxa=$(basename ${distfile})
	if [[ ! -d ${local_DBs}/aniDB/all_test/${taxa} ]]; then
		mkdir -p ${local_DBs}/aniDB/all_test/${taxa}/localANIDB
	fi
	counter=0
	max_ani_samples=30
	echo "${distfile}-${taxa}"
	> "${local_DBs}/aniDB/all_test/${taxa}/thirty_closest_dists.txt"
	if [[ ! -d ${local_DBs}/aniDB/all_test/${taxa}/localANIDB ]]; then
		mkdir "${local_DBs}/aniDB/all_test/${taxa}/localANIDB"
	fi
	while IFS= read -r line;  do
			echo "${counter}:-:-:${line}"
			if [[ ${counter} -eq 0 ]]; then
				ref_path=$(echo "${line}" | cut -d'	' -f2)
				echo "rp-${ref_path}"
				echo "${ref_path}" >> "${local_DBs}/aniDB/all_test/${taxa}/thirty_closest_dists.txt"
				fasta=$(basename ${ref_path})
				fasta=${fasta:0:-3}"fasta"
				echo "cp ${ref_path} ${local_DBs}/aniDB/all_test/${taxa}/localANIDB/${fasta}"
				cp ${ref_path} ${local_DBs}/aniDB/all_test/${taxa}/localANIDB/${fasta}
			fi
			if [[ ${counter} -gt ${max_ani_samples} ]]; then
				break
			else
				source_path=$(echo "${line}" | cut -d'	' -f1)
				echo "sp-${source_path}"
				echo "${source_path}" >> "${local_DBs}/aniDB/all_test/${taxa}/thirty_closest_dists.txt"
				fasta=$(basename ${source_path})
				fasta=${fasta:0:-3}"fasta"
				echo "cp ${source_path} ${local_DBs}/aniDB/all_test/${taxa}/localANIDB/${fasta}"
				cp ${source_path} ${local_DBs}/aniDB/all_test/${taxa}/localANIDB/${fasta}
			fi
			counter=$(( counter + 1 ))
	done < ${distfile}
done
echo ${counter}

exit

sub_counter=0
max_subs=50
samples=()
main_dir="/scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/mass_subs/ani_TEST"
mkdir ${main_dir}

for localANIDB in ${local_DBs}/aniDB/Single_ANI_Test/*; do
	echo "${localANIDB}"
	sample=$(basename ${localANIDB} | rev | cut -d'/' -f1 | rev)
	if [[ ! -d ${localANIDB} ]]; then
		break
	else
		samples[sub_counter]=${sample}
	fi
	echo "${sample}"
	if [[ ${sub_counter} -lt ${max_subs} ]]; then
		echo  "Index is below max submissions, submitting"
		echo -e "#!/bin/bash -l\n" > "${main_dir}/aniB_${sample}_${start_time}.sh"
		echo -e "#$ -o aniB_${sample}.out" >> "${main_dir}/aniB_${sample}_${start_time}.sh"
		echo -e "#$ -e aniB_${sample}.err" >> "${main_dir}/aniB_${sample}_${start_time}.sh"
		echo -e "#$ -N aniB_${sample}"   >> "${main_dir}/aniB_${sample}_${start_time}.sh"
		echo -e "#$ -cwd"  >> "${main_dir}/aniB_${sample}_${start_time}.sh"
		echo -e "#$ -q short.q\n"  >> "${main_dir}/aniB_${sample}_${start_time}.sh"
		echo -e "module load pyani/1.0" "${main_dir}/aniB_${sample}_${start_time}.sh"
		echo -e "average_nucleotide_identity.py -i \"${localANIDB}/localANIDB\" -o \"${local_DBs}/aniDB/all_test/${taxa}/aniB\" \"--write_excel\"" >> "${main_dir}/aniB_${sample}_${start_time}.sh"
		echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_aniB_complete.txt\"" >> "${main_dir}/aniB_${sample}_${start_time}.sh"
		qsub "${main_dir}/aniB_${sample}_${start_time}.sh"

		echo -e "#!/bin/bash -l\n" > "${main_dir}/aniM_${sample}_${start_time}.sh"
		echo -e "#$ -o aniM_${sample}.out" >> "${main_dir}/aniM_${sample}_${start_time}.sh"
		echo -e "#$ -e aniM_${sample}.err" >> "${main_dir}/aniM_${sample}_${start_time}.sh"
		echo -e "#$ -N aniM_${sample}"   >> "${main_dir}/aniM_${sample}_${start_time}.sh"
		echo -e "#$ -cwd"  >> "${main_dir}/aniM_${sample}_${start_time}.sh"
		echo -e "#$ -q short.q\n"  >> "${main_dir}/aniM_${sample}_${start_time}.sh"
		echo -e "module load pyani/1.0" >> "${main_dir}/aniM_${sample}_${start_time}.sh"
		echo -e "average_nucleotide_identity.py -i \"${localANIDB}/localANIDB\" -o \"${local_DBs}/aniDB/all_test/${taxa}/aniM\" \"--write_excel\"" >> "${main_dir}/aniM_${sample}_${start_time}.sh"
		echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_aniM_complete.txt\"" >> "${main_dir}/aniM_${sample}_${start_time}.sh"
		qsub "${main_dir}/aniM_${sample}_${start_time}.sh"
else
	waiting_for_index=$(( counter - max_subs ))
	waiting_sample=$(echo "${samples[${waiting_for_index}]}")
	timer=0
	echo "Index is above max submissions, waiting for index ${waiting_for_index}:${waiting_sample} to complete"
	while :
	do
		if [[ ${timer} -gt 1800 ]]; then
			echo "Timer exceeded limit of 1800 seconds 30 minutes"
			break
		fi
		if [[ -f "${main_dir}/complete/${waiting_sample}_aniB_complete.txt" ]] && [[ -f "${main_dir}/complete/${waiting_sample}_aniM_complete.txt" ]]; then
			echo  "Index is below max submissions, submitting"
			echo -e "#!/bin/bash -l\n" > "${main_dir}/ani_${sample}_${start_time}.sh"
			echo -e "#$ -o aniB_${sample}.out" >> "${main_dir}/aniB_${sample}_${start_time}.sh"
			echo -e "#$ -e aniB_${sample}.err" >> "${main_dir}/aniB_${sample}_${start_time}.sh"
			echo -e "#$ -N aniB_${sample}"   >> "${main_dir}/aniB_${sample}_${start_time}.sh"
			echo -e "#$ -cwd" >> "${main_dir}/aniB_${sample}_${start_time}.sh"
			echo -e "#$ -q short.q\n"  >> "${main_dir}/aniB_${sample}_${start_time}.sh"
			echo -e "module load pyani/1.0" "${main_dir}/aniB_${sample}_${start_time}.sh"
			echo -e "average_nucleotide_identity.py -i \"${localANIDB}/localANIDB\" -o \"${local_DBs}/aniDB/all_test/${taxa}/aniB\" \"--write_excel\"" >> "${main_dir}/aniB_${sample}_${start_time}.sh"
			echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_aniB_complete.txt\"" >> "${main_dir}/aniB_${sample}_${start_time}.sh"
			qsub "${main_dir}/aniB_${sample}_${start_time}.sh"

			echo -e "#!/bin/bash -l\n" > "${main_dir}/aniM_${sample}_${start_time}.sh"
			echo -e "#$ -o aniM_${sample}.out" >> "${main_dir}/aniM_${sample}_${start_time}.sh"
			echo -e "#$ -e aniM_${sample}.err" >> "${main_dir}/aniM_${sample}_${start_time}.sh"
			echo -e "#$ -N aniM_${sample}"   >> "${main_dir}/aniM_${sample}_${start_time}.sh"
			echo -e "#$ -cwd"  >> "${main_dir}/aniM_${sample}_${start_time}.sh"
			echo -e "#$ -q short.q\n"  >> "${main_dir}/aniM_${sample}_${start_time}.sh"
			echo -e "module load pyani/1.0" >> "${main_dir}/aniM_${sample}_${start_time}.sh"
			echo -e "average_nucleotide_identity.py -i \"${localANIDB}/localANIDB\" -o \"${local_DBs}/aniDB/all_test/${taxa}/aniM\" \"--write_excel\"" >> "${main_dir}/aniM_${sample}_${start_time}.sh"
			echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_aniM_complete.txt\"" >> "${main_dir}/aniM_${sample}_${start_time}.sh"
			qsub "${main_dir}/aniM_${sample}_${start_time}.sh"
	fi
	done
fi
	sub_counter=$(( counter + 1 ))
done





#Extracts the query sample info line for percentage identity from the percent identity file
while IFS='' read -r line;
do
#	echo "!-${line}"
	if [[ ${line:0:7} = "sample_" ]]; then
		sampleline=${line}
#		echo "found it-"$sampleline
		break
	fi
done < "${local_DBs}/aniDB/all_test/${taxa}/aniM/ANIm_percentage_identity.tab"

#Extracts the top line from the %id file to get all the sample names used in analysis (they are tab separated along the top row)
if [[ -s "${local_DBs}/aniDB/all_test/${taxa}/aniM/ANIm_percentage_identity.tab" ]]; then
	firstline=$(head -n 1 "${local_DBs}/aniDB/all_test/${taxa}/aniM/ANIm_percentage_identity.tab")
else
	echo "No "${local_DBs}/aniDB/all_test/${taxa}/aniM/ANIm_percentage_identity.tab" file, exiting"
	exit 1
fi

#Extracts the query sample info line for percentage identity from the percent identity file
while IFS='' read -r line;
do
#	echo "!-${line}"
	if [[ ${line:0:7} = "sample_" ]]; then
		sampleline=${line}
#		echo "found it-"$sampleline
		break
	fi
done < "${local_DBs}/aniDB/all_test/${taxa}/aniM/ANIm_alignment_coverage.tab"

#Extracts the top line from the %id file to get all the sample names used in analysis (they are tab separated along the top row)
if [[ -s "${local_DBs}/aniDB/all_test/${taxa}/aniM/ANIm_alignment_coverage.tab" ]]; then
	firstline=$(head -n 1 "${local_DBs}/aniDB/all_test/${taxa}/aniM/ANIm_alignment_coverage.tab")
else
	echo "No "${local_DBs}/aniDB/all_test/${taxa}/aniM/ANIm_alignment_coverage.tab" file, exiting"
	exit 1
fi

#Arrays to read sample names and the %ids for the query sample against those other samples
IFS="	" read -r -a samples_aniM_coverage <<< "${firstline}"
IFS="	" read -r -a percents_aniM_coverage <<< "${sampleline}"
declare -A coverage_array
counter=0
for isolate in "${firstline[@]}"; do
	coverage_array[isolate]=${sampleline[${counter}]}
	counter=$(( counter + 1 ))
done


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
	definition=$(head -1 "${local_DBs}/aniDB/all_test/${taxa}/localANIDB/${samples[i]}.fna")
	# Prints all matching samples to file (Except the self comparison) by line as percent_match  sample_name  fasta_header
	echo "${percents[i+1]} ${coverage_array[${samples[i]}]}	${samples[i]} ${definition}" >> "${local_DBs}/aniDB/all_test/${taxa}/best_hits_aniM.txt"
done

#Sorts the list in the file based on %id (best to worst)
sort -nr -t' ' -k1 -o "${local_DBs}/aniDB/all_test/${taxa}/best_hits__aniM_ordered.txt" "${local_DBs}/aniDB/all_test/${taxa}/best_hits_aniM.txt"

#Extracts the query sample info line for percentage identity from the percent identity file
cov_counter=0
total_lines
while IFS='' read -r line;
do
	coverage=$(echo ${line} | cut -d '	' -f2)
	if [[ coverage -eq 1 ]]; then
		coverage=100
	elif [[ coverage -eq 0 ]]; then
		coverage=0
	else
		coverage=$(echo ${coverage} | cut '.' -f2)
		coverage=${coverage:0:2}
	fi
	if [[ ${coverage} -gt 40]]; then
		counter=$(( counter + 1))
		total_lines=$(( total_lines + 1 ))
	else
		total_lines=$(( total_lines + 1 ))
	fi
done < "${local_DBs}/aniDB/all_test/${taxa}/best_hits_aniM.txt"

high_cov_limit=$(( total_lines - cov_counter ))
low_cov=$(head -n${cov_counter} "${local_DBs}/aniDB/all_test/${taxa}/best_hits_aniM.txt")
high_cov=$(tail -n${high_cov_limit} "${local_DBs}/aniDB/all_test/${taxa}/best_hits_aniM.txt")
mv "${local_DBs}/aniDB/all_test/${taxa}/best_hits_aniM.txt" "${local_DBs}/aniDB/all_test/${taxa}/best_hits_aniM_sorted.txt"
cat ${high_cov} ${low_cov} > "${local_DBs}/aniDB/all_test/${taxa}/best_hits_aniM_filtered.txt"




#Extracts the query sample info line for percentage identity from the percent identity file
while IFS='' read -r line;
do
#	echo "!-${line}"
	if [[ ${line:0:7} = "sample_" ]]; then
		sampleline=${line}
#		echo "found it-"$sampleline
		break
	fi
done < "${local_DBs}/aniDB/all_test/${taxa}/aniB/ANIb_percentage_identity.tab"

#Extracts the top line from the %id file to get all the sample names used in analysis (they are tab separated along the top row)
if [[ -s "${local_DBs}/aniDB/all_test/${taxa}/aniB/ANIb_percentage_identity.tab" ]]; then
	firstline=$(head -n 1 "${local_DBs}/aniDB/all_test/${taxa}/aniB/ANIb_percentage_identity.tab")
else
	echo "No "${local_DBs}/aniDB/all_test/${taxa}/aniB/ANIb_percentage_identity.tab" file, exiting"
	exit 1
fi

#Extracts the query sample info line for percentage identity from the percent identity file
while IFS='' read -r line;
do
#	echo "!-${line}"
	if [[ ${line:0:7} = "sample_" ]]; then
		sampleline=${line}
#		echo "found it-"$sampleline
		break
	fi
done < "${local_DBs}/aniDB/all_test/${taxa}/aniB/ANIb_alignment_coverage.tab"

#Extracts the top line from the %id file to get all the sample names used in analysis (they are tab separated along the top row)
if [[ -s "${local_DBs}/aniDB/all_test/${taxa}/aniB/ANIb_alignment_coverage.tab" ]]; then
	firstline=$(head -n 1 "${local_DBs}/aniDB/all_test/${taxa}/aniB/ANIb_alignment_coverage.tab")
else
	echo "No "${local_DBs}/aniDB/all_test/${taxa}/aniB/ANIb_alignment_coverage.tab" file, exiting"
	exit 1
fi

#Arrays to read sample names and the %ids for the query sample against those other samples
IFS="	" read -r -a samples_aniB_coverage <<< "${firstline}"
IFS="	" read -r -a percents_aniB_coverage <<< "${sampleline}"
declare -A coverage_array
counter=0
for isolate in "${firstline[@]}"; do
	coverage_array[isolate]=${sampleline[${counter}]}
	counter=$(( counter + 1 ))
done


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
	definition=$(head -1 "${local_DBs}/aniDB/all_test/${taxa}/localANIDB/${samples[i]}.fna")
	# Prints all matching samples to file (Except the self comparison) by line as percent_match  sample_name  fasta_header
	echo "${percents[i+1]} ${coverage_array[${samples[i]}]}	${samples[i]} ${definition}" >> "${local_DBs}/aniDB/all_test/${taxa}/best_hits_aniB.txt"
done

#Sorts the list in the file based on %id (best to worst)
sort -nr -t' ' -k1 -o "${local_DBs}/aniDB/all_test/${taxa}/best_hits__aniM_ordered.txt" "${local_DBs}/aniDB/all_test/${taxa}/best_hits_aniB.txt"

#Extracts the query sample info line for percentage identity from the percent identity file
cov_counter=0
total_lines
while IFS='' read -r line;
do
	coverage=$(echo ${line} | cut -d '	' -f2)
	if [[ coverage -eq 1 ]]; then
		coverage=100
	elif [[ coverage -eq 0 ]]; then
		coverage=0
	else
		coverage=$(echo ${coverage} | cut '.' -f2)
		coverage=${coverage:0:2}
	fi
	if [[ ${coverage} -gt 40]]; then
		counter=$(( counter + 1))
		total_lines=$(( total_lines + 1 ))
	else
		total_lines=$(( total_lines + 1 ))
	fi
done < "${local_DBs}/aniDB/all_test/${taxa}/best_hits_aniB.txt"

high_cov_limit=$(( total_lines - cov_counter ))
low_cov=$(head -n${cov_counter} "${local_DBs}/aniDB/all_test/${taxa}/best_hits_aniB.txt")
high_cov=$(tail -n${high_cov_limit} "${local_DBs}/aniDB/all_test/${taxa}/best_hits_aniB.txt")
mv "${local_DBs}/aniDB/all_test/${taxa}/best_hits_aniM.txt" "${local_DBs}/aniDB/all_test/${taxa}/best_hits_aniB_sorted.txt"
cat ${high_cov} ${low_cov} > "${local_DBs}/aniDB/all_test/${taxa}/best_hits_aniB_filtered.txt"

end_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")
echo "ENDed ANI at ${end_time}"

#Script exited gracefully (unless something else inside failed)
exit 0
