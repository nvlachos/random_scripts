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

working_dir="all_test"

# Sets the output folder to the sample_name folder in processed samples
OUTDATADIR="${local_DBs}/aniDB"

# counter=0
# for ref in ${local_DBs}/aniDB/all/*.fna; do
# 	echo ${ref}
# 	counter=$(( counter + 1 ))
# 	filename=$(basename ${ref} | cut -d'_' -f1,2)
# 	mash dist "${local_DBs}/aniDB/all/all_sketch.msh" "${ref}" > "${OUTDATADIR}/${working_dir}/dists/${filename}_unsorted.dists"
# 	sort -k3 -n -o "${OUTDATADIR}/${working_dir}/dists/${filename}.dists" "${OUTDATADIR}/${working_dir}/dists/${filename}_unsorted.dists"
# 	rm -r "${OUTDATADIR}/${working_dir}/dists/${filename}_unsorted.dists"
# done
#
# for distfile in ${local_DBs}/aniDB/${working_dir}/dists/*.dists; do
# 	[ -f "$distfile" ] || break
# 	taxa=$(basename ${distfile})
# 	if [[ ! -d ${local_DBs}/aniDB/${working_dir}/${taxa} ]]; then
# 		mkdir -p ${local_DBs}/aniDB/${working_dir}/${taxa}/localANIDB
# 	fi
# 	counter=0
# 	max_ani_samples=30
# 	echo "${distfile}-${taxa}"
# 	> "${local_DBs}/aniDB/${working_dir}/${taxa}/thirty_closest_dists.txt"
# 	if [[ ! -d ${local_DBs}/aniDB/${working_dir}/${taxa}/localANIDB ]]; then
# 		mkdir "${local_DBs}/aniDB/${working_dir}/${taxa}/localANIDB"
# 	fi
# 	while IFS= read -r line;  do
# 			echo "${counter}:-:-:${line}"
# 			if [[ ${counter} -eq 0 ]]; then
# 				ref_path=$(echo "${line}" | cut -d'	' -f2)
# 				echo "rp-${ref_path}"
# 				echo "${ref_path}" >> "${local_DBs}/aniDB/${working_dir}/${taxa}/thirty_closest_dists.txt"
# 				fasta=$(basename ${ref_path})
# 				fasta=${fasta:0:-3}"fasta"
# 				echo "cp ${ref_path} ${local_DBs}/aniDB/${working_dir}/${taxa}/localANIDB/${fasta}"
# 				cp ${ref_path} ${local_DBs}/aniDB/${working_dir}/${taxa}/localANIDB/${fasta}
# 			fi
# 			if [[ ${counter} -gt ${max_ani_samples} ]]; then
# 				break
# 			else
# 				source_path=$(echo "${line}" | cut -d'	' -f1)
# 				echo "sp-${source_path}"
# 				echo "${source_path}" >> "${local_DBs}/aniDB/${working_dir}/${taxa}/thirty_closest_dists.txt"
# 				fasta=$(basename ${source_path})
# 				fasta=${fasta:0:-3}"fasta"
# 				echo "cp ${source_path} ${local_DBs}/aniDB/${working_dir}/${taxa}/localANIDB/${fasta}"
# 				cp ${source_path} ${local_DBs}/aniDB/${working_dir}/${taxa}/localANIDB/${fasta}
# 			fi
# 			counter=$(( counter + 1 ))
# 	done < ${distfile}
# done
# echo ${counter}



sub_counter=0
max_subs=50
samples=()
main_dir="/scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/mass_subs/ani_TEST"
mkdir ${main_dir}


for ref_tax in ${local_DBs}/aniDB/${working_dir}/*; do
	echo "${ref_tax}"
	#sub_counter=$(( sub_counter + 1 ))
	#continue
	sample=$(basename ${ref_tax} | rev | cut -d'/' -f1 | rev)
	echo "sample: ${sample}"
	echo "Looking in ${ref_tax}/localANIDB for ${sample}*.fasta"
	reference=$(find ${ref_tax}/localANIDB -name "${sample}*.fasta" -type f)
	echo "reference: ${reference}"
	#mv "${reference}" "${ref_tax}"
	closest=$(wc -l ${ref_tax}/thirty_closest_dists.txt | cut -d ' ' -f1)
	if [[ ${closest} -gt 30 ]]; then
		echo "$(tail -30 ${ref_tax}/thirty_closest_dists.txt)" > ${ref_tax}/thirty_closest_dists_30.txt
		mv ${ref_tax}/thirty_closest_dists_30.txt ${ref_tax}/thirty_closest_dists.txt
	fi
	#exit
	if [[ ! -d ${ref_tax} ]]; then
		break
	else
		samples[sub_counter]=${sample}
	fi
	echo "${sample}"
	if [[ ${sub_counter} -lt ${max_subs} ]]; then
		echo  "Index is below max submissions, submitting"
		# echo -e "#!/bin/bash -l\n" > "${main_dir}/aniB_${sample}_${start_time}.sh"
		# echo -e "#$ -o aniB_${sample}.out" >> "${main_dir}/aniB_${sample}_${start_time}.sh"
		# echo -e "#$ -e aniB_${sample}.err" >> "${main_dir}/aniB_${sample}_${start_time}.sh"
		# echo -e "#$ -N aniB_${sample}"   >> "${main_dir}/aniB_${sample}_${start_time}.sh"
		# echo -e "#$ -cwd"  >> "${main_dir}/aniB_${sample}_${start_time}.sh"
		# echo -e "#$ -q short.q\n"  >> "${main_dir}/aniB_${sample}_${start_time}.sh"
		# echo -e "module load pyani/1.0" >> "${main_dir}/aniB_${sample}_${start_time}.sh"
		# echo -e "module load Python/3.5.2" >> "${main_dir}/aniB_${sample}_${start_time}.sh"
		# echo -e "average_nucleotide_identity.py -i \"${ref_tax}/localANIDB\" -o \"${ref_tax}/aniB\" -m \"ANIb\" -g \"--write_excel\"" >> "${main_dir}/aniB_${sample}_${start_time}.sh"
		# echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_aniB_complete.txt\"" >> "${main_dir}/aniB_${sample}_${start_time}.sh"
		# qsub "${main_dir}/aniB_${sample}_${start_time}.sh"

		echo -e "#!/bin/bash -l\n" > "${main_dir}/aniM_${sample}_${start_time}.sh"
		echo -e "#$ -o aniM_${sample}.out" >> "${main_dir}/aniM_${sample}_${start_time}.sh"
		echo -e "#$ -e aniM_${sample}.err" >> "${main_dir}/aniM_${sample}_${start_time}.sh"
		echo -e "#$ -N aniM_${sample}"   >> "${main_dir}/aniM_${sample}_${start_time}.sh"
		echo -e "#$ -cwd"  >> "${main_dir}/aniM_${sample}_${start_time}.sh"
		echo -e "#$ -q short.q\n"  >> "${main_dir}/aniM_${sample}_${start_time}.sh"
		echo -e "module load pyani/1.0" >> "${main_dir}/aniM_${sample}_${start_time}.sh"
		echo -e "module load Python/3.5.2" >> "${main_dir}/aniM_${sample}_${start_time}.sh"
		echo -e "average_nucleotide_identity.py -i \"${ref_tax}/localANIDB\" -o \"${ref_tax}/aniM\" \"--write_excel\"" >> "${main_dir}/aniM_${sample}_${start_time}.sh"
		echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_aniM_complete.txt\"" >> "${main_dir}/aniM_${sample}_${start_time}.sh"
		qsub "${main_dir}/aniM_${sample}_${start_time}.sh"

		echo -e "#!/bin/bash -l\n" > "${main_dir}/Fani_${sample}_${start_time}.sh"
		echo -e "#$ -o Fani_${sample}.out" >> "${main_dir}/Fani_${sample}_${start_time}.sh"
		echo -e "#$ -e Fani_${sample}.err" >> "${main_dir}/Fani_${sample}_${start_time}.sh"
		echo -e "#$ -N Fani_${sample}"   >> "${main_dir}/Fani_${sample}_${start_time}.sh"
		echo -e "#$ -cwd"  >> "${main_dir}/Fani_${sample}_${start_time}.sh"
		echo -e "#$ -q short.q\n"  >> "${main_dir}/Fani_${sample}_${start_time}.sh"
		echo -e "${shareScript}/fastANI --refList \"${ref_tax}/thirty_closest_dists.txt\" --query \"${reference}\" -t \"${procs}\" -o \"${ref_tax}/${sample}.fani\""  >> "${main_dir}/Fani_${sample}_${start_time}.sh"
		echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_Fani_complete.txt\"" >> "${main_dir}/Fani_${sample}_${start_time}.sh"
		qsub "${main_dir}/Fani_${sample}_${start_time}.sh"
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
			if [[ -f "${main_dir}/complete/${waiting_sample}_aniM_complete.txt" ]] && [[ -f "${main_dir}/complete/${waiting_sample}_Fani_complete.txt" ]]; then
				echo  "Index is below max submissions, submitting"
				# echo -e "#!/bin/bash -l\n" > "${main_dir}/ani_${sample}_${start_time}.sh"
				# echo -e "#$ -o aniB_${sample}.out" >> "${main_dir}/aniB_${sample}_${start_time}.sh"
				# echo -e "#$ -e aniB_${sample}.err" >> "${main_dir}/aniB_${sample}_${start_time}.sh"
				# echo -e "#$ -N aniB_${sample}"   >> "${main_dir}/aniB_${sample}_${start_time}.sh"
				# echo -e "#$ -cwd" >> "${main_dir}/aniB_${sample}_${start_time}.sh"
				# echo -e "#$ -q short.q\n"  >> "${main_dir}/aniB_${sample}_${start_time}.sh"
				# echo -e "module load pyani/1.0" "${main_dir}/aniB_${sample}_${start_time}.sh"
				# echo -e "average_nucleotide_identity.py -i \"${ref_tax}/localANIDB\" -o \"${ref_tax}/aniB\" -m \"ANIb\" -g \"--write_excel\"" >> "${main_dir}/aniB_${sample}_${start_time}.sh"
				# echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_aniB_complete.txt\"" >> "${main_dir}/aniB_${sample}_${start_time}.sh"
				# qsub "${main_dir}/aniB_${sample}_${start_time}.sh"

				echo -e "#!/bin/bash -l\n" > "${main_dir}/aniM_${sample}_${start_time}.sh"
				echo -e "#$ -o aniM_${sample}.out" >> "${main_dir}/aniM_${sample}_${start_time}.sh"
				echo -e "#$ -e aniM_${sample}.err" >> "${main_dir}/aniM_${sample}_${start_time}.sh"
				echo -e "#$ -N aniM_${sample}"   >> "${main_dir}/aniM_${sample}_${start_time}.sh"
				echo -e "#$ -cwd"  >> "${main_dir}/aniM_${sample}_${start_time}.sh"
				echo -e "#$ -q short.q\n"  >> "${main_dir}/aniM_${sample}_${start_time}.sh"
				echo -e "module load pyani/1.0" >> "${main_dir}/aniM_${sample}_${start_time}.sh"
				echo -e "module load Python/3.5.2" >> "${main_dir}/aniM_${sample}_${start_time}.sh"
				echo -e "average_nucleotide_identity.py -i \"${ref_tax}/localANIDB\" -o \"${ref_tax}/aniM\" \"--write_excel\"" >> "${main_dir}/aniM_${sample}_${start_time}.sh"
				echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_aniM_complete.txt\"" >> "${main_dir}/aniM_${sample}_${start_time}.sh"
				qsub "${main_dir}/aniM_${sample}_${start_time}.sh"

				echo -e "#!/bin/bash -l\n" > "${main_dir}/aniM_${sample}_${start_time}.sh"
				echo -e "#$ -o Fani_${sample}.out" >> "${main_dir}/Fani_${sample}_${start_time}.sh"
				echo -e "#$ -e Fani_${sample}.err" >> "${main_dir}/Fani_${sample}_${start_time}.sh"
				echo -e "#$ -N Fani_${sample}"   >> "${main_dir}/Fani_${sample}_${start_time}.sh"
				echo -e "#$ -cwd"  >> "${main_dir}/Fani_${sample}_${start_time}.sh"
				echo -e "#$ -q short.q\n"  >> "${main_dir}/Fani_${sample}_${start_time}.sh"
				echo -e "${shareScript}/fastANI --refList \"${ref_tax}/thirty_closest_dists.txt\" --query \"${reference}\" -t \"${procs}\" -o \"${ref_tax}/${sample}.fani\""  >> "${main_dir}/Fani_${sample}_${start_time}.sh"
				echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_Fani_complete.txt\"" >> "${main_dir}/Fani_${sample}_${start_time}.sh"
				qsub "${main_dir}/Fani_${sample}_${start_time}.sh"
			fi
		done
	fi
	sub_counter=$(( counter + 1 ))
 done
#echo "Count-${sub_counter}"
#exit
timer=0
waiting_sample=$(echo "${samples[${counter}]}")
while :
do
 	if [[ ${timer} -gt 1800 ]]; then
 		echo "Timer exceeded limit of 1800 seconds 30 minutes"
 		break
 	fi
 	if [[ -f "${main_dir}/complete/${waiting_sample}_aniM_complete.txt" ]] && [[ -f "${main_dir}/complete/${waiting_sample}_Fani_complete.txt" ]]; then
 		break
	else
		echo  "Waited ${timer}s on ${waiting_sample} to finish ANI's"
		sleep 5
		timer=$(( timer + 5 ))
	fi
done

#Extracts the query sample info line for percentage identity from the percent identity file
while IFS='' read -r line;
do
#	echo "!-${line}"
	if [[ ${line} == ${sample}* ]]; then
		sampleIMline=${line}
#		echo "found it-"$sampleline
		break
	fi
done < "${local_DBs}/aniDB/${working_dir}/${sample}/aniM/ANIm_percentage_identity.tab"

#Extracts the top line from the %id file to get all the sample names used in analysis (they are tab separated along the top row)
if [[ -s "${local_DBs}/aniDB/${working_dir}/${sample}/aniM/ANIm_percentage_identity.tab" ]]; then
	firstIMline=$(head -n 1 "${local_DBs}/aniDB/${working_dir}/${sample}/aniM/ANIm_percentage_identity.tab")
else
	echo "No "${local_DBs}/aniDB/${working_dir}/${sample}/aniM/ANIm_percentage_identity.tab" file, exiting"
	exit 1
fi

#Extracts the query sample info line for percentage identity from the percent identity file
while IFS='' read -r line;
do
#	echo "!-${line}"
	if [[ ${line} == ${sample}* ]]; then
		sampleCMline=${line}
#		echo "found it-"$sampleline
		break
	fi
done < "${local_DBs}/aniDB/${working_dir}/${sample}/aniM/ANIm_alignment_coverage.tab"

#Extracts the top line from the %id file to get all the sample names used in analysis (they are tab separated along the top row)
if [[ -s "${local_DBs}/aniDB/${working_dir}/${sample}/aniM/ANIm_alignment_coverage.tab" ]]; then
	firstCMline=$(head -n 1 "${local_DBs}/aniDB/${working_dir}/${sample}/aniM/ANIm_alignment_coverage.tab")
else
	echo "No "${local_DBs}/aniDB/${working_dir}/${sample}/aniM/ANIm_alignment_coverage.tab" file, exiting"
	exit 1
fi

#Arrays to read sample names and the %ids for the query sample against those other samples
IFS="	" read -r -a samples_aniM_identity <<< "${firstIMline}"
IFS="	" read -r -a percents_aniM_identity <<< "${sampleIMline}"
IFS="	" read -r -a samples_aniM_coverage <<< "${firstCMline}"
IFS="	" read -r -a percents_aniM_coverage <<< "${sampleCMline}"

echo "Making pyani_coverage_array"
declare -A pyani_coverage_array
counter=1
for isolate in "${samples_aniM_coverage[@]}"; do
	#echo "${isolate}"
	#temp_isolate=$(echo ${isolate} | cut -d'.' -f1)
	temp_isolate=${isolate//./_dot_}
	temp_percent=${percents_aniM_coverage[${counter}]}
	#echo "${temp_isolate}=${percents_aniM_coverage[${counter}]}=${temp_percent}"
	if [[ "${temp_percent}" == "1.0" ]]; then
		temp_percent=100
	else
		temp_percent_digits=$(echo ${temp_percent} | cut -d'.' -f2)
		temp_percent="${temp_percent_digits:0:2}.${temp_percent_digits:2}"
	fi
	pyani_coverage_array[${temp_isolate}]=${temp_percent}
	counter=$(( counter + 1 ))
done
for x in "${!pyani_coverage_array[@]}"; do printf "[%s]=%s\n" "$x" "${pyani_coverage_array[$x]}" ; done
echo "Making pyani_identity_array"
declare -A pyani_identity_array
counter=1
for isolate in "${samples_aniM_identity[@]}"; do
	#echo "${isolate}"
	#temp_isolate=$(echo ${isolate} | cut -d'.' -f1)
	temp_isolate=${isolate//./_dot_}
	temp_percent=${percents_aniM_identity[${counter}]}
	#echo "${temp_isolate}=${percents_aniM_identity[${counter}]}=${temp_percent}"
	if [[ "${temp_percent}" == "1.0" ]]; then
		temp_percent=100
	else
		temp_percent_digits=$(echo ${temp_percent} | cut -d'.' -f2)
		temp_percent="${temp_percent_digits:0:2}.${temp_percent_digits:2}"
	fi
	pyani_identity_array[${temp_isolate}]=${temp_percent}
	counter=$(( counter + 1 ))
done
for x in "${!pyani_idenity_array[@]}"; do printf "[%s]=%s\n" "$x" "${pyani_identity_array[$x]}" ; done
echo "Making fastANI_identity_array"
#Extracts the query sample info line for ANI
declare -A fastANI_identity_array
while IFS='' read -r line;
do
	temp_isolate=$(echo ${line} | cut -d' ' -f2 | rev | cut -d'/' -f1 | cut -d'.' -f2- | rev)
	temp_percent=$(echo ${line} | cut -d' ' -f3)

	#echo "Tax:${temp_isolate}"
	#echo "Temp_tax:${temp_isolate_tax}"
	#echo "Sample:${sample}"
	#echo "%:${temp_percent}"
	#temp_isolate=$(echo ${tax} | cut -d'.' -f1)
	#echo "${temp_isolate}-${temp_percent}"
	temp_isolate=${temp_isolate//./_dot_}
	fastANI_identity_array[${temp_isolate}]=${temp_percent}
done < "${local_DBs}/aniDB/${working_dir}/${sample}/${sample}.fani"
for x in "${!fastANI_identity_array[@]}"; do printf "[%s]=%s\n" "$x" "${fastANI_identity_array[$x]}" ; done
#echo "4"
declare -A mash_dist_array
declare -A mash_kmers_array
if [[ ! -f "${local_DBs}/aniDB/${working_dir}/${sample}/${sample}.dists" ]]; then
	cp "${local_DBs}/aniDB/all_test/dists/${sample}.dists" "${local_DBs}/aniDB/${working_dir}/${sample}"
	head -n 31 ${local_DBs}/aniDB/${working_dir}/${sample}/${sample}.dists > ${local_DBs}/aniDB/${working_dir}/${sample}/${sample}_top30.dists
fi
while IFS='' read -r line;
do
	isolate=$(echo ${line} | cut -d' ' -f1 | rev | cut -d'/' -f1 | cut -d'.' -f2- | rev)
	temp_isolate=${isolate//./_dot_}
	mash_dist=$(echo ${line} | cut -d' ' -f3)
	mash_kmer=$(echo ${line} | cut -d' ' -f5)
	#echo "Tax:${temp_isolate}"
	#echo "dist:${mash_dist}"
	#echo "kmer:${mash_kmer}"
	#temp_isolate=$(echo ${tax} | cut -d'.' -f1)
	mash_dist_array[${temp_isolate}]=${mash_dist}
	mash_kmers_array[${temp_isolate}]=${mash_kmer}
done < "${local_DBs}/aniDB/${working_dir}/${sample}/${sample}_top30.dists"
for x in "${!mash_dist_array[@]}"; do printf "[%s]=%s\n" "$x" "${mash_dist_array[$x]}" ; done
for x in "${!mash_kmers_array[@]}"; do printf "[%s]=%s\n" "$x" "${mash_kmers_array[$x]}" ; done
#echo "5"
if [[ -f ${local_DBs}/aniDB/${working_dir}/${sample}/${sample}_ani_summary.tsv ]]; then
	rm -r ${local_DBs}/aniDB/${working_dir}/${sample}/${sample}_ani_summary.tsv
fi

for isolate in "${samples_aniM_identity[@]}"; do
	#temp_isolate=$(echo ${isolate} | rev | cut -d'.' -f2 | rev)
	#echo "A"
	temp_isolate=${isolate//./_dot_}
	temp_isolate_tax=$(echo ${temp_isolate} | cut -d'_' -f1,2)


	pyani_percent_ID=${pyani_identity_array[${temp_isolate}]}
	#echo "B"
	pyani_coverage=${pyani_coverage_array[${temp_isolate}]}
	#echo "C"
	fastANI_percent_ID=${fastANI_identity_array[${temp_isolate}]}
	if [[ "${temp_isolate_tax}" == "${sample}" ]]; then
		fastANI_percent_ID=100
	fi
	mash_dist=${mash_dist_array[${temp_isolate}]}
	mash_kmer=${mash_kmers_array[${temp_isolate}]}
	#echo "D"
	#echo "${isolate}:${temp_isolate}:${pyani_percent_ID}:${pyani_coverage}:${fastANI_percent_ID}:${mash_dist}:${mash_kmer}"

	if [[ -z ${fastANI_percent_ID} ]]; then
		fastANI_percent_ID="<<80"
	fi
	echo -e "${isolate}	${pyani_percent_ID}	${pyani_coverage}	${fastANI_percent_ID}	${mash_dist}	${mash_kmer}" >> ${local_DBs}/aniDB/${working_dir}/${sample}/${sample}_ani_summary.tsv
done
sort -k2 -nr -o ${local_DBs}/aniDB/${working_dir}/${sample}/${sample}_ani_summary.tsv ${local_DBs}/aniDB/${working_dir}/${sample}/${sample}_ani_summary.tsv
tail -n +2 "${local_DBs}/aniDB/${working_dir}/${sample}/${sample}_ani_summary.tsv" > "${local_DBs}/aniDB/${working_dir}/${sample}/${sample}_ani_summary.tsv.tmp" && mv "${local_DBs}/aniDB/${working_dir}/${sample}/${sample}_ani_summary.tsv.tmp" "${local_DBs}/aniDB/${working_dir}/${sample}/${sample}_ani_summary.tsv"
echo -e "ANI summary for ${sample}\nreference	pyani_%_ID	pyani_coverage	fastANI_%_ID	mash_dist mash_kmers\n $(cat ${local_DBs}/aniDB/${working_dir}/${sample}/${sample}_ani_summary.tsv)" > ${local_DBs}/aniDB/${working_dir}/${sample}/${sample}_ani_summary.tsv
