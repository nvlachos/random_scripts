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

working_dir="Single_ANI_Test"

# Sets the output folder to the sample_name folder in processed samples
OUTDATADIR="${local_DBs}/aniDB"
#if [[ ! -d ${OUTDATADIR}/${working_dir}/dists ]]; then
#	mkdir -p ${OUTDATADIR}/${working_dir}/dists
#fi

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
#	mash dist "${local_DBs}/aniDB/all/all_sketch.msh" "${ref}" > "${OUTDATADIR}/${working_dir}/dists/${filename}_unsorted.dists"
#	sort -k3 -n -o "${OUTDATADIR}/${working_dir}/dists/${filename}.dists" "${OUTDATADIR}/${working_dir}/dists/${filename}_unsorted.dists"
#	rm -r "${OUTDATADIR}/${working_dir}/dists/${filename}_unsorted.dists"
#done

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
#
# exit

sub_counter=0
max_subs=50
samples=()
main_dir="/scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/mass_subs/ani_TEST"
mkdir ${main_dir}


for ref_tax in ${local_DBs}/aniDB/${working_dir}/*; do
	echo "${ref_tax}"
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
# 	if [[ ${sub_counter} -lt ${max_subs} ]]; then
# 		echo  "Index is below max submissions, submitting"
# 		# echo -e "#!/bin/bash -l\n" > "${main_dir}/aniB_${sample}_${start_time}.sh"
# 		# echo -e "#$ -o aniB_${sample}.out" >> "${main_dir}/aniB_${sample}_${start_time}.sh"
# 		# echo -e "#$ -e aniB_${sample}.err" >> "${main_dir}/aniB_${sample}_${start_time}.sh"
# 		# echo -e "#$ -N aniB_${sample}"   >> "${main_dir}/aniB_${sample}_${start_time}.sh"
# 		# echo -e "#$ -cwd"  >> "${main_dir}/aniB_${sample}_${start_time}.sh"
# 		# echo -e "#$ -q short.q\n"  >> "${main_dir}/aniB_${sample}_${start_time}.sh"
# 		# echo -e "module load pyani/1.0" >> "${main_dir}/aniB_${sample}_${start_time}.sh"
# 		# echo -e "module load Python/3.5.2" >> "${main_dir}/aniB_${sample}_${start_time}.sh"
# 		# echo -e "average_nucleotide_identity.py -i \"${ref_tax}/localANIDB\" -o \"${ref_tax}/aniB\" -m \"ANIb\" -g \"--write_excel\"" >> "${main_dir}/aniB_${sample}_${start_time}.sh"
# 		# echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_aniB_complete.txt\"" >> "${main_dir}/aniB_${sample}_${start_time}.sh"
# 		# qsub "${main_dir}/aniB_${sample}_${start_time}.sh"
#
# 		echo -e "#!/bin/bash -l\n" > "${main_dir}/aniM_${sample}_${start_time}.sh"
# 		echo -e "#$ -o aniM_${sample}.out" >> "${main_dir}/aniM_${sample}_${start_time}.sh"
# 		echo -e "#$ -e aniM_${sample}.err" >> "${main_dir}/aniM_${sample}_${start_time}.sh"
# 		echo -e "#$ -N aniM_${sample}"   >> "${main_dir}/aniM_${sample}_${start_time}.sh"
# 		echo -e "#$ -cwd"  >> "${main_dir}/aniM_${sample}_${start_time}.sh"
# 		echo -e "#$ -q short.q\n"  >> "${main_dir}/aniM_${sample}_${start_time}.sh"
# 		echo -e "module load pyani/1.0" >> "${main_dir}/aniM_${sample}_${start_time}.sh"
# 		echo -e "module load Python/3.5.2" >> "${main_dir}/aniM_${sample}_${start_time}.sh"
# 		echo -e "average_nucleotide_identity.py -i \"${ref_tax}/localANIDB\" -o \"${ref_tax}/aniM\" \"--write_excel\"" >> "${main_dir}/aniM_${sample}_${start_time}.sh"
# 		echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_aniM_complete.txt\"" >> "${main_dir}/aniM_${sample}_${start_time}.sh"
# 		qsub "${main_dir}/aniM_${sample}_${start_time}.sh"
#
# 		echo -e "#!/bin/bash -l\n" > "${main_dir}/aniM_${sample}_${start_time}.sh"
# 		echo -e "#$ -o Fani_${sample}.out" >> "${main_dir}/Fani_${sample}_${start_time}.sh"
# 		echo -e "#$ -e Fani_${sample}.err" >> "${main_dir}/Fani_${sample}_${start_time}.sh"
# 		echo -e "#$ -N Fani_${sample}"   >> "${main_dir}/Fani_${sample}_${start_time}.sh"
# 		echo -e "#$ -cwd"  >> "${main_dir}/Fani_${sample}_${start_time}.sh"
# 		echo -e "#$ -q short.q\n"  >> "${main_dir}/Fani_${sample}_${start_time}.sh"
# 		echo -e "${shareScript}/fastANI --refList \"${ref_tax}/thirty_closest_dists.txt\" --query \"${reference}\" -t \"${procs}\" -o \"${ref_tax}/${sample}.fani\""  >> "${main_dir}/Fani_${sample}_${start_time}.sh"
# 		echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_Fani_complete.txt\"" >> "${main_dir}/Fani_${sample}_${start_time}.sh"
# 		qsub "${main_dir}/Fani_${sample}_${start_time}.sh"
# 	else
# 		waiting_for_index=$(( counter - max_subs ))
# 		waiting_sample=$(echo "${samples[${waiting_for_index}]}")
# 		timer=0
# 		echo "Index is above max submissions, waiting for index ${waiting_for_index}:${waiting_sample} to complete"
# 		while :
# 		do
# 			if [[ ${timer} -gt 1800 ]]; then
# 				echo "Timer exceeded limit of 1800 seconds 30 minutes"
# 				break
# 			fi
# 			if [[ -f "${main_dir}/complete/${waiting_sample}_aniM_complete.txt" ]] && [[ -f "${main_dir}/complete/${waiting_sample}_Fani_complete.txt" ]]; then
# 				echo  "Index is below max submissions, submitting"
# 				# echo -e "#!/bin/bash -l\n" > "${main_dir}/ani_${sample}_${start_time}.sh"
# 				# echo -e "#$ -o aniB_${sample}.out" >> "${main_dir}/aniB_${sample}_${start_time}.sh"
# 				# echo -e "#$ -e aniB_${sample}.err" >> "${main_dir}/aniB_${sample}_${start_time}.sh"
# 				# echo -e "#$ -N aniB_${sample}"   >> "${main_dir}/aniB_${sample}_${start_time}.sh"
# 				# echo -e "#$ -cwd" >> "${main_dir}/aniB_${sample}_${start_time}.sh"
# 				# echo -e "#$ -q short.q\n"  >> "${main_dir}/aniB_${sample}_${start_time}.sh"
# 				# echo -e "module load pyani/1.0" "${main_dir}/aniB_${sample}_${start_time}.sh"
# 				# echo -e "average_nucleotide_identity.py -i \"${ref_tax}/localANIDB\" -o \"${ref_tax}/aniB\" -m \"ANIb\" -g \"--write_excel\"" >> "${main_dir}/aniB_${sample}_${start_time}.sh"
# 				# echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_aniB_complete.txt\"" >> "${main_dir}/aniB_${sample}_${start_time}.sh"
# 				# qsub "${main_dir}/aniB_${sample}_${start_time}.sh"
#
# 				echo -e "#!/bin/bash -l\n" > "${main_dir}/aniM_${sample}_${start_time}.sh"
# 				echo -e "#$ -o aniM_${sample}.out" >> "${main_dir}/aniM_${sample}_${start_time}.sh"
# 				echo -e "#$ -e aniM_${sample}.err" >> "${main_dir}/aniM_${sample}_${start_time}.sh"
# 				echo -e "#$ -N aniM_${sample}"   >> "${main_dir}/aniM_${sample}_${start_time}.sh"
# 				echo -e "#$ -cwd"  >> "${main_dir}/aniM_${sample}_${start_time}.sh"
# 				echo -e "#$ -q short.q\n"  >> "${main_dir}/aniM_${sample}_${start_time}.sh"
# 				echo -e "module load pyani/1.0" >> "${main_dir}/aniM_${sample}_${start_time}.sh"
# 				echo -e "average_nucleotide_identity.py -i \"${ref_tax}/localANIDB\" -o \"${ref_tax}/aniM\" \"--write_excel\"" >> "${main_dir}/aniM_${sample}_${start_time}.sh"
# 				echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_aniM_complete.txt\"" >> "${main_dir}/aniM_${sample}_${start_time}.sh"
# 				qsub "${main_dir}/aniM_${sample}_${start_time}.sh"
#
# 				echo -e "#!/bin/bash -l\n" > "${main_dir}/aniM_${sample}_${start_time}.sh"
# 				echo -e "#$ -o Fani_${sample}.out" >> "${main_dir}/Fani_${sample}_${start_time}.sh"
# 				echo -e "#$ -e Fani_${sample}.err" >> "${main_dir}/Fani_${sample}_${start_time}.sh"
# 				echo -e "#$ -N Fani_${sample}"   >> "${main_dir}/Fani_${sample}_${start_time}.sh"
# 				echo -e "#$ -cwd"  >> "${main_dir}/Fani_${sample}_${start_time}.sh"
# 				echo -e "#$ -q short.q\n"  >> "${main_dir}/Fani_${sample}_${start_time}.sh"
# 				echo -e "${shareScript}/fastANI --refList \"${ref_tax}/thirty_closest_dists.txt\" --query \"${reference}\" -t \"${procs}\" -o \"${ref_tax}/${sample}.fani\""  >> "${main_dir}/Fani_${sample}_${start_time}.sh"
# 				echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_Fani_complete.txt\"" >> "${main_dir}/Fani_${sample}_${start_time}.sh"
# 				qsub "${main_dir}/Fani_${sample}_${start_time}.sh"
# 			fi
# 		done
# 	fi
# 	sub_counter=$(( counter + 1 ))
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
declare -A coverage_array
counter=0
for isolate in "${samples_aniM_coverage[@]}"; do
	#echo "${isolate}"
	temp_isolate=$(echo ${isolate} | cut -d'.' -f1)
	echo "${temp_isolate}-${percents_aniM_coverage[${counter}]}"
	pyani_coverage_array[${temp_isolate}]=${percents_aniM_coverage[${counter}]}
	counter=$(( counter + 1 ))
done
echo "Making pyani_identity_array"
declare -A identity_array
counter=0
for isolate in "${samples_aniM_identity[@]}"; do
	#echo "${isolate}"
	temp_isolate=$(echo ${isolate} | cut -d'.' -f1)
	echo "${temp_isolate}-${percents_aniM_identity[${counter}]}"
	pyani_identity_array[${temp_isolate}]=${percents_aniM_identity[${counter}]}
	counter=$(( counter + 1 ))
done
echo "Making fastANI_identity_array"
#Extracts the query sample info line for ANI
declare -A fastANI_identity_array
while IFS='' read -r line;
do
	current_tax=$(echo ${line} | cut -d'	' -f2 | cut -d' ' -f1 | rev | cut -d'/' -f1 | cut -d'.' -f2- | rev)
	current_percent=$(echo ${line} | cut -d'	' -f3 | cut -d' ' -f1)
	echo "Tax:${current_tax}"
	echo "%:${current_percent}"
	temp_isolate=$(echo ${tax} | cut -d'.' -f1)
	echo "${temp_isolate}-${current_percent}"
	fastANI_identity_array[${temp_isolate}]=${current_percent}
done < "${local_DBs}/aniDB/${working_dir}/${sample}/${sample}.fani"
echo "4"
if [[ -f ${local_DBs}/aniDB/${working_dir}/${sample}/${sample}_ani_summary.tsv ]]; then
	rm -r ${local_DBs}/aniDB/${working_dir}/${sample}/${sample}_ani_summary.tsv
fi
echo "ANI summary for ${sample}" >> ${local_DBs}/aniDB/${working_dir}/${sample}/${sample}_ani_summary.tsv
echo -e "reference	pyani_%_ID	pyani_coverage	fastANI_%_ID" >> ${local_DBs}/aniDB/${working_dir}/${sample}/${sample}_ani_summary.tsv
for isolate in "${samples_aniM_identity[@]}"; do
	temp_isolate=$(echo ${isolate} | cut -d'.' -f1)
	pyani_percent_ID=${pyani_identity_array[${temp_isolate}]}
	pyani_coverage=${pyani_coverage_array[${temp_isolate}]}
	fastANI_percent_ID=${fastANI_identity_array[${temp_isolate}]}
	if [[ -z ${fastANI_percent_ID} ]]; then
		fastANI_percent_ID="<<80"
	fi
	echo -e "${isolate}	${pyani_percent_ID}	${pyani_coverage}	${fastANI_percent_ID}" >> ${local_DBs}/aniDB/${working_dir}/${sample}/${sample}_ani_summary.tsv
done
