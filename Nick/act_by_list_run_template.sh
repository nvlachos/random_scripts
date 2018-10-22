#!/bin/sh -l

#$ -o act_by_list.out
#$ -e act_by_list.err
#$ -N act_by_list
#$ -cwd
#$ -q all.q

#Import the config file with shortcuts and settings
. /scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/scripts/config.sh
#Import the module file that loads all necessary mods
. "${mod_changers}/pipeline_mods"
. "${mod_changers}/list_modules.sh"

#
# Usage ./act_by_list.sh list_name(currently has to be placed in /scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR folder)
#
# script changes depending on what needs to be run through the list
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to act_by_list.sh, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./act_by_list.sh path_to_list_file(single sample ID per line, e.g. B8VHY/1700128 (it must include project id also))"
	echo "Output location varies depending on which tasks are performed but will be found somewhere under ${processed}"
	exit 0
fi

# Loop through and act on each sample name in the passed/provided list
while IFS= read -r var; do
		sample_name=$(echo "${var}" | awk -F/ '{ print $2}' | tr -d '[:space:]')
		project=$(echo "${var}" | awk -F/ '{ print $1}' | tr -d '[:space:]')
		OUTDATADIR="${processed}/${project}/${sample_name}"
		current_time=$(date "+%m-%d-%Y @ %Hh_%Mm_%Ss")
		echo "Starting ${project}/${sample_name} through loop at ${current_time}"
	##### Download yourself (it can be done, but a pain to find and get single samples for both accessing and storing without messing up durrent configuration)
	##### Unzip files (Also more complicated due to stripped S numbers, will try again later???)
	### Check QC counts on raw reads
	#if [ ! -d "${OUTDATADIR}/preQCcounts" ]; then
	#	mkdir -p "${OUTDATADIR}/preQCcounts"
	#fi
	#	python $shareScript/Fastq_Quality_Printer.py "${OUTDATADIR}/FASTQs/${sample_name}_R1_001.fastq" "${OUTDATADIR}/FASTQs/${sample_name}_R2_001.fastq" > "${OUTDATADIR}/preQCcounts/${sample_name}_counts.txt"
	### Run BBDuk
	#if [ ! -d "${OUTDATADIR}/removedAdapters" ]; then 
	#	mkdir -p "${OUTDATADIR}/removedAdapters"
	#fi
	#	bbduk.sh -$bbduk_mem threads=${procs} in="${OUTDATADIR}/FASTQs/${sample_name}_R1_001.fastq" in2="${OUTDATADIR}/FASTQs/${sample_name}_R2_001.fastq" out="${OUTDATADIR}/removedAdapters/${sample_name}-noPhiX-R1.fsq" out2="${OUTDATADIR}/removedAdapters/${sample_name}-noPhiX-R2.fsq" ref="$shareScript/phiX.fasta k=$bbduk_k hdist=$bbduk_hdist"
	### Run Trimmomatic
	#if [ ! -d "${OUTDATADIR}/trimmed" ]; then
	#	mkdir -p "${OUTDATADIR}/trimmed"
	#fi
	#	trimmomatic $trim_endtype -$trim_phred -threads $procs ${OUTDATADIR}/removedAdapters/${sample_name}-noPhiX-R1.fsq ${OUTDATADIR}/removedAdapters/${sample_name}-noPhiX-R2.fsq ${OUTDATADIR}/trimmed/${sample_name}_R1_001.paired.fq ${OUTDATADIR}/trimmed/${sample_name}_R1_001.unpaired.fq ${OUTDATADIR}/trimmed/${sample_name}_R2_001.paired.fq ${OUTDATADIR}/trimmed/${sample_name}_R2_001.unpaired.fq ILLUMINACLIP:$trim_adapter_location:$trim_seed_mismatch:$trim_palindrome_clip_threshold:$trim_simple_clip_threshold:$trim_min_adapt_length:$trim_complete_palindrome SLIDINGWINDOW:$trim_window_size:$trim_window_qual LEADING:$trim_leading TRAILING:$trim_trailing MINLEN:$trim_min_length
	### Concatenate read files after trimmomatic
	#	cat ${OUTDATADIR}/trimmed/${sample_name}_R1_001.paired.fq ${OUTDATADIR}/trimmed/${sample_name}_R2_001.paired.fq > ${OUTDATADIR}/trimmed/${sample_name}.paired.fq
	#	cat ${OUTDATADIR}/trimmed/${sample_name}_R1_001.unpaired.fq ${OUTDATADIR}/trimmed/${sample_name}_R2_001.unpaired.fq > ${OUTDATADIR}/trimmed/${sample_name}.single.fq
	### Check QC counts on trimmed reads
	#	python $shareScript/Fastq_Quality_Printer.py "${OUTDATADIR}/trimmed/${sample_name}.paired.fq" "${OUTDATADIR}/trimmed/${sample_name}.single.fq" > "${OUTDATADIR}/preQCcounts/${sample_name}_trimmed_counts.txt"
	### Run Kraken on reads
	#   ${shareScript}/run_kraken.sh ${sample_name} pre paired ${project}
	#   ${shareScript}/best_hit_from_kraken.sh ${sample_name} pre paired ${project}
	### Run Gottcha on reads
	#   ${shareScript}/run_gottcha.sh ${sample_name} ${project}
	### Get best hits from gottcha
	#   ${shareScript}/best_hit_from_gottcha1.sh ${sample_name} ${project}
	### Assemble with SPAdes
	#   ${shareScript}/run_SPAdes.sh ${sample_name} normal ${project}
	#	echo "${sample_name} plasmid ${project}"
	#if [ ! -f ${OUTDATADIR}/plasmidAssembly/scaffolds.fasta ]; then
	#   ${shareScript}/run_SPAdes.sh ${sample_name} plasmid ${project}
	#fi
	### Remove short contigs and rename assembly file
	#   perl ${shareScript}/removeShortContigs.pl ${OUTDATADIR}/Assembly/scaffolds.fasta
	#   mv ${OUTDATADIR}/Assembly/scaffolds.fasta.TRIMMED.fasta ${OUTDATADIR}/Assembly/${sample_name}_scaffolds_trimmed.fasta
	### Remove short contigs and rename plasmid assembly file
	#   perl ${shareScript}/removeShortContigs.pl ${OUTDATADIR}/plasmidAssembly/scaffolds.fasta
	#   mv ${OUTDATADIR}/plasmidAssembly/scaffolds.fasta.TRIMMED.fasta ${OUTDATADIR}/plasmidAssembly/${sample_name}_plasmid_scaffolds_trimmed.fasta
	### Run Kraken on assembly
	#   ${shareScript}/run_kraken.sh ${sample_name} post assembled ${project}
	#   ${shareScript}/best_hit_from_kraken.sh ${sample_name} post assembled ${project}
	#	"${shareScript}/best_hit_from_kraken.sh" "${sample_name}" "post" "assembled_BP_data" "${project}"
	### FIND TAXONOMY in kraken output
						while IFS= read -r line;
						do
							first=${line::1}
							if [ "${first}" = "S" ]
							then
								species=$(echo "${line}" | awk -F ' ' '{print $4}')
							elif [ "${first}" = "G" ]
							then
								genus=$(echo "${line}" | awk -F ' ' '{print $4}')
	#						elif [ "${first}" = "F" ]
	#						then
	#							family=$(echo "${line}" | awk -F ' ' '{print $4}')
	#						elif [ "${first}" = "O" ]
	#						then
	#							order=$(echo "${line}" | awk -F ' ' '{print $4}')
	#						elif [ "${first}" = "C" ]
	#						then
	#							class=$(echo "${line}" | awk -F ' ' '{print $4}')
	#						then
	#							phylum=$(echo "${line}" | awk -F ' ' '{print $4}')
	#						elif [ "${first}" = "K" ]
	#						then
	#							kingdom=$(echo "${line}" | awk -F ' ' '{print $4}')
	#						elif [ "${first}" = "D" ]
	#						then
	#							domain=$(echo "${line}" | awk -F ' ' '{print $4}')
							fi
						done < "${OUTDATADIR}/kraken/postAssembly/${sample_name}_kraken_summary_assembled_BP_data.txt"
	### Check Assembly QC
	#   ${shareScript}/run_Assembly_Quality_Check.sh ${sample_name} ${project}
	### prokka
	#   ${shareScript}/run_prokka.sh ${sample_name} ${project}
	### busco
	#buscoDB="bacteria_odb9"
	# Iterate through taxon levels (species to domain) and test if a match occurs to entry in database. If so, compare against it
	#busco_found=0
	#for tax in $species $genus $family $order $class $phylum $kingdom $domain
	#do
	#	if [ -d "${share}/DBs/BUSCO/${tax,}_odb9" ]
	#	then
	#		buscoDB="${tax}_odb9"
	#		busco_found=1
	#		break
	#	fi
	#done
	# Report an unknown sample to the maintenance file to look into
	#if [[ "${busco_found}" -eq 0 ]]; then
	#	global_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")
	#	echo "BUSCO: ${domain} ${kingdom} ${phylum} ${class} ${order} ${family} ${genus} ${species} - ${project}/${sample_name} at ${global_time}"
	#else
	#	# Show which database entry will be used for comparison
	#	echo "buscoDB:${buscoDB}"
	#	# Run busco
	#	"${shareScript}/do_busco.sh" "${sample_name}" "${buscoDB}" "${project}"
	#fi
	### ANI
	#if [ ! -f ${OUTDATADIR}/ANI/best_ANI_hits_ordered(${sample_name}_vs_${genus}).txt ]; then
	#   "${shareScript}/run_ANI.sh" "${sample_name}" "${genus}" "${species}" "${project}"
	#fi
	### Run csstar in gapped mode at high similarity (98%) on full assembly
	#	`python -V`
	#	rm -r "${processed}/${project}/${sample_name}/c-sstar"
	#	"${shareScript}/run_c-sstar_on_single.sh" "${sample_name}" "${csstar_gapping}" "${csstar_identity}" "${project}"
	### Run csstar in gapped mode at other similarity (40%) on plasmid assembly
	#	rm -r "${processed}/${project}/${sample_name}/c-sstar_plasmid"
	#	"${shareScript}/run_c-sstar_on_single.sh" "${sample_name}" "${csstar_gapping}" "${csstar_plasmid_identity}" "${project}" "--plasmid"
	### Run SRST2
	#	rm -r "${processed}/${project}/${sample_name}/srst2"
	#	"${shareScript}/run_srst2_on_singleDB.sh" "${sample_name}" "${project}"
	#mv "${processed}/${project}/${sample_name}/srst2/resgannot__${sample_name}.ResGANNOT_20180409_srst2.pileup" "${processed}/${project}/${sample_name}/srst2/ResGANNOT__${sample_name}.ResGANNOT_20180409_srst2.pileup"
	#mv "${processed}/${project}/${sample_name}/srst2/resgannot__${sample_name}.ResGANNOT_20180409_srst2.sorted.bam" "${processed}/${project}/${sample_name}/srst2/ResGANNOT__${sample_name}.ResGANNOT_20180409_srst2.sorted.bam"
	#mv "${processed}/${project}/${sample_name}/srst2/resgannot__fullgenes__ResGANNOT_20180409_srst2__results.txt" "${processed}/${project}/${sample_name}/srst2/ResGANNOT__fullgenes__ResGANNOT_20180409_srst2__results.txt"
	#mv "${processed}/${project}/${sample_name}/srst2/resgannot__genes__ResGANNOT_20180409_srst2__results.txt" "${processed}/${project}/${sample_name}/srst2/ResGANNOT__genes__ResGANNOT_20180409_srst2__results.txt"

	### Get MLST profile
	#if [[ ! -f ${processed}/${project}/${sample_name}/MLST/${sample_name}.mlst_oxf ]]; then
	#	mv "${processed}/${project}/${sample_name}/MLST/${sample_name}.mlst" "${processed}/${project}/${sample_name}/MLST/${sample_name}.mlst2"
	#	${shareScript}/run_MLST.sh "${sample_name}" "${project}" "-f" "abaumannii"
	#	mv "${processed}/${project}/${sample_name}/MLST/${sample_name}.mlst" "${processed}/${project}/${sample_name}/MLST/${sample_name}.mlst_oxf"
	#	mv "${processed}/${project}/${sample_name}/MLST/${sample_name}.mlst2" "${processed}/${project}/${sample_name}/MLST/${sample_name}.mlst"
	#fi
	### Get ID from 16s
	# "${shareScript}/16s_blast.sh" "${sample_name}" "${project}"
	# check for plasmids
	#if [[ -d "${processed}/${project}/${sample_name}/plasmid" ]]; then
	#	echo "Removing ${processed}/${project}/${sample_name}/plasmid"
	#	rm -r "${processed}/${project}/${sample_name}/plasmid"
	#fi
	#"${shareScript}/run_plasmidFinder.sh" "${sample_name}" "${project}" "plasmid"
	# check for plasmids on plasmid assembly
	#if [[ -d "${processed}/${project}/${sample_name}/plasmid_on_plasmidAssembly" ]]; then
	#	echo "${processed}/${project}/${sample_name}/plasmid_on_plasmidAssembly"
	#	rm -r "${processed}/${project}/${sample_name}/plasmid_on_plasmidAssembly"
	#fi
	#"${shareScript}/run_plasmidFinder.sh" "${sample_name}" "${project}" "plasmid_on_plasmidAssembly"
	### validate sample
	#   ${shareScript}/validate_piperun.sh "${sample_name}" ${project} > ${OUTDATADIR}/${sample_name}_pipeline_stats.txt
	### Add sample validation to run validation
	#	${shareScript}/validate_piperun.sh "${sample_name}" ${project} >> ${share}/total_check.txt
		current_time=$(date "+%m-%d-%Y @ %Hh_%Mm_%Ss")
done < "${share}/${1}"
echo "All samples completed"
global_end_time=$(date "+%m-%d-%Y @ %Hh_%Mm_%Ss")
#Script exited gracefully (unless something else inside failed)
printf "%s %s" "Act_by_list.sh has completed whatever it was doing at" "${global_end_time}" | mail -s "act_by_list complete" nvx4@cdc.gov
exit 0
