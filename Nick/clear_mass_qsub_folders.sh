#!/bin/sh -l

#$ -o clear_mqs.out
#$ -e clear_mqs.err
#$ -N clear_mqs
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
. ./config.sh
#Import the module file that loads all necessary mods
. "${mod_changers}/pipeline_mods"

#
# Usage ./clear_mass_qsub_fodlers [1,2,3] (1-qsub folders, 2-qsub outs/errs, 3-Both)
#
# script changes naming structure of SPAdes output to include isolate name for every contig and removes coverage info
#
if [[ ${1} -eq 1 ]] || [[ ${1} -eq 3 ]]; then
	for folder in ${share}/mass_subs/*
	do
		if [[ -d ${folder} ]]; then
			for folder1 in ${folder}
			do
				echo "Deleting every sh script in ${folder1}"
				rm -r ${folder1}/*.sh
				echo "Deleting every sh script in ${folder1}/complete"
				rm -r ${folder1}/complete/*.txt
			done
		fi
	done
fi

if [[ ${1} -eq 2 ]] || [[ ${1} -eq 3 ]]; then
	rm -r ${shareScript}/blast16sID_*.out
	rm -r ${shareScript}/blast16sID_*.err
	rm -r ${shareScript}/blast16s_*.out
	rm -r ${shareScript}/blast16s_*.err
	rm -r ${shareScript}/ani_*.out
	rm -r ${shareScript}/ani_*.err
	rm -r ${shareScript}/ANI_*.out
	rm -r ${shareScript}/ANI_*.err
	rm -r ${shareScript}/BTQC_*.out
	rm -r ${shareScript}/BTQC_*.err
	rm -r ${shareScript}/BUSCO_*.out
	rm -r ${shareScript}/BUSCO_*.err
	rm -r ${shareScript}/getFASTQR1_*.out
	rm -r ${shareScript}/getFASTQR1_*.err
	rm -r ${shareScript}/getFASTQR2_*.out
	rm -r ${shareScript}/getFASTQR2_*.err
	rm -r ${shareScript}/csstn_*.out
	rm -r ${shareScript}/csstn_*.err
	rm -r ${shareScript}/csstp_*.out
	rm -r ${shareScript}/csstp_*.err
	rm -r ${shareScript}/kraka_*.out
	rm -r ${shareScript}/kraka_*.err
	rm -r ${shareScript}/krakr_*.out
	rm -r ${shareScript}/krakr_*.err
	rm -r ${shareScript}/gott_*.out
	rm -r ${shareScript}/gott_*.err
	rm -r ${shareScript}/mlst_*.out
	rm -r ${shareScript}/mlst_*.err
	rm -r ${shareScript}/pFinf_*.out
	rm -r ${shareScript}/pFinf_*.err
	rm -r ${shareScript}/pFinp_*.out
	rm -r ${shareScript}/pFinp_*.err
	rm -r ${shareScript}/SPAdn_*.out
	rm -r ${shareScript}/SPAdn_*.err
	rm -r ${shareScript}/SPAdp_*.out
	rm -r ${shareScript}/SPAdp_*.err
	rm -r ${shareScript}/plasFlow_*.out
	rm -r ${shareScript}/plasFlow_*.err
	rm -r ${shareScript}/pFinf_*.out
	rm -r ${shareScript}/pFinf_*.err
	rm -r ${shareScript}/PROKK_*.out
	rm -r ${shareScript}/PROKK_*.err
	rm -r ${shareScript}/PROKK_*.e*
	rm -r ${shareScript}/srst2AR_*.out
	rm -r ${shareScript}/srst2AR_*.err
	rm -r ${shareScript}/srst2MLST_*.out
	rm -r ${shareScript}/srst2MLST_*.err
	rm -r ${shareScript}/srst22MLST_*.out
	rm -r ${shareScript}/srst22MLST_*.err
	rm -r ${shareScript}/QUAST_*.out
	rm -r ${shareScript}/QUAST_*.err
	rm -r ${shareScript}/QC_*.out
	rm -r ${shareScript}/QC_*.err
	rm -r ${shareScript}/MLST_*.out
	rm -r ${shareScript}/MLST_*.err
	rm -r ${shareScript}/taxID_*.out
	rm -r ${shareScript}/taxID_*.err
	rm -r ${shareScript}/validate_*.out
	rm -r ${shareScript}/validate_*.err
	rm -r ${shareScript}/sum_*.out
	rm -r ${shareScript}/sum_*.err
	rm -r ${shareScript}/pFn_*.out
	rm -r ${shareScript}/pFn_*.err
	rm -r ${shareScript}/pFp_*.out
	rm -r ${shareScript}/pFp_*.err
fi
