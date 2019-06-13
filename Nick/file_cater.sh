for filename in ${1}/*.blast.best; do
	sample=$(basename $filename .blast.best)
	#echo $sample
	while IFS=read -r var  || [ -n "$var" ]; do
		node=$(echo "${var}" | cut -d'	' -f2)
		ref_start=$(echo "${var}" | cut -d'	' -f7)
		ref_stop=$(echo "${var}" | cut -d'	' -f8)
		query_start=$(echo "${var}" | cut -d'	' -f9)
		query_stop=$(echo "${var}" | cut -d'	' -f10)
		ref_size=$(echo "${var}" | cut -d'	' -f13)
		rev_comp="false"
		echo "rsta=${ref_start};rsto=${ref_stop};qsta=${query_start};qsto=${query_stop};rsiz=${ref_size}"
		if	[[ $query_start -lt $query_stop ]]; then
			lookup_start=$(( query_start - (ref_start - 1) ))
			lookup_stop=$(( query_stop + (ref_size - ref_stop) ))
		else
			lookup_start=$(( query_stop - (ref_size - ref_stop) ))
			lookup_stop=$(( query_start + (ref_start - 1) ))
			rev_comp="true"
		fi
		echo "Node:$node	Start:${lookup_start}   Stop:${lookup_stop} REV:${rev_comp}"
		sequence=$(samtools faidx ${1}/Assemblies/${sample}_scaffolds_trimmed.fasta ${node}:${lookup_start}-${lookup_stop} | tail -n1)
		#sequence=$(echo "${sequence}" | cut -d'\n' -f2)
		if [[ "${rev_comp}" = "true" ]]; then
			echo "Reversing: ${sequence}"
			sequence=$(echo "$sequence" | tr "[ATGCatgc]" "[TACGtacg]" | rev)
		fi
		echo "Now: ${sequence}"
		echo "${sample}	${var}	${sequence}" >> "${1}/concat_${2}.txt"
	done < ${filename}
done
