import sys
import glob
import math
import itertools as it

# main function that looks if all MLST types are defined for an outptu mlst file
def do_MLST_check(input_MLST_file, MLST_filetype):
	types=""
	schemes=[]
	MLST_file=open(input_MLST_file,'r')
	MLST_line=MLST_file.readline().strip()
	MLST_items=MLST_line.split("	")
	MLST_file.close()
	#print("\n".join(MLST_items))
	if MLST_filetype == "standard":
		sample=MLST_items[0]
		db_name=MLST_items[1]
		mlst_temp_type=MLST_items[2].replace("/", ",")
		if "," not in mlst_temp_type:
			mlstype=[MLST_items[2]]
			for i in range(0, len(mlstype)):
				if mlstype[i] != '-':
					mlstype[i] = int(mlstype[i])
			mlstype.sort()
		else:
			mlstype=MLST_items[2].split(",")
		print("Current MLST type:", mlstype, "\n")
		allele_list=[]
		allele_names=[]
		allele_count=len(MLST_items)
		for allele in range(3, allele_count):
			#print(MLST_items[allele])
			allele_Identifier=MLST_items[allele].split("(")[0]
			alleles=MLST_items[allele].split("(")[1].split(")")[0].split(",")
			allele_names.append(allele_Identifier)
			allele_list.append(alleles)
			list_size=len(allele_list)
		#allele_list=[['1'], ['3'], ['189','3'], ['2'], ['2'], ['96','107'], ['3']]
		print("Allele_names:", allele_names)
		print("Alleles_found:", allele_list, "\n")
		#for allele_index in range(0,len(allele_list)):
		#	allele_list[allele_index]=allele_list[allele_index].sort()
		if list_size == 7:
			schemes = it.product(allele_list[0], allele_list[1], allele_list[2], allele_list[3], allele_list[4], allele_list[5], allele_list[6])
		elif list_size == 8:
			schemes = it.product(allele_list[0], allele_list[1], allele_list[2], allele_list[3], allele_list[4], allele_list[5], allele_list[6], allele_list[7])
		else:
			print("Unknown size "+str(list_size)+" of allele_list")
		schemes=(list(schemes))
		print("All possible schemes:")
		print(*schemes, sep = "\n")
		print()

		checking=False
		for profile_index in range(0, len(schemes)):
			#print(profile_index, schemes[profile_index])
			temp_scheme=[]
			for temp_allele in schemes[profile_index]:
				temp_scheme.append(temp_allele)
			schemes[profile_index]=temp_scheme
			#print(profile_index, schemes[profile_index])
		if len(schemes) == 0:
			print("No schemes found???")
		elif len(schemes) == 1:
			if mlstype[0] != "-":
		 		print("This sample is singular and defined\n")
			else:
				print("This sample is singular and UNdefined\n")
				new_types=get_type(schemes, allele_names, db_name)
				checking=True
		elif len(schemes) > 1:
			if "-" not in mlstype:
				if len(schemes) == len(mlstype):
					print("This sample is a multiple and defined\n")
				elif len(schemes) > len(mlstype):
					print("Not enough types to match schemes")
					new_types=get_type(schemes, allele_names, db_name)
					checking=True
				elif len(schemes) < len(mlstype):
					print("Not enough schemes to match types")
					new_types=get_type(schemes, allele_names, db_name)
					checking=True
			else:
				print("This sample is a multiple and something is UNdefined")
				new_types=get_type(schemes, allele_names, db_name)
				checking=True
		print("Old types:", mlstype, "\n")
		if checking:
			print("New types:", new_types, "\n")
			if mlstype != new_types:
				for i in range(0, len(new_types)):
					print(new_types[i])
					new_types[i] = str(new_types[i])
				#new_types.sort()
				new_types=','.join(new_types)
				print("Updating MLST types in", input_MLST_file, "from", mlstype, "to", new_types)
				MLST_items[2]=new_types
				new_info='	'.join(MLST_items)
				MLST_file=open(input_MLST_file,'w')
				MLST_file.write(new_info)
				MLST_file.close()
		else:
			print(input_MLST_file, "is as good as it gets with type", mlstype)
		if '-' in MLST_items[2]:
			filepath=input_MLST_file[::-1].split("/")[2:4]
			print(filepath)
			for i in range(0, len(filepath)):
				print(filepath[i])
				filepath[i]=filepath[i][::-1]
			filepath=filepath[::-1]
			filepath="/".join(filepath)
			print("Must try srst2 on input:", filepath)
	elif MLST_filetype == "srst2":
		print("Not implemented yet")
	else:
		print("Unknown MLST filetype, can not continue")
		exit()
	counter=0


def get_type(list_of_profiles, list_of_allele_names, DB_file):
	types=["Not_initialized"]
	full_db_path="/scicomp/groups/OID/NCEZID/DHQP/CEMB/databases/pubmlsts/"+DB_file+"/"+DB_file+".txt"
	with open(full_db_path,'r') as f:
		profile_size=0
		types = ["-"] * len(list_of_profiles)
		print("Size:", len(types), " &  contents:", types)
		for line in f:
			db_line=line.strip()
			db_items=db_line.split("	")
			if db_items[0] == "ST":
				for item in db_items:
					if item != "clonal_complex" and item != "species":
						profile_size+=1
					else:
						break
				print(db_items[1:profile_size])
				print(list_of_allele_names)
				if db_items[1:profile_size] == list_of_allele_names:
					print("Allele names match, yay!")
				else:
					print("We'll have to fix this if it ever comes up")
					print("db: "+db_items)
					print("list:"+allele_names)
			else:
				for index in range(0,len(types)):
					current_profile=db_items[1:profile_size]
					type(current_profile)
					type(list_of_profiles)
					#print(current_profile)
					#print(list_of_profiles[index])
					if current_profile == list_of_profiles[index]:
						print("Match-"+str(db_items[0]), current_profile)
						types[index] = int(db_items[0])
						break
	types.sort()
	return types





































	# #print("Start")
	# while csstar_line != '':
	# 	#print(counter, csstar_line)
	# 	#print("Start csstar loop")
	# 	csstar_line_sections=csstar_line.split("	")
	# 	ar_list=csstar_line_sections[7].split(",")
	# 	ar_dict={}
	# 	for ar_gene in ar_list:
	# 		gene_name=ar_gene.split("[")[0]
	# 		if gene_name == "No AR genes discovered":
	# 			gene_stats="[0/0:#-]C"
	# 		else:
	# 			gene_stats="["+ar_gene.split("[")[1]+"C"
	# 		ar_dict[gene_name]=gene_stats
	# 		if gene_name not in all_ARs_in_file:
	# 			all_ARs_in_file.append(gene_name)
	# 			#print("Adding", gene_name)
	# 	srst2_file=open(input_srst2_AR,'r')
	# 	srst2_line=srst2_file.readline().strip()
	# 	counter=0
	# 	while srst2_line != '':
	# 		#print("Start srst2 loop")
	# 		spot_count=0
	# 		#for k, v in ar_dict.items():
	# 		#	print(spot_count)
	# 		#	print(k, v)
	# 		#	spot_count+=1
	# 		#print("Checking", srst2_line)
	# 		srst2_line_sections=srst2_line.split("	")
	# 		if csstar_line_sections[0] == srst2_line_sections[0] and csstar_line_sections[1] == srst2_line_sections[1]:
	# 			#print("Found",  csstar_line_sections[1], "in srst2 summary file")
	# 			if srst2_line_sections[2] == "No AR genes discovered":
	# 				gene_name="No AR genes discovered"
	# 				gene_stats="[0/0]S"
	# 				#print("Looking up", gene_name, "in csstar dic")
	# 				if ar_dict.get(gene_name):
	# 					ar_dict[gene_name]=""+ar_dict.get(gene_name)+":"+gene_stats
	# 				else:
	# 					print("New AR-less isolate found in srst2")
	# 					ar_dict[gene_name]=gene_stats
	# 				if gene_name not in all_ARs_in_file:
	# 					all_ARs_in_file.append(gene_name)
	# 			else:
	# 				srst2_ar_list=srst2_line_sections[2].split(",")
	# 				for srst2_ar_gene in srst2_ar_list:
	# 					gene_name=srst2_ar_gene.split("[")[0]
	# 					gene_stats="["+srst2_ar_gene.split("[")[1]+"S"
	# 					#print("Looking up", gene_name, "in csstar dic")
	# 					if ar_dict.get(gene_name):
	# 						if ar_dict.get(gene_name) != "No Other AR genes":
	# 							#print("Found", gene_name, "in both outputs")
	# 							#print("New value: "+ar_dict.get(gene_name)+":"+gene_stats)
	# 							ar_dict[gene_name]=""+ar_dict.get(gene_name)+":"+gene_stats
	# 						#else:
	# 							#print("No AR found in csstar for", gene_name)
	# 						#	:
	# 					else:
	# 						print("New gene", gene_name,"found in srst2")
	# 						ar_dict[gene_name]=gene_stats
	# 					if gene_name not in all_ARs_in_file:
	# 						all_ARs_in_file.append(gene_name)
	# 			break
	# 		#else:
	# 		#	#print(csstar_line_sections[1], "C does not equal S", srst2_line_sections[1])
	# 		#	:
	# 		srst2_line=srst2_file.readline().strip()
	# 		counter+=1
	# 	srst2_file.close()
	# 	#for k, v in ar_dict.items():
	# 	#	print(k, v)
	# 	#print("1:",csstar_line_sections[0])
	# 	#print("0:", csstar_line_sections[0], "1:", csstar_line_sections[1],"2:" , csstar_line_sections[2], "3:", csstar_line_sections[3])
	# 	samples.append([csstar_line_sections[0], csstar_line_sections[1], csstar_line_sections[2], csstar_line_sections[3],  csstar_line_sections[4], csstar_line_sections[5], csstar_line_sections[6], ar_dict])
	# 	#print("Total AR genes in sample set:", len(all_ARs_in_file)-1)
	# 	csstar_line = csstar_file.readline().strip()
	# csstar_file.close
	# all_ARs_in_file.sort()
	# if len(all_ARs_in_file) == 0:
	# 	print("\n")
	# 	print("Total AR genes in sample set: 0")
	# else:
	# 	print("Total AR genes in sample set:",len(all_ARs_in_file))
	# 	#print(*all_ARs_in_file, sep = "\n")
	# 	for gene in all_ARs_in_file:
	# 		if gene != "No other AR genes":
	# 			print (gene)
	# print()
	#
	#
	#
# 	#Parse plasmid summary file
# 	all_plasmids_in_file=[]
# 	plas_file=open(input_plas, 'r')
# 	line = plas_file.readline().strip()
# 	sample_p_plasmids_dict={}
# 	sample_f_plasmids_dict={}
# 	current_id=""
# 	counter=0
# 	while line != '':
# 		#print(counter, line)
# 		plasmid_line_sections=line.split("	")
# 		#print("Current id:", current_id, ":", plasmid_line_sections[0])
# 		if current_id == "":
# 			current_id = plasmid_line_sections[0]+"/"+plasmid_line_sections[1]
# 		if plasmid_line_sections[0]+"/"+plasmid_line_sections[1] != current_id:
# 			#print("New name!")
# 			for sample_index in range(0,len(samples)):
# 				#print("Looking for", current_id, ", found", samples[sample_index][0].strip())
# 				if current_id == samples[sample_index][0]+"/"+samples[sample_index][1].strip():
# 					#print(samples[sample_index], "adding", sample_f_plasmids, "and", sample_p_plasmids)
# 					samples[sample_index].append(sample_f_plasmids_dict)
# 					samples[sample_index].append(sample_p_plasmids_dict)
# 					break
# 				#else:
# 					#print("Sample", current_id, "does not exist")
# 			current_id=plasmid_line_sections[0]+"/"+plasmid_line_sections[1].strip()
# 			sample_f_plasmids_dict={}
# 			sample_p_plasmids_dict={}
# 		source_assembly=plasmid_line_sections[2]
# 		#print("Test:"+plasmid_line_sections[4]+":")
# 		#if plasmid_line_sections[4].find("_contigs-") >= 0:
# 		#	line = plas_file.readline().strip()
# 		#	continue
#
# 		if plasmid_line_sections[3] == "No_Plasmids_Found":
# 			plas_perc_id="-"
# 			plas_perc_length="-"
# 		else:
# 			plas_perc_id=math.floor(float(plasmid_line_sections[4]))
# 			#print("testing:", plasmid_line_sections[5].split("/")[0], plasmid_line_sections[5].split("/")[1])
# 			plas_perc_length=(100*int(plasmid_line_sections[5].split("/")[1])//int(plasmid_line_sections[5].split("/")[0]))
# 		#plas_match_info="["+plas_perc_id+"/"+plas_percpercent_length+"]"
# 		if source_assembly == "full_assembly":
# 			#print("Adding:", plasmid_line_sections[3], "to sample_f_plasmids")
# 			sample_f_plasmids_dict[plasmid_line_sections[3]]="["+str(plas_perc_id)+"/"+str(plas_perc_length)+"]"
# 		elif source_assembly == "plasmid_assembly":
# 			#print("Adding:", plasmid_line_sections[3], "to sample_p_plasmids")
# 			sample_p_plasmids_dict[plasmid_line_sections[3]]="["+str(plas_perc_id)+"/"+str(plas_perc_length)+"]"
# 		if len(plasmid_line_sections) > 1:
# 			if plasmid_line_sections[3] not in all_plasmids_in_file:
# 				all_plasmids_in_file.append(plasmid_line_sections[3])
# 				#print("Adding to project list", plasmid_line_sections[3])
# 		#else:
# 			#print("Line length:", len(csstar_plasmid_line_sections))
# 		#print()
# 		line = plas_file.readline().strip()
# 		counter=counter+1
# 	for sample_index in range(0,len(samples)):
# 		#print("Looking for", current_id, ", found", samples[sample_index][0].strip())
# 		if current_id == samples[sample_index][0]+"/"+samples[sample_index][1].strip():
# 			#print(samples[sample_index], "adding", sample_f_plasmids, "and", sample_p_plasmids)
# 			samples[sample_index].append(sample_f_plasmids_dict)
# 			samples[sample_index].append(sample_p_plasmids_dict)
# 			break
# 		#else:
# 			#print("Sample", current_id, "does not exist")
# 	all_plasmids_in_file.sort()
# 	if len(all_plasmids_in_file) == 0:
# 		#print("\n")
# 		print("Total plasmid replicons in sample set: 0")
# 	else:
# 		print("Total plasmid replicons in sample set:", len(all_plasmids_in_file)-1)
# 		print(*all_plasmids_in_file, sep= "\n")
# 	print()
# 	plas_file.close
# 	all_ar_and_plasmids=all_ARs_in_file+["|"]+all_plasmids_in_file
# 	#all_AR_to_write=all_ARs_in_file
# 	#all_AR_to_write.insert(0,",")
# 	#all_AR_to_write.insert(0,",")
# 	#all_AR_to_write=','.join(map(str, all_AR_to_write))
# 	header="id, Project__autocolour, Species__autocolour, Species_determinant__autocolour, Species_Support__autocolour , MLST__autocolour, ALT_MLST__autocolour,"
# 	for thing in all_ar_and_plasmids:
# 		header = header + " " + thing + "__autocolour,"
# 	header = header[:-1]
# 	summary_out=open(output_file, 'w')
# 	summary_out.write(header+'\n')
# 	#all_AR_to_write=all_AR_to_write[2:]
# 	#print("List:", all_ar_and_plasmids)
# 	#for sample in samples:
# 	#	print ("2:",sample[0])
# 	#return
# 	for sample in samples:
# 		sample_details=[sample[1], sample[0], sample[2], sample[3], sample[4], sample[5], sample[6]]
# 		#print("pre:",sample)
# 		for gene in all_ar_and_plasmids:
# 			status=" "
# 			if gene == "|":
# 				sample_details.append(gene)
# 				continue
# 			if sample[7].get(gene):
# 				status=sample[7].get(gene)
# 			elif sample[8].get(gene):
# 				if sample[9].get(gene):
# 					status="F:"+sample[8].get(gene)+";P:"+sample[9].get(gene)
# 				else:
# 					status="F:"+sample[8].get(gene)
# 			elif sample[8].get(gene):
# 				status="P:"+sample[9].get(gene)
# 			sample_details.append(status)
# 		#print("Post Sample check", sample_details)
# 		sample_details=','.join(map(str, sample_details))
# 		summary_out.write(sample_details+"\n")
# 	summary_out.close
# #print (sys.argv[1:])



print("Parsing MLST file ...\n")
do_MLST_check(sys.argv[1], sys.argv[2]) #, "/scicomp/groups/OID/NCEZID/DHQP/CEMB/databases/mlst/abaumannii_Pasteur.txt") #sys.argv[3])
