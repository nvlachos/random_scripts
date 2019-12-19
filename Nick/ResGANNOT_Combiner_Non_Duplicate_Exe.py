import sys
import Bio
import glob
import fileinput
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet.IUPAC import unambiguous_dna
from Bio.Seq import MutableSeq
from Bio.Alphabet import IUPAC
from Bio import SeqIO
from decimal import *
getcontext().prec = 3

def Gene_Counter(input_fasta):
    Gene_List = list(SeqIO.parse(input_fasta, 'fasta'))
    Gene_Count = len(Gene_List)
    return Gene_Count

def check_alphabet(sequence,  code="ATGCatgc"):
    for base in sequence:
        if base not in code:
            return False
    return True

def Duplicate_Gene_Remover(input_fasta, output_file, output_copy_file):
    Gene_List = list(SeqIO.parse(input_fasta, 'fasta'))
    Gene_List_out = [Gene_List[0]]
    Gene_copies_out = []
    Gene_Adder = 0
    output_fasta = open(output_file, 'w')
    copies_fasta = open(output_copy_file, 'w')
    bad_seqs=[]
    for genes in range(1, len(Gene_List)):
        if genes > 0:
            for genes_2 in range(len(Gene_List_out)):
                if str(Gene_List[genes].seq).upper() == str(Gene_List_out[genes_2].seq).upper():
                    Gene_Adder = 1
                    #print("Wrong1")
                    #print(Gene_List[genes].description)
                    #print(Gene_List_out[genes_2].description)
                    #Gene_copies_out.append(Gene_List[genes].description)
                    #Gene_copies_out.append(Gene_List[genes].seq)
                    #Gene_copies_out.append(Gene_List_out[genes_2].description)
                    #Gene_copies_out.append(Gene_List_out[genes_2].seq)
                    Gene_copies_out.append(Gene_List[genes])
                    Gene_copies_out.append(Gene_List_out[genes_2])
                    break
                elif str(Gene_List[genes].seq.reverse_complement()).upper() == str(Gene_List_out[genes_2].seq).upper():
                    Gene_Adder = 1
                    #print("Wrong2")
                    #print(Gene_List[genes].description)
                    #print(Gene_List_out[genes_2].description)
                    #Gene_copies_out.append(Gene_List[genes].description)
                    #Gene_copies_out.append(Gene_List[genes].seq)
                    #Gene_copies_out.append(Gene_List_out[genes_2].description)
                    #Gene_copies_out.append(Gene_List_out[genes_2].seq)
                    Gene_copies_out.append(Gene_List[genes])
                    Gene_copies_out.append(Gene_List_out[genes_2])
                    break
                else:
                    continue
        if Gene_Adder == 0:
            # check that sequence is bases only
			#temp_seq=str(Gene_List[genes].seq).replace("\\", "")
			#Gene_List[genes] = Seq(temp_seq)
            if (check_alphabet(str(Gene_List[genes].seq))):
##              print(Gene_List[genes].format('fasta'))
                Gene_List_out.append(Gene_List[genes])
            else:
                errors=""
                for base in str(Gene_List[genes].seq):
                    if base not in ("A","C","G","T","a","c","g","t"):
                        if (errors == ""):
                            errors = base
                        else:
                            errors = errors + "," + base
                Gene_List[genes].description = Gene_List[genes].description + " (Non-standard Base(s) - " + errors + ")"
                bad_seqs.append(Gene_List[genes])
        else:
            Gene_Adder = 0
    for output in Gene_copies_out:
        #print(output)
        SeqIO.write(output, copies_fasta, 'fasta')
    copies_fasta.close()
    for output_genes in Gene_List_out:
        SeqIO.write(output_genes, output_fasta, 'fasta')
    output_fasta.close()
    bad_output_fasta = open(output_file+".bad", "w")
    for output_genes in bad_seqs:
        SeqIO.write(output_genes, bad_output_fasta, 'fasta')
    bad_output_fasta.close()

def Tri_Fasta_Combiner(file_1, file_2, file_3, out_file,t1,t2,t3):
    """Reads in Two fasta files and writes them out to a single new file, if using an ARG it must always be second"""
    title_one=">["+t1+"]"
    title_two=">["+t2+"]"
    title_three=">["+t3+"]"
    print(file_1+"\n"+file_2+"\n"+file_3)
    f1 = open(file_1, 'r')
    fo = open(out_file, 'w')
    for line in f1:
        fo.write(line.replace('>', title_one).replace('/', '-'))
    f1.close()
    f2 = open(file_2, 'r')
    for line in f2:
        if (line[0] == ">") and (line[1] != "(") and t2 == "ARG":
            resist=line[1:4]
            oldline=line[4:]
            newline=title_two+"("+resist+")"+oldline
            fo.write(newline)
        else:
            fo.write(line.replace('>', title_two).replace('/','-'))
    f2.close()
    f3 = open(file_3, 'r')
    for line in f3:
        if line[0] == ">" and t3 == "NCBI":
            ncbi_list=line.split("|")
            ngene=ncbi_list[5]
            ncount=ncbi_list[3]
            naccession=ncbi_list[2]
            nrange=line.split(" ")[1].split(":")[1]
            newline=title_three+ngene+"_"+ncount+"_"+naccession+"_"+nrange
            fo.write(newline.replace('/', '-'))
        else
            fo.write(line.replace('>', title_three).replace('/','-'))
    f3.close()
    fo.close()



Tri_Fasta_Combiner(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[6], sys.argv[7], sys.argv[8])
Duplicate_Gene_Remover(sys.argv[4], sys.argv[4], sys.argv[5])
