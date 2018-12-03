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
        for genes_2 in range(len(Gene_List_out)):
            if str(Gene_List[genes].seq) == str(Gene_List_out[genes_2].seq):
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
            elif str(Gene_List[genes].seq.reverse_complement()) == str(Gene_List_out[genes_2].seq):
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
                    if base not in ("A","C","G","T"):
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

def Dual_Fasta_Combiner(res_File, arg_file, new_file):
    """Reads in Two fasta files and writes them out to a single new file"""
    print(res_File,"and",arg_file)
    f1 = open(res_File, 'r')
    f2 = open(new_file, 'w')
    for line in f1:
        f2.write(line.replace('>', '>[RES]').replace('/', '-'))
    f1.close()
    f2.close()
    f1 = open(arg_file, 'r')
    f2 = open(new_file, 'a')
    for line in f1:
        if (line[0] == ">") and (line[1] != "("):
            resist=line[1:4]
            oldline=line[4:]
            newline=">[ARG]("+resist+")"+oldline
            f2.write(newline)
        else:
            f2.write(line.replace('>', '>[ARG]').replace('/','-'))

Dual_Fasta_Combiner(sys.argv[1], sys.argv[2], sys.argv[3])
Duplicate_Gene_Remover(sys.argv[3], sys.argv[3], sys.argv[4])
