import sys
import Bio
import glob
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet.IUPAC import unambiguous_dna
from Bio.Alphabet import IUPAC
from Bio import SeqIO

def Gene_Overlap_Finder(fasta_1, fasta_2):
    """Reads in multi gene fasta files and returns list of genes absent from second fasta file"""
    Genome_1 = list(SeqIO.parse(fasta_1, 'fasta'))
    Genome_2 = list(SeqIO.parse(fasta_2, 'fasta'))
    Missing_Genes = []
    for genes_1 in Genome_1:
        Gene_Adder = 1
        for genes_2 in Genome_2:
            if str(genes_1.seq) == str(genes_2.seq):
                Gene_Adder = 0
                break
        if Gene_Adder == 1:
            Missing_Genes.append(genes_1.id)
    return Missing_Genes

def Genome_Difference(fasta_1, fasta_2, new_fasta):
    """Makes a fasta of novel genes present only in fasta_2"""
    Gene_List = (Gene_Overlap_Finder(fasta_2, fasta_1))
    Genome = SeqIO.to_dict(SeqIO.parse(fasta_2, 'fasta'))
    output_handle = open(new_fasta, 'w')
    for genes in Gene_List:
        SeqIO.write(Genome[genes], output_handle, 'fasta')
    output_handle.close()

## Usage: python Genome_Difference_Exe.py Fasta_1 Fasta_2 Output_Fasta
## Output_Fasta contains all genes from Fasta_2 absent from Fasta_1

Genome_Difference(sys.argv[1], sys.argv[2], sys.argv[3])
