import sys
import Bio
import glob
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet.IUPAC import unambiguous_dna
from Bio.Seq import MutableSeq
from Bio.Alphabet import IUPAC
from Bio import SeqIO

def Contig_Writer(input_fasta, contig, output_fasta):
    """Writes out a fasta file of just the contig in question"""
    Genome = SeqIO.to_dict(SeqIO.parse(input_fasta, 'fasta'))
    Output = open(output_fasta, 'w')
    SeqIO.write(Genome[contig], Output, 'fasta')
    Output.close()

Contig_Writer(sys.argv[1], sys.argv[2], sys.argv[3])
