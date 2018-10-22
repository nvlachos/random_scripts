from __future__ import division
import sys
import Bio
import glob
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet.IUPAC import unambiguous_dna
from Bio.Alphabet import IUPAC 
from Bio.Alphabet import generic_dna
from Bio.Seq import MutableSeq
from Bio.Alphabet import IUPAC
from Bio import Seq
from Bio import SeqIO
from Bio import pairwise2
from operator import itemgetter
from itertools import product

def String_Repeat_Remover(input_string):
        """returns a string with any repeated characters removed"""
        Out_String = ''
        for characters in input_string:
                if (characters in Out_String) == False:
                        Out_String = Out_String + characters
        return Out_String

##def In_Line(value, line_string):
##    if len(str(value)) > len(line_string):
##        return False
##    else:
##        for characters in range(0, (len(line_string)-len(value))):
##            if value == line_string[characters:(characters + len(value))]:
##                return True
##    return False    


def Find_Match(value, line_string):
    """Returns the start and stop positions of the value string in line_string"""
    if (value in line_string) == False:
        return "No match!"
    else:
        Out_List = []
        for characters in range(0, (len(line_string)-len(value))):
            if value == line_string[characters:(characters + len(value))]:
                Out_List.append(characters)
                Out_List.append(characters + len(value))
                return Out_List

def Find_Rev_Match(value, line_string):
        """Returns the start and stop positions of the reverse compliment of the value string in line_string"""
        Value_Seq = SeqRecord(Seq(value, generic_dna))
        Rev_Seq = str(Value_Seq.seq.reverse_complement())
        return Find_Match(Rev_Seq, line_string)
 
def Amplicon_Length(gene, primer_f, primer_r):
    """Takes in a gene and 2 primers and returns the length of the amplicon product"""
    primer_rev_seq = SeqRecord(Seq(primer_r, generic_dna))
    primer_rev = str(primer_rev_seq.seq.reverse_complement())
    if (primer_f in gene) == False:
        return "Forward primer doesn't match!"
    elif (primer_rev in gene) == False:
        return "Reverse primer doesn't match!"
    else:
        Length_List = []
        Length_List.append(Find_Match(primer_f, gene)[0])
        Length_List.append(Find_Match(primer_rev, gene)[1])
        Length = Length_List[1] - Length_List[0]
        return Length

def Amplicon_RF_Length(gene, primer_f, primer_r):
    """Takes in a gene and 2 primers in the same RF and returns the length of the amplicon product"""
    if (primer_f in gene) == False:
        return "Forward primer doesn't match!"
    elif (primer_r in gene) == False:
        return "Reverse primer doesn't match!"
    else:
        Length_List = []
        Length_List.append(Find_Match(primer_f, gene)[0])
        Length_List.append(Find_Match(primer_r, gene)[1])
        Length = Length_List[1] - Length_List[0]
        return Length

def Amplicon_RF_Finder(gene, primer_f, primer_r):
    """Takes in a gene and 2 primers and returns the amplicon product"""
    if (primer_f in gene) == False:
        return "Forward primer doesn't match!"
    elif (primer_r in gene) == False:
        return "Reverse primer doesn't match!"
    else:
        Length_List = []
        Length_List.append(Find_Match(primer_f, gene)[0])
        Length_List.append(Find_Match(primer_r, gene)[1])
        Amplicon = gene[Length_List[0]:Length_List[1]]
        return Amplicon

def Amplicon_Finder(gene, primer_f, primer_r):
    """Takes in a gene and 2 primers and returns the amplicon product"""
    primer_rev_seq = SeqRecord(Seq(primer_r, generic_dna))
    primer_rev = str(primer_rev_seq.seq.reverse_complement())
    if (primer_f in gene) == False:
        return "Forward primer doesn't match!"
    elif (primer_rev in gene) == False:
        return "Reverse primer doesn't match!"
    else:
        Length_List = []
        Length_List.append(Find_Match(primer_f, gene)[0])
        Length_List.append(Find_Match(primer_rev, gene)[1])
        Amplicon = gene[Length_List[0]:Length_List[1]]
        return Amplicon

def Amplicon_Finder_Assembly(assembly, primer_f, primer_r):
        """Takes in an assembly and returns the contig/gene name"""
        Genome = list(SeqIO.parse(assembly, 'fasta'))
        primer_rev_seq = SeqRecord(Seq(primer_r, generic_dna))
        primer_rev = str(primer_rev_seq.seq.reverse_complement())
        for gene in Genome:
                if (primer_f in str(gene.seq)) == False:
                        continue
                elif (primer_rev in str(gene.seq)) == False:
                        continue
                else:
                        return gene.id
        else:
                return "No matches!"

def GC_Content(input_string):
	CG = 0
	for bases in input_string:
		if bases == 'C' or bases == 'G':
			CG = CG + 1
	GC_total = float(CG) / len(input_string)
	return GC_total

def T_Melting(input_string):
	AT = 0
	CG = 0
	for bases in input_string:
		if bases == 'C' or bases == 'G':
			CG = CG + 1
		else:
			AT = AT + 1
	Tmelt = 2*AT + 4*CG
	return Tmelt

def GC_Repeats(input_string):
        """returns True if G/C repeats more than 3 bases, otherwise false"""
        for characters in range(0, len(input_string) - 3):
                if GC_Content(input_string[characters:(characters + 4)]) == 1:
                        return True
        else:
                return False

def Long_Repeats(input_string):
        """returns True if there is a string of repeats >4 of any base, otherwise false"""
        for characters in range(0, len(input_string) - 4):
                if len(String_Repeat_Remover(input_string[characters:(characters + 5)])) == 1:
                        return True
        else:
                return False
        
def Conserved_N_List(input_fasta, number):
    Out_List = []
    Gene_List = list(SeqIO.parse(input_fasta, 'fasta'))
    Shortest_Length = len(str(Gene_List[0].seq))
    Shortest_Gene = Gene_List[0]
    for genes in Gene_List:
        if len(str(genes.seq)) < Shortest_Length:
            Shortest_Length = len(str(genes.seq))
            Shortest_Gene = genes
    Shortest_String = str(Shortest_Gene.seq)
    if len(Shortest_String) > number:
        List_Number = []
        for numbers in range(0, len(Shortest_String) - number):
            List_Number.append(Shortest_String[numbers:(numbers + number)])
        Counter = 1
        for strings in List_Number:
            Gene_Match = 1
            for genes in Gene_List:
                Gene_String = str(genes.seq)
                if (strings in Gene_String) == False:
                    Gene_Match = 0
                    break
                else:
                    continue
            if Gene_Match == 1:
                Out_List.append(strings)
            Counter = Counter + 1
        return Out_List

def Conserved_Primers(input_fasta, gene_name, number):
        """Takes in input multi-fasta of gene variants and makes a list of possible primers"""
        Mer_List =  Conserved_N_List(input_fasta, number)
        Genes = SeqIO.to_dict(SeqIO.parse(input_fasta, 'fasta'))
        Gene = str(Genes[gene_name].seq)
        List1 = []
        for mers in Mer_List:
                for mers_2 in Mer_List:
                        Length = Amplicon_RF_Length(Gene, mers, mers_2)
                        if Length <= 150 and Length >= 75:
                                Pair = []
                                Pair.append(mers)
                                Pair.append(mers_2)
                                Pair.append(Length)
                                List1.append(Pair)
        List2 = []
        for Pairs in List1:
                Amplicon = Amplicon_RF_Finder(Gene, Pairs[0], Pairs[1])
                if (GC_Content(Amplicon) >= 0.5 and GC_Content(Amplicon) <= 0.6) and Long_Repeats == False:
                        List2.append(Pairs)
        return List2
        List3 = []
        for Pairs in List2:
                if (T_Melting(Pairs[0]) <=65 and T_Melting(Pairs[0]) >=50) and (T_Melting(Pairs[1]) <=65 and T_Melting(Pairs[1]) >=50):
                        Pairs.append(T_Melting(Pairs[0]))
                        Pairs.append(T_Melting(Pairs[1]))
                        Pairs.append(abs(T_Melting(Pairs[0]) - T_Melting(Pairs[1])))
                        List3.append(Pairs)
        List4 = []
        for Pairs in List3:
                if ((Pairs[0][0] == 'G' or Pairs[0][0] == 'C') and (Pairs[0][-1] == 'G' or Pairs[0][-1] == 'C')) and ((Pairs[1][0] == 'G' or Pairs[1][0] == 'C') and (Pairs[1][-1] == 'G' or Pairs[1][-1] == 'C')):
                        List4.append(Pairs)
        List5 = []
        for Pairs in List4:
                if ((GC_Content(Pairs[0]) >= 0.5 and GC_Content(Pairs[0]) <= 0.6) and (GC_Content(Pairs[1]) >= 0.5 and GC_Content(Pairs[1]) <= 0.6)):
                        List5.append(Pairs)
        List6 = []
        for Pairs in List5:
                if GC_Repeats(Pairs[0]) == False and GC_Repeats(Pairs[1]) == False:
                        List6.append(Pairs)
        List6 = sorted(List6, key=itemgetter(5))
        return List6
        
                
def Conserved_Amplicon(input_fasta, gene_name, number):
        """Takes in input multi-fasta of gene variants and makes a list of possible primers"""
        Mer_List =  Conserved_N_List(input_fasta, number)
        Genes = SeqIO.to_dict(SeqIO.parse(input_fasta, 'fasta'))
        Gene = str(Genes[gene_name].seq)
        List1 = []
        for mers in Mer_List:
                for mers_2 in Mer_List:
                        Length = Amplicon_RF_Length(Gene, mers, mers_2)
                        if Length <= 150 and Length >= 75:
                                Pair = []
                                Pair.append(mers)
                                Pair.append(mers_2)
                                Pair.append(Length)
                                List1.append(Pair)
        List2 = []
        for Pairs in List1:
                Amplicon = Amplicon_RF_Finder(Gene, Pairs[0], Pairs[1])
                if (GC_Content(Amplicon) >= 0.5 and GC_Content(Amplicon) <= 0.6) and Long_Repeats == False:
                        List2.append(Pairs)
        return List2
        

def In_List(item, list1):
    """determines if item is in list1"""
    for x in list1:
        if x == item:
            return True
    else: return False

def extend_ambiguous_dna(seq):
   """return list of all possible sequences given an ambiguous DNA input"""
   d = Seq.IUPAC.IUPACData.ambiguous_dna_values
   r = []
   for i in product(*[d[j] for j in seq]):
      r.append("".join(i))
   return r 
                
##def Ambiguous_Amplicon_Assembly(assembly, primer_f, primer_r):
##        """Finds contig/gene of amplicon with ambiguous base primers"""
##        primer_f_list = extend_ambiguous_dna(primer_f)
##        primer_rev_seq = SeqRecord(Seq(primer_r, IUPACAmbiguousDNA()))
##        primer_rev = str(primer_rev_seq.seq.reverse_complement())
##        primer_r_list = extend_ambiguous_dna(primer_rev)


def Amplicon_Finder_Assembly(assembly, primer_f, primer_r):
        """Takes in an assembly and returns the contig/gene name"""
        Genome = list(SeqIO.parse(assembly, 'fasta'))
        primer_rev_seq = SeqRecord(Seq(primer_r, generic_dna))
        primer_rev = str(primer_rev_seq.seq.reverse_complement())
        for gene in Genome:
                if (primer_f in str(gene.seq)) == False:
                        continue
                elif (primer_rev in str(gene.seq)) == False:
                        continue
                else:
                        return gene.id
        else:
                return "No matches!"    
