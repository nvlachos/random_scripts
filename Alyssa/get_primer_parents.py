#!/usr/bin/python

# USAGE python get_primer_parents.py -i reference.fasta -p primer.fasta -o outfile.txt

from Bio import SeqIO
from Bio import Seq
from argparse import ArgumentParser
import sys

def parse_args():
	parser=ArgumentParser(description='Get Primer affiliations from a reference')
	parser.add_argument('-i',dest='infile',required=True,help = 'reference fasta')
	parser.add_argument('-p',dest='primers', required=True,help='primer fasta')
	parser.add_argument('-o',dest='outfile', required=True, help='outfile')
	return parser.parse_args()

def main():
	args = parse_args()
	primerdict = {}
	#for primer in primers output if it is in the database
	with open(args.outfile,'w') as outfile:
		for rec in SeqIO.parse(args.primers,'fasta'):
			seq = rec.seq
			primerdict[rec.id] = []
			for refseq in SeqIO.parse(args.infile,'fasta'):
				if seq in refseq.seq:
					primerdict[rec.id].append(refseq.id)
				if seq.reverse_complement() in refseq.seq:
					primerdict[rec.id].append(refseq.id)
			pp = ','.join(primerdict[rec.id])
			if pp == '':
				pp = 'NA'
			primer_name='_'.join(rec.id.split('_')[:-1])
			outfile.write('{0}\t{1}\t{2}\n'.format(rec.id, primer_name, pp))




if __name__ == '__main__':
	sys.exit(main())
