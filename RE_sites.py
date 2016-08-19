from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

import re 
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-a","--assembly",help="assembled contigs")
    #parser.add_argument("-m","--mapping", help="mapping of read to contigs in bam format")
    #parser.add_argument("-d","--dir",help="output directory for results",default='out')
    args = parser.parse_args()
    RF = 'AAGCTT'
    for record in SeqIO.parse(args.assembly, "fasta"):
        
        id,seq = record.id, str(record.seq)

        pos  = [m.start(0) for m in re.finditer(RF,seq)]
     
        length = len(seq)

        left_count = 0
        rigt_count = 0
        for each in pos:
        	if each < length/2:
        		left_count += 1
        	else:
        		rigt_count += 1

        print id, left_count, rigt_count

    

if __name__ == '__main__':
    main()
