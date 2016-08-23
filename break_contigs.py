from pbcore.io import FastaReader
from pbcore.io import FastaWriter
import argparse

def main():
	id2seq = {}
	parser = argparse.ArgumentParser()
	parser.add_argument("-b","--breakpoint",help="file containing breakpoints")
	parser.add_argument("-a","--assembly",help="fasta file containing contigs")
	parser.add_argument("-o","--outfile",help="new assembly file")
	parser.add_argument("-l","--lenfile",help="length of contigs")

	args = parser.parse_args()

	lenfile = open(args.lenfile,'w')

	lenmap = {}
	f = FastaReader(args.assembly)
    for record in f:
        id = record.id
        id2seq[id] = record.sequence[0:-10]


	new_seq = {}
	f = open(args.breakpoint,'r')
	lines = f.readlines()
	for line in lines:
		attrs = line.split()
		if len(attrs) == 1:
			curr_contig = attrs[0]
			seq = id2seq[curr_contig]
		else:
			start = long(attrs[0])
			end = long(attrs[1])
			new_id = curr_contig + '_' + attrs[0] + '_' + attrs[1]
			new_seq[new_id] = seq[start:end]
			lenmap[new_id] = end - start + 1

	rec_list = []
	writer = FastaWriter(args.scaffold)
	for key in new_seq:
		writer.writeRecord(key,new_seq[key])

	for key in lenmap:
		lenfile.write(key + "\t" + str(lenmap[key]) + '\n')

if __name__ == '__main__':
	main()
