import sys

#reads input assembly, breakpoints given by the method and outputs new contig file with lengths
#offsets bed file as well

def parse_fasta(fh):
    fa = {}
    current_short_name = None
    # Part 1: compile list of lines per sequence
    for ln in fh:
        if ln[0] == '>':
            # new name line; remember current sequence's short name
            long_name = ln[1:].rstrip()
            current_short_name = long_name.split()[0]
            fa[current_short_name] = []
        else:
            # append nucleotides to current sequence
            fa[current_short_name].append(ln.rstrip())
    # Part 2: join lists into strings
    for short_name, nuc_list in fa.iteritems():
        # join this sequence's lines into one long string
        fa[short_name] = ''.join(nuc_list)
    return fa



#read fasta first
input_seqs = parse_fasta(open(sys.argv[1],'r'))

#read breakpoints

contig2breakpoints = {}
with open(sys.argv[2],'r') as f:
    for line in f:
        attrs = line.split()
        contig2breakpoints[attrs[0]] = int(attrs[1])



#first breake the input assembly and store the mapping of old to new names of contigs
contig2new = {}
contig2newseq = {}
contig_id = 1
for seq in input_seqs:
    if seq not in contig2breakpoints:
        contig2new[seq] = seq
        contig2newseq[seq] = input_seqs[seq]
        contig_id += 1
    else:
        first = input_seqs[seq][:contig2breakpoints[seq]]
        second = input_seqs[seq][contig2breakpoints[seq]:]
        first_id = seq+'_1'
        contig_id += 1
        second_id = seq+'_2'
        contig_id += 1
        contig2new[seq] = [first_id,second_id]
        contig2newseq[first_id] = first
        contig2newseq[second_id] = second


#now update the bed file in streaming fashion
oline = ""
count = 0
ofile = open(sys.argv[4]+'/alignment_iteration_1.tmp.bed','w')
with open(sys.argv[3],'r') as f:
    for line in f:
        attrs = line.split()
        if attrs[0] not in contig2breakpoints:
            oline += str(contig2new[attrs[0]] +'\t'+attrs[1]+'\t'+attrs[2]+'\t'+attrs[3]+'\n')
            count += 1
        else:
            pos1 = int(attrs[1])
            pos2 = int(attrs[2])
            breakpoint_loc = contig2breakpoints[attrs[0]]
            if pos1 < breakpoint_loc and pos2 > breakpoint_loc:
                continue
            else:
                if pos2 < breakpoint_loc:
                    oline += str(contig2new[attrs[0]][0]+'\t'+attrs[1]+'\t'+attrs[2]+'\t'+attrs[3]+'\n')
                    count += 1
                else:
                    new_start = pos1 - breakpoint_loc
                    new_end = pos2 - breakpoint_loc
                    oline += str(contig2new[attrs[0]][1]+'\t'+str(new_start)+'\t'+str(new_end)+'\t'+attrs[3]+'\n')
                    count += 1
        if count >= 1000000:
            ofile.write(oline)
            oline = ""
            count = 0

ofile.close()
#write fasta file
ofasta = open(sys.argv[4]+'/asm.cleaned.fasta','w')
for seq in contig2newseq:
    contig_seq = contig2newseq[seq]
    chunks = [contig_seq[i:i+80] for i in xrange(0,len(contig_seq),80)]
    ofasta.write('>'+seq+'\n')
    for chunk in chunks:
        ofasta.write(chunk+'\n')


ofasta.close()

#write lengths
olens = open(sys.argv[4]+'/scaffold_length_iteration_1','w')
for seq in contig2newseq:
    olens.write(seq+'\t'+str(len(contig2newseq[seq]))+'\n')

olens.close()

