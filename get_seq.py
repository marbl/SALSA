#!/usr/bin/env python2
import sys
import re
import pickle
import argparse

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

# return value is a dictionary. keys=str, values=int.
# the keys are tig names and the values are break points.
# It is the result of parsing the input_breaks file generated
# from cleaning the input assembly to SALSA using break_contigs_start.cpp.
# input param: fn is the filename (a string) of the input_breaks file.
def parse_breaks(fn):
    if fn is None:
        return None
    breaks = {}
    with open(fn, 'r') as fh:
        for line in fh:
            fields = line.rstrip('\n').split('\t')
            tigname = fields[0]          # original tig name (i.e., before _1 & _2 were appended)
            break_point = int(fields[1]) # first position of next segment (0-based). Alternately, last position of first segment (1-based, inclusive).
            assert not tigname in breaks # if a contig is broken multiple times (reported in input_breaks 2+ times), we need to handle things differently here
            #                            # (e.g., by using a list of positions instead of a single position), which will have implications elsewhere too (obviously)
            breaks[tigname] = break_point
    return breaks

# returns two values, a string and an integer.
# input params: name is the tig names from the cleaned fasta (i.e., not scaffold_* (unless the input to salsa was also named that way))
#               position is the start position we're about to write to an agp file
#               breaks_map is the dictionary generated from the parse_breaks function (cannot be None)
def convertToPreCleanedNamesAndCoords(name, position, breaks_map):
    assert not breaks_map is None # If breaks_map were None, we wouldn't have the necessary information to possibly update name and position

    # initialize the return values
    orig_name = name
    orig_pos = position

    if not re.search(r"_[12]$", name) is None:
        # name pattern indicates it may be a broken/cleaned tig
        possible_orig_name = name[:-2]

        if possible_orig_name in breaks_map:
            # the un-suffixed name was in breaks_map as a key, so it is a broken/cleaned tig

            orig_name = possible_orig_name # save the actual original name
            break_pos = breaks_map[orig_name] # extract the break point

            if not re.search(r"_[2]$", name) is None:
                # the name pattern indicates that we're dealing with the part of the tig after
                # the break point. If it had been the first portion (i.e., _1), we wouldn't
                # update the position because it would already be in the correct coordinate space.
                # Since we aren't in the first portion, we do need to add the length of the first
                # portion to our position so we are in the correct coordinate space.
                orig_pos +=  break_pos

    # return the (possibly updated) name and position
    return orig_name, orig_pos

parser = argparse.ArgumentParser()
parser.add_argument("-a","--cleaned",help="cleaned assembly (input, required)")
parser.add_argument("-b","--breaks", help="input_breaks file resulting from cleaning input assembly. The -G option must be supplied for this to work. (input, optional)")
parser.add_argument("-f","--scaffold",help="final scaffold file (output, required)")
parser.add_argument("-g","--agp",help="agp file (output, required)")
parser.add_argument("-G","--agp-orig-coords", dest="agp_oc", help="agp file with the original names and coordinates from the input assembly before cleaning. The -b option must be supplied for this to work. (output, optional)")
parser.add_argument("-p","--map",help="pickle map of scaffolds (input, required)")

args = parser.parse_args()

revcompl = lambda x: ''.join([{'A':'T','B':'N','C':'G','G':'C','T':'A','N':'N','R':'N','M':'N','Y':'N','S':'N','W':'N','K':'N','a':'t','c':'g','g':'c','t':'a','n':'n',' ':'',}[B] for B in x][::-1])
scaff_map = pickle.load(open(args.map,'r'))

contig_length = {}
id2seq = {}

if bool(args.breaks) ^ bool(args.agp_oc):
    sys.stderr.write("ERROR: Both -b, --breaks and -G, --agp-orig-coords must be supplied together. None or both are fine, but not one only.\n")
    sys.exit(1)
write_agp_orig_coords = bool(args.agp_oc) # true if we should write the special agp file, false otherwise.
breaks = parse_breaks(args.breaks)

id2seq = parse_fasta(open(args.cleaned,'r'))
for key in id2seq:
    id2seq[key] = id2seq[key]
    contig_length[key] = len(id2seq[key])

#first sort scaffolds in decreasing order of length

scaff2length = {}
for scaffold in scaff_map:
    path = scaff_map[scaffold]
    length = 0
    for i in xrange(0,len(path)-1,2):
        length += contig_length[path[i].split(':')[0]]
    scaff2length[scaffold] = length

sorted_scaffolds = sorted(scaff2length.items(), key=lambda x: x[1],reverse=True)

c_id = 1
with open(args.agp, 'w') as agp_output:
    with open(args.scaffold, 'w') as ofile:
        agp_origcoords_output = open(args.agp_oc, 'w') if write_agp_orig_coords else None
        for key in sorted_scaffolds:
            #print 'scaffold_'+str(c_id) + '\t' + key
            key = key[0]
            start = 1
            local_comp = 1
            #if len(scaff_map[key]) >= 4:
            path = scaff_map[key]
            scaff_len = 0
            curr_contig = ""
            #print c_id
            for i in range(0,len(path)-1,2):

                curr = path[i].split(':')
                next = path[i+1].split(':')
                curr_len = contig_length[curr[0]]
                scaff_len += curr_len
                end = curr_len + start - 1
                direction = ''
                if curr[1] == 'B' and next[1] == 'E':
                    curr_contig += id2seq[curr[0]]
                    direction = '+'
                if curr[1] == 'E' and next[1] == 'B':
                    curr_contig += revcompl(id2seq[curr[0]])
                    direction = '-'
                #             0............................................................................................0  1.....1  2.2  3...........3  4..............4
                line_parts = ["scaffold_" + str(c_id) + '\t' + str(start) + '\t' + str(end) + '\t' + str(local_comp) + "\tW", curr[0], '1', str(curr_len), direction + '\n']
                agp_output.write('\t'.join(line_parts))
                if write_agp_orig_coords:
                    orig_name, orig_start = convertToPreCleanedNamesAndCoords(curr[0], 1, breaks)
                    orig_end = curr_len + orig_start - 1
                    line_parts[1] = orig_name
                    line_parts[2] = str(orig_start)
                    line_parts[3] = str(orig_end)
                    agp_origcoords_output.write('\t'.join(line_parts))
                    # note: line_parts is unsafe to use after this if block w/o reassigning the start/end and name:
                    #line_parts[1] = curr[0]
                    #line_parts[2] = '1'
                    #line_parts[3] = str(curr_len)

                start = end + 1
                local_comp += 1

                if i != len(path) - 2:
                    curr_contig += 'N' * 500

                    end = 500 + start - 1
                    line = "scaffold_" + str(c_id) + '\t' + str(start) + '\t' + str(end) + '\t' + str(local_comp) + '\t' + "N\t500\tscaffold\tyes\tna\n"
                    agp_output.write(line)
                    if write_agp_orig_coords:
                        agp_origcoords_output.write(line)
                    start = end + 1
                    local_comp += 1
            # rec = SeqRecord(Seq(curr_contig,generic_dna),id='scaffold_'+str(c_id))
            # recs.append(rec)
            # print c_id
            chunks = [curr_contig[i:i+80] for i in xrange(0,len(curr_contig),80)]
            ofile.write('>scaffold_'+str(c_id)+'\n')
            for chunk in chunks:
                ofile.write(chunk+'\n')
            c_id += 1

        if write_agp_orig_coords:
            agp_origcoords_output.close()
