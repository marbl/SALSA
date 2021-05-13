#!/usr/bin/env python2
import os,sys
import argparse
import pickle

parser = argparse.ArgumentParser()
parser.add_argument('-b','--bed',help='Original bed file of alignments')
parser.add_argument('-a','--agp',help='AGP file for the final scaffolds')
parser.add_argument('-l','--length',help="Length of input unitigs")

args = parser.parse_args()



scaffolds_current = {}

with open(args.agp,'r') as f:
    for line in f:
        attrs = line.split()
        if attrs[4] == 'N':
            continue
        if attrs[0] not in scaffolds_current:
            scaffolds_current[attrs[0]] = []
        contig = attrs[5]
        if attrs[-1] == '+':
            scaffolds_current[attrs[0]].append(contig+':B')
            scaffolds_current[attrs[0]].append(contig+':E')
        else:
            scaffolds_current[attrs[0]].append(contig+':E')
            scaffolds_current[attrs[0]].append(contig+':B')



#print breakpoints
scaff_id = 1
scaffolds_new = {}
unitig_length = {}
old2new = {}

def update_bed(expanded_scaffold):
    contig2scaffold = {}
    contig2info = {}
    scaffold_length = {}

    #print re_counts

    for key in expanded_scaffold:
        path = expanded_scaffold[key]
        scaffold_length[key] = 0
        offset = 0
        for i in xrange(0,len(path)-1,2):
            contig = path[i].split(':')[0]
            contig2scaffold[contig] = key
            ori = path[i].split(':')[1] + path[i+1].split(':')[1]
            if ori == 'BE':
                contig2info[contig] = (offset,offset+unitig_length[contig],'FOW')
            else:
                contig2info[contig] = (offset,offset+unitig_length[contig],'REV')
            offset += unitig_length[contig]
            scaffold_length[key] += unitig_length[contig]


    o_lines = ""
    count = 0
    prev_line = ""
    #print contig2scaffold
    #sys.exit(1)
    if True:
        #output = open(args.directory+'/alignment_iteration_'+str(iteration)+'.bed','w')
        with open(args.bed,'r') as f:
            olines = ""
            prev_scaffold = ""
            for line in f:
                if prev_line == "":
                    prev_line = line
                    continue
                else:
                    prev_attrs = prev_line.split()
                    curr_attrs = line.split()
                    try:
                        prev_contig = prev_attrs[0]
                        curr_contig = curr_attrs[0]
                        prev_read = prev_attrs[3].split('/')[0]
                        curr_read = curr_attrs[3].split('/')[0]
                        first = prev_attrs[3].split('/')[1]
                        second = curr_attrs[3].split('/')[1]
                    except:
                        continue
                    if prev_contig in contig2scaffold and curr_contig in contig2scaffold:
                        prev_scaffold = contig2scaffold[prev_contig]
                        curr_scaffold = contig2scaffold[curr_contig]
                        if prev_read == curr_read and first == '1' and second == '2':
                       # if prev_read == curr_read and first == '1' and second == '2':
                            prev_info = contig2info[prev_contig]
                            prev_start = int(prev_attrs[1])
                            prev_end = int(prev_attrs[2])
                            if prev_info[2] == 'FOW':
                                new_prev_start = prev_info[0] + prev_start
                                new_prev_end = prev_info[0] + prev_end
                            else:
                                new_prev_start = prev_info[0] + unitig_length[prev_contig] - prev_end
                                new_prev_end = prev_info[0] + unitig_length[prev_contig] - prev_start
                            olines += "0\t"+prev_scaffold+'\t'+str(new_prev_start)+"\t0\t"
                            #o_lines += prev_scaffold+'\t'+str(new_prev_start)+'\t'+str(new_prev_end)+'\t'+prev_attrs[3]+'\n'
                            count += 1

                            curr_info = contig2info[curr_contig]
                            curr_start = int(curr_attrs[1])
                            curr_end = int(curr_attrs[2])
                            if curr_info[2] == 'FOW':
                                new_curr_start = curr_info[0] + curr_start
                                new_curr_end = curr_info[0] + curr_end
                            else:
                                new_curr_start = curr_info[0] + unitig_length[curr_contig] - curr_end
                                new_curr_end = curr_info[0] + unitig_length[curr_contig] - curr_start
                            olines += "1\t"+curr_scaffold+'\t'+str(new_curr_start)+"\t1\n"
                            #o_lines += curr_scaffold+'\t'+str(new_curr_start)+'\t'+str(new_curr_end)+'\t'+curr_attrs[3]+'\n'
                            count += 1
                            if count == 1000000:
                                print olines
                                #output.write(o_lines)
                                count = 0
                                olines = ""

                prev_line = line

            #write remaining lines
            print olines
            #output.write(o_lines)
            #output.close()

with open(args.length,'r') as f:
    for line in f:
        attrs = line.split()
        unitig_length[attrs[0]] = int(attrs[1])


update_bed(scaffolds_current)

