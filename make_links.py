#!/usr/bin/env python2
import operator
import math
import sys
import argparse

contig_lengths = {}
RF_counts_left = {}
RF_counts_right = {}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-b","--mapping", help="mapping of read to contigs in bed format")
    parser.add_argument("-i","--iteration", help="iteration number")
    parser.add_argument("-d",'--directory',help="output directory")
    parser.add_argument("-x","--dup",help="Duplicated links")

    args = parser.parse_args()

    iteration = str(args.iteration)
    dup_links = {}
    if args.dup != 'abc' and iteration  == 1:
        with open(args.dup,'r') as f:
            for line in f:
                attrs = line.split()
                dup_links[attrs[0]+'$'+attrs[1]] = 1
                dup_links[attrs[1]+'$'+attrs[0]] = 1

    avoid_links = {}
    if int(iteration) > 1:
        with open(args.directory+'/links_avoid_iteration_'+str(iteration),'r') as f:
            for line in f:
                attrs = line.split()
                avoid_links[attrs[0]+'$'+attrs[1]] = 1
                avoid_links[attrs[1]+'$'+attrs[0]] = 1

    with open(args.directory+'/scaffold_length_iteration_'+str(iteration),'r') as cfile:
        lines = cfile.readlines()
        for line in lines:
            attrs = line.split()
            contig_lengths[attrs[0]] = float(attrs[1])

    with open(args.directory+'/re_counts_iteration_'+str(iteration),'r') as f:
    #lines = f.readlines()
        for line in f:
            attrs = line.split()
            RF_counts_left[attrs[0]] = int(attrs[1])
            RF_counts_right[attrs[0]] = int(attrs[2])

    contig_links = {}
    exp_links = {}
    norm_score = {}
    all_counts = {}

    coordinates_map = {}
    first_in_pair = {}
    second_in_pair = {}
    print >> sys.stderr, 'bedfile started'
    prev_line = ''
    with open(args.mapping,'r') as f:
        for line in f:
            attrs = line.split()
            # if len(attrs) < 6:
            #     continue
            # try:
            if prev_line == '':
                prev_line = line
                continue

            prev_attrs = prev_line.split()
            if prev_attrs[0] != attrs[0]:
                ks = prev_attrs[0] + '$' + attrs[0]
                if ks in dup_links or ks in avoid_links:
                    continue
                prev_read = prev_attrs[3]
                curr_read = attrs[3]
                #print prev_read.split('/')[0]
                if prev_read.split('/')[0] == curr_read.split('/')[0]:
                    #print 'here'

                    pos1 = (long(prev_attrs[1]) + long(prev_attrs[2]))/2.0
                    pos2 = (long(attrs[1]) + long(attrs[2]))/2.0
                    # l1 = 300000.0
                    # l2 = 300000.0
                    #print pos1, pos2
                    key = ''
                    #print key
                    first = ''
                    second = ''
                    if prev_attrs[0] < attrs[0]:
                        first = prev_attrs[0]
                        second = attrs[0]
                    else:
                        first = attrs[0]
                        second = prev_attrs[0]
                        pos1,pos2 = pos2, pos1
                    t_key = first+'$'+second
                    if t_key not in all_counts:
                        all_counts[t_key] = 0
                    all_counts[t_key] += 1
                    len1 = contig_lengths[first]
                    len2 = contig_lengths[second]
                    #print len1,len2
                    l1 = len1/2
                    l2 = len2/2
                    r1 = l1 / len1
                    r2 = l2 / len2
                    if pos1 <= l1 and pos2 <= l2:
                        #print 'here1'
                        key = first + ":B$" + second + ':B'
                        #key1 = rec2[0] + ":B$" + rec1[0] + ':B'
                        #norm_score[key] = RF_counts_left[first]*r1 + RF_counts_left[second]*r2
                        norm_score[key] = RF_counts_left[first]*r1 + RF_counts_left[second]*r2
                        #norm_score[key1] = RF_counts_left[rec1[0]] + RF_counts_left[rec2[0]]
                    if pos1 <= l1 and pos2 > len2 - l2:
                        #print 'here2'
                        key = first + ":B$" + second + ':E'
                        #key1 = rec2[0] + ":E$" + rec1[0] + ':B'
                        #norm_score[key] = RF_counts_left[first]*r1 + RF_counts_right[second]*r2
                        norm_score[key] = RF_counts_left[first]*r1 + RF_counts_right[second]*r2
                        #norm_score[key1] = RF_counts_left[rec1[0]] + RF_counts_left[rec2[0]]
                        #exp_links[key] = exp_num
                    if pos1 > len1 - l1 and pos2 <= l2:
                        #print 'here3'
                        key = first + ":E$" + second + ':B'
                        #print key
                        #key1 = rec2[0] + ":B$" + rec1[0] + ':E'
                        #norm_score[key] = RF_counts_right[first]*r1 + RF_counts_left[second]*r2
                        norm_score[key] = RF_counts_right[first]*r1 + RF_counts_left[second]*r2
                        #norm_score[key1] = RF_counts_left[rec1[0]] + RF_counts_left[rec2[0]]
                        #exp_links[key] = exp_num
                    if pos1 > len1 - l1 and pos2 > len2 - l2:
                        #print 'here4'
                        key = first + ":E$" + second + ':E'
                        #key1 = rec2[0] + ":E$" + rec1[0] + ':E'
                        #norm_score[key] = RF_counts_right[first]*r1 + RF_counts_right[second]*r2
                        norm_score[key] = RF_counts_right[first]*r1 + RF_counts_right[second]*r2
                        #norm_score[key1] = RF_counts_left[rec1[0]] + RF_counts_left[rec2[0]]
                        #exp_links[key] = exp_num
                    if key not in contig_links and key != '':
                        contig_links[key] = 0
                    if key != '':
                        #print 'adding'
                        contig_links[key] += 1

            prev_line = line


                #print >> sys.stderr, len(first_in_pair)
            # except:
            #     print 'here'
            #     continue


    print >> sys.stderr, 'bedfile loaded'



    print len(contig_links)
    ofile = open(args.directory+'/contig_links_iteration_'+str(iteration),'w')
    for key in contig_links:
        if key != '':
            edge = key.split('$')
            #rf_sum = RF_counts[edge[0].split(':')[0]] + RF_counts[edge[1].split(':')[0]]
            if norm_score[key] == 0:
                score = 0
            else:
                score = contig_links[key]*1.0/norm_score[key]
            ofile.write(edge[0]+'\t'+edge[1]+'\t'+str(score)+'\t'+str(contig_links[key])+'\n')

    #write all count links
    ofile.close()

    # ofile =  open(args.directory+'/contig_links_raw_iteration_1','w')
    # for key in all_counts:
    #     contigs = key.split('$')
    #     ofile.write(contigs[0]+'\t'+contigs[1]+'\t'+str(all_counts[key])+'\n')

if __name__ == '__main__':
    main()
