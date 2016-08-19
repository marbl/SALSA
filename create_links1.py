import operator
import math
import sys
import argparse

contig_lengths = {}
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-m","--mapping", help="mapping of read to contigs in bam format")
    parser.add_argument("-c",'--counts',help="RE counts")
    parser.add_argument("-l",'--contiglen',help="file containing length of contigs")

    args = parser.parse_args()

    with open(args.contiglen,'r') as cfile:
        lines = cfile.readlines()
        for line in lines:
            attrs = line.split()
            contig_lengths[attrs[0]] = float(attrs[1])


    coordinates_map = {}
    first_in_pair = {}
    second_in_pair = {}
    print >> sys.stderr, 'bedfile started'
    with open(args.mapping,'r') as f:
        for line in f:
            attrs = line.split()
            try:
                pos = (long(attrs[1]) + long(attrs[2]))/2.0
                read = attrs[3]
                if read[-1] == '1':
                    read = read[:-2]
                    rec = (attrs[0],pos)
                    first_in_pair[read] = rec
                else:
                    read = read[:-2]
                    rec = (attrs[0],pos)
                    second_in_pair[read] = rec
                #print >> sys.stderr, len(first_in_pair)
            except:
                continue


    print >> sys.stderr, 'bedfile loaded'
    RF_counts_left = {}
    RF_counts_right = {}

    with open(args.counts,'r') as f:
        #lines = f.readlines()
        for line in f:
            attrs = line.split()
            RF_counts_left[attrs[0]] = int(attrs[1])
            RF_counts_right[attrs[0]] = int(attrs[2])

    contig_links = {}
    exp_links = {}
    norm_score = {}

    for read in first_in_pair:
        if read in second_in_pair:
            rec1 = first_in_pair[read]
            rec2 = second_in_pair[read]
            if rec1[0] != rec2[0]:
                print >> sys.stderr, 'here'
                len1 = contig_lengths[rec1[0]]
                len2 = contig_lengths[rec2[0]]
                pos1 = rec1[1]
                pos2 = rec2[1]
                key = ''
                key1 = {}
                tot_len = len1/2 + len2/2
                if pos1 <= len1/2 and pos2 <= len2/2:
                    key = rec1[0] + ":B$" + rec2[0] + ':B'
                    norm_score[key] = RF_counts_left[rec1[0]] + RF_counts_left[rec2[0]]
                if pos1 <= len1/2 and pos2 > len2/2:
                    key = rec1[0] + ":B$" + rec2[0] + ':E'
                    norm_score[key] = RF_counts_left[rec1[0]] + RF_counts_right[rec2[0]]
                    #exp_links[key] = exp_num
                if pos1 > len1/2 and pos2 <= len2/2:
                    key = rec1[0] + ":E$" + rec2[0] + ':B'
                    norm_score[key] = RF_counts_right[rec1[0]] + RF_counts_left[rec2[0]]
                    #exp_links[key] = exp_num
                if pos1 > len1/2 and pos2 > len2/2:
                    key = rec1[0] + ":E$" + rec2[0] + ':E'
                    norm_score[key] = RF_counts_right[rec1[0]] + RF_counts_right[rec2[0]]
                    #exp_links[key] = exp_num
                if key not in contig_links:
                    contig_links[key] = 1
                contig_links[key] += 1

    #print len(contig_links)

    for key in contig_links:
        edge = key.split('$')
        #rf_sum = RF_counts[edge[0].split(':')[0]] + RF_counts[edge[1].split(':')[0]]
        if norm_score[key] == 0:
            score = 0
        else:
            score = contig_links[key]*1.0/norm_score[key]
        print edge[0],edge[1],score, contig_links[key]

if __name__ == '__main__':
    main()
