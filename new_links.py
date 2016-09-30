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
    parser.add_argument("-c",'--counts',help="RE counts")
    parser.add_argument("-l",'--contiglen',help="file containing length of contigs")

    args = parser.parse_args()

    with open(args.contiglen,'r') as cfile:
        lines = cfile.readlines()
        for line in lines:
            attrs = line.split()
            contig_lengths[attrs[0]] = float(attrs[1])

    with open(args.counts,'r') as f:
    #lines = f.readlines()
        for line in f:
            attrs = line.split()
            RF_counts_left[attrs[0]] = int(attrs[1])
            RF_counts_right[attrs[0]] = int(attrs[2])

    contig_links = {}
    exp_links = {}
    norm_score = {}


    coordinates_map = {}
    first_in_pair = {}
    second_in_pair = {}
    print >> sys.stderr, 'bedfile started'
    prev_line = ''
    with open(args.mapping,'r') as f:
        for line in f:
            attrs = line.split()
            try:
                if prev_line == '':
                    prev_line = line
                    continue

                prev_attrs = prev_line.split()
                if prev_attrs[0] != attrs[0]:
                    prev_read = prev_attrs[3]
                    curr_read = attrs[3]
                    #print prev_read.split('/')[0]
                    if prev_read.split('/')[0] == curr_read.split('/')[0]:
                        #print 'here'
                        pos1 = (long(prev_attrs[1]) + long(prev_attrs[2]))/2.0
                        pos2 = (long(attrs[1]) + long(attrs[2]))/2.0
                        l1 = 300000.0
                        l2 = 300000.0
                       
                        first = ''
                        second = ''
                        if prev_attrs[0] < attrs[0]:
                            first = prev_attrs[0]
                            second = attrs[0]
                        else:
                            first = attrs[0]
                            second = prev_attrs[0]
                            pos1,pos2 = pos2, pos1

                        len1 = contig_lengths[first]
                        len2 = contig_lengths[second]
			#l1 = len1/4
			#l2 = len2/4
                        r1 = l1 / len1
                        r2 = l2 / len2
                        if pos1 <= l1 and pos2 <= l2:
                            key = first + ":B$" + second + ':B'
                            #key1 = rec2[0] + ":B$" + rec1[0] + ':B'
                            norm_score[key] = RF_counts_left[first]*r1 + RF_counts_left[second]*r2
                            #norm_score[key1] = RF_counts_left[rec1[0]] + RF_counts_left[rec2[0]]
                        if pos1 <= l1 and pos2 > len2 - l2:
                            key = first + ":B$" + second + ':E'
                            #key1 = rec2[0] + ":E$" + rec1[0] + ':B'
                            norm_score[key] = RF_counts_left[first]*r1 + RF_counts_right[second]*r2
                            #norm_score[key1] = RF_counts_left[rec1[0]] + RF_counts_left[rec2[0]]
                            #exp_links[key] = exp_num
                        if pos1 > len1 - l1 and pos2 <= l2:
                            key = first + ":E$" + second + ':B'
                            #key1 = rec2[0] + ":B$" + rec1[0] + ':E'
                            norm_score[key] = RF_counts_right[first]*r1 + RF_counts_left[second]*r2
                            #norm_score[key1] = RF_counts_left[rec1[0]] + RF_counts_left[rec2[0]]
                            #exp_links[key] = exp_num
                        if pos1 > len1 - l1 and pos2 > len2 - l2:
                            key = first + ":E$" + second + ':E'
                            #key1 = rec2[0] + ":E$" + rec1[0] + ':E'
                            norm_score[key] = RF_counts_right[first]*r1 + RF_counts_right[second]*r2
                            #norm_score[key1] = RF_counts_left[rec1[0]] + RF_counts_left[rec2[0]]
                            #exp_links[key] = exp_num
                        if key not in contig_links and key != '':
                            contig_links[key] = 0
                        if key != '':
                            contig_links[key] += 1

                prev_line = line
                        

                #print >> sys.stderr, len(first_in_pair)
            except:
                continue


    print >> sys.stderr, 'bedfile loaded'
  


    #print len(contig_links)

    for key in contig_links:
        if key != '':
            edge = key.split('$')
            #rf_sum = RF_counts[edge[0].split(':')[0]] + RF_counts[edge[1].split(':')[0]]
            if norm_score[key] == 0:
                score = 0
            else:
                score = contig_links[key]*1.0/norm_score[key]
            print edge[0],edge[1],score, contig_links[key]

if __name__ == '__main__':
    main()
