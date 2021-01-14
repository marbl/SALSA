#!/usr/bin/env python2
import argparse,pickle

parser = argparse.ArgumentParser()
parser.add_argument('-b','--bed',help="unitig to contig bed file")
parser.add_argument('-c','--contigs',help="contigs fasta file")
parser.add_argument('-u','--unitigs',help="unitigs fasta file")
parser.add_argument('-p','--pickle',help='pickle file for scaffolds')
parser.add_argument('-o','--output',help='output file to write scaffolds')

args = parser.parse_args()

revcompl = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A','N':'N','R':'N','M':'N','Y':'N','S':'N','W':'N','K':'N','a':'t','c':'g','g':'c','t':'a','n':'n',' ':'',}[B] for B in x][::-1])

def parse_fasta(fh):
    fa = {}
    current_short_name = None
    for ln in fh:
        if ln[0] == '>':
            long_name = ln[1:].rstrip()
            current_short_name = long_name.split()[0]
            fa[current_short_name] = []
        else:
            fa[current_short_name].append(ln.rstrip())
    for short_name, nuc_list in fa.iteritems():
        fa[short_name] = ''.join(nuc_list)
    return fa


scaffolds = pickle.load(open(args.pickle,'r'))

contig_seqs = parse_fasta(open(args.contigs,'r'))
unitig_seqs = parse_fasta(open(args.unitigs,'r'))


ofile = open(args.output,'w')

scaffold_id = 1
'''
Load the bed file and store the unitig to contig mapping
'''

unitig2contig = {}
contig2path = {}
contig2unitig = {}
unitig2coord = {}
with open(args.bed,'r') as f:
    for line in f:
        attrs = line.split()
        #contig = 'tig'+attrs[0][3:]
        contig = attrs[0]
        unitig = attrs[3]
        #unitig = 'tig'+attrs[3][3:]
        unitig2contig[unitig] = contig
        if contig not in contig2path:
            contig2path[contig] = []
        contig2path[contig].append(unitig)
        unitig2coord[unitig] = [int(attrs[1]),int(attrs[2])]
        contig2unitig[contig] = unitig

'''
Now, iterate through the scaffold paths to see if we find maximal contig paths in a merge sort style
'''

for scaffold in scaffolds:
    path = scaffolds[scaffold]
    contig_start = 0
    contig_end = 0
    path_normal = []
    path_with_ori = []
    for i in xrange(0,len(path),2):
       path_normal.append(path[i].split(':')[0])
       if path[i].split(':')[1] + path[i+1].split(':')[1] == 'BE':
        path_with_ori.append(path[i].split(':')[0]+':FOW')
       else:
        path_with_ori.append(path[i].split(':')[0]+':REV')

    '''
    Now check in merge sort style
    '''
    new_path = []
    i = 0
    if len(path_normal) == 1:
        new_path = path_with_ori
        continue
    while i < len(path_normal):
        curr_utg = path_normal[i]
        #print curr_utg
        if curr_utg not in unitig2contig:
            new_path.append(path_with_ori[i])
            i += 1
        else:
            '''
            Check for forward orientation of unitig to contig tiling
            '''
            contig_path = contig2path[unitig2contig[curr_utg]]
            if len(contig_path) == 1:
                i += 1
                new_path.append(path_with_ori[i-1])
                continue
            utg_index_f = contig_path.index(curr_utg)
            curr_ind_f = i+1
            utg_index_f += 1
            to_check_f = True

            if utg_index_f >= len(contig_path):
                to_check_f = False
            if curr_ind_f >= len(path_normal):
                new_path.append(path_with_ori[i])
                break

            start_coord_f = unitig2coord[curr_utg][0]
            end_coord_f = unitig2coord[curr_utg][1]

            fcount = 1
            if to_check_f:
                while (contig_path[utg_index_f] == path_normal[curr_ind_f]):
                    end_coord = unitig2coord[contig_path[utg_index_f]][1]
                    utg_index_f += 1
                    curr_ind_f += 1
                    fcount += 1
                    if utg_index_f >= len(contig_path) or curr_ind_f >= len(path_normal):
                        break

                if fcount > 1:
                    new_path.append([unitig2contig[curr_utg],start_coord_f,end_coord_f,'FOW'])
                    i = curr_ind_f
                    continue

            rpath = contig_path[::-1]
            rcount = 1
            utg_index_r = rpath.index(curr_utg)
            #print 'curr utg = ' + str(curr_utg) + ' ' + str(utg_index_r)
            curr_ind_r = i + 1
            utg_index_r += 1

            #print 'fpath = ' + str(contig_path)
            #print 'rpath = ' + str(rpath)
            if utg_index_r >= len(contig_path):
                new_path.append(path_with_ori[i])
                i += 1
                continue
            if curr_ind_r >= len(path_normal):
                new_path.append(path_with_ori[i])
                break
            #print 'prev scaf = ' + str(rpath[utg_index_r - 1])
            #print 'curr scaf = ' + str(rpath[utg_index_r])
            #print 'prev tiling = ' + str(path_normal[curr_ind_r - 1])
            #print 'curr tiling = ' + str(path_normal[curr_ind_r])
            start_coord_r = unitig2coord[curr_utg][1]
            end_coord_r = unitig2coord[curr_utg][0]
            while (rpath[utg_index_r] == path_normal[curr_ind_r]):
                end_coord_r = unitig2coord[rpath[utg_index_r]][0]
                utg_index_r += 1
                curr_ind_r += 1
                rcount += 1
                if utg_index_r >= len(rpath)  or curr_ind_r >= len(path_normal):
                    break

            if rcount > 1:
                new_path.append([unitig2contig[curr_utg],start_coord_r,end_coord_r,'REV'])
                i = curr_ind_r
            else:
                new_path.append(path_with_ori[i])
                i += 1

    scaffold_seq = ''
    gap = ''
    for i in xrange(500):
        gap += 'N'
    for i in xrange(len(new_path)):
        if len(new_path[i]) != 4:
            utg,ori = new_path[i].split(':')
            if ori == 'FOW':
                scaffold_seq += unitig_seqs[utg]
            else:
                scaffold_seq += revcompl(unitig_seqs[utg])

            if i != len(new_path) - 1:
                scaffold_seq += gap
        else:
            contig = new_path[i][0]
            start = int(new_path[i][1])
            end = int(new_path[i][2])
            if start < end:
                scaffold_seq += contig_seqs[contig][start:end]
            else:
                scaffold_seq += revcompl(contig_seqs[contig][end:start])

            if i != len(new_path) - 1:
                scaffold_seq += gap
    ofile.write('>scaffold_'+str(scaffold_id)+'\n')
    scaffold_id += 1
    chunks = [scaffold_seq[i:i+80] for i in xrange(0,len(scaffold_seq),80)]
    for each in chunks:
        ofile.write(each+'\n')

    print "=============================="
    print path_with_ori
    print new_path
    print "============================="




