import os,sys
import argparse
import pickle

parser = argparse.ArgumentParser()
parser.add_argument('-i','--iteration',help='Number of current iteration in the pipeline')
parser.add_argument('-d','--directory',help='Project Directory')

args = parser.parse_args()

iteration = int(args.iteration)
breakpoints = {}
scaffolds_current = pickle.load(open(args.directory+'/scaffolds_iteration_'+str(iteration-1)+'.p','r'))

with open(args.directory+'/misasm_iteration_'+args.iteration+'.report','r') as f:
    for line in f:
        attrs = line.strip().split()
#        print attrs
        if line[0] != 's':
            continue
        try:
            if len(attrs) < 2:
                continue
            breakpoints[attrs[0]] = [int(x) for x in attrs[1:]]
        except:
            continue
#print breakpoints
scaff_id = 1
scaffolds_new = {}
unitig_length = {}
old2new = {}

def update_bed(expanded_scaffold):
    contig2scaffold = {}
    contig2info = {}
    scaffold_length = {}

    re_counts = {}
    with open(args.directory+'/re_counts_iteration_1','r') as f:
        for line in f:
            attrs = line.split()
            re_counts[attrs[0]] = (int(attrs[1]),int(attrs[2]))

    #print re_counts

    for key in expanded_scaffold:
        path = expanded_scaffold[key]
        scaffold_length[key] = 0
        offset = 0
        for i in xrange(0,len(path)-1,2):
            contig = path[i].split(':')[0]
            contig2scaffold[contig] = key
            #print path[i],path[i+1]
            ori = path[i].split(':')[1] + path[i+1].split(':')[1]
            if ori == 'BE':
                contig2info[contig] = (offset,'FOW')
            else:
                contig2info[contig] = (offset,'REV')
            offset += unitig_length[contig]
            scaffold_length[key] += unitig_length[contig]

    scaffold_re = {}
    for key in expanded_scaffold:
        path = expanded_scaffold[key]
        length = scaffold_length[key]
        offset = 0
        s_left = 0
        s_right = 0
        if len(path) == 2:
            contig = path[0].split(':')[0]
            scaffold_re[key] = re_counts[contig]
        else:
            for i in xrange(0,len(path),2):
                contig = path[i].split(':')[0]
                contig2scaffold[contig] = key
                left,right = re_counts[contig]
                if offset <= length/2 and i+2 < len(path):
                    if contig2info[path[i+2].split(':')[0]][0] <= length/2:
                        s_left += (left + right)
                    else:
                        contig_next = path[i+2].split(':')[0]
                        if contig2info[contig_next][0] >= length/2:
                            left_part = length/2 - offset
                            right_part = contig2info[path[i+2].split(':')[0]][0] - length/2
                            s_left += left*left_part/unitig_length[contig]
                            s_right += right*right_part/unitig_length[contig]

                if offset <= length/2 and i + 2 >= len(path):
                    left_part = length/2 - offset
                    right_part = length/2
                    s_left += left*left_part/unitig_length[contig]
                    s_right += right*right_part/unitig_length[contig]

                if offset >= length/2:
                    s_right += (left+right)
                offset += unitig_length[contig]
                #scaffold_length[key] += contig_length[contig]
            scaffold_re[key] = (s_left,s_right)

    o_lines = ""
    count = 0
    prev_line = ""
    #print contig2scaffold
    #sys.exit(1)
    if True:
        output = open(args.directory+'/alignment_iteration_'+str(iteration)+'.bed','w')
        with open(args.directory+'/alignment_iteration_1.bed','r') as f:
            for line in f:
                if prev_line == "":
                    prev_line = line
                    continue
                else:
                    prev_attrs = prev_line.split()
                    curr_attrs = line.split()
                    prev_contig = prev_attrs[0]
                    curr_contig = curr_attrs[0]
                    if prev_contig == curr_contig:
                        continue
                    prev_read = prev_attrs[3].split('/')[0]
                    curr_read = curr_attrs[3].split('/')[0]
                    first = prev_attrs[3].split('/')[1]
                    second = curr_attrs[3].split('/')[1]
                    if prev_contig in contig2scaffold and curr_contig in contig2scaffold:
                        prev_scaffold = contig2scaffold[prev_contig]
                        curr_scaffold = contig2scaffold[curr_contig]
                        if prev_read == curr_read and prev_scaffold != curr_scaffold and first == '1' and second == '2':
                       # if prev_read == curr_read and first == '1' and second == '2':
                            prev_info = contig2info[prev_contig]
                            prev_start = int(prev_attrs[1])
                            prev_end = int(prev_attrs[2])
                            new_prev_start = prev_start + prev_info[0]
                            new_prev_end = prev_end + prev_info[0]
                            o_lines += prev_scaffold+'\t'+str(new_prev_start)+'\t'+str(new_prev_end)+'\t'+prev_attrs[3]+'\n'
                            count += 1

                            curr_info = contig2info[curr_contig]
                            curr_start = int(curr_attrs[1])
                            curr_end = int(curr_attrs[2])
                            new_curr_start = curr_start + curr_info[0]
                            new_curr_end = curr_end + curr_info[0]
                            o_lines += curr_scaffold+'\t'+str(new_curr_start)+'\t'+str(new_curr_end)+'\t'+curr_attrs[3]+'\n'
                            count += 1
                            if count == 1000000:
                                output.write(o_lines)
                                count = 0
                                o_lines = ""

                prev_line = line

            #write remaining lines
            output.write(o_lines)
            output.close()
    len_output = open(args.directory+'/scaffold_length_iteration_'+str(iteration),'w')
    for key in scaffold_length:
        len_output.write(key+'\t'+str(scaffold_length[key])+'\n')
    len_output.close()


    re_out = open(args.directory+'/re_counts_iteration_'+str(iteration),'w')
    for key in scaffold_re:
        re_out.write(key+'\t'+str(scaffold_re[key][0])+'\t'+str(scaffold_re[key][1])+'\n')
    re_out.close()

avoid_file = open(args.directory+'/links_avoid_iteration_'+args.iteration,'w')

with open(args.directory+'/scaffold_length_iteration_1','r') as f:
    for line in f:
        attrs = line.split()
        unitig_length[attrs[0]] = int(attrs[1])

for key in scaffolds_current:
    if key not in breakpoints:
        scaffolds_new['scaffold_'+str(scaff_id)] = scaffolds_current[key]
        #old2new[key] = ['scaffold_'+str(scaff_id)]
        scaff_id += 1
    else:
        #print key
        cum_len = 0
        breaks_list = []
        for i in xrange(1,len(scaffolds_current[key]),2):
            utg = scaffolds_current[key][i].split(':')[0]
            cum_len += unitig_length[utg]
            #print cum_len
            if cum_len in breakpoints[key]:
                breaks_list.append(i)
                #TODO: Fix this to consider all breakpoints
        #print breaks_list
        if len(breaks_list) == 1:
            #print 'here'
            first_part = scaffolds_current[key][:breaks_list[0]+1]
            second_part = scaffolds_current[key][breaks_list[0]+1:]
            first_id = 'scaffold_'+str(scaff_id)
            scaff_id += 1
            second_id = 'scaffold_'+str(scaff_id)
            #old2new[key] = [first_id,second_id,break_list[0]]
            scaff_id += 1
            avoid_file.write(first_id+'\t'+second_id+'\n')
            scaffolds_new[first_id] = first_part
            scaffolds_new[second_id] = second_part
        else:
            #print 'here'
            scaffolds_new['scaffold_'+str(scaff_id)] = scaffolds_current[key]
            scaff_id += 1
            #prev = 0
            #prev_scaffold = ''
            #for i in xrange(len(breaks_list)):
            #    if i == len(breaks_list) - 1:
            #        scaff = scaffolds_current[key][i:]
            #        scaffolds_new['scaffold_'+str(scaff_id)] = scaff
            #        avoid_file.write(prev_scaffold+'\t'+'scaffold_'+str(scaff_id)+'\n')
            #        scaff_id += 1
            #        continue
            #    scaff = scaffolds_current[key][prev:i+1]
            #    prev = i + 1
            #    scaffolds_new['scaffold_'+str(scaff_id)] = scaff
            #    if prev_scaffold == '':
            #        prev_scaffold = 'scaffold_'+str(scaff_id)
            #    else:
            #        avoid_file.write(prev_scaffold+'\t'+'scaffold_'+str(scaff_id)+'\n')
            #        prev_scaffold = 'scaffold_'+str(scaff_id)
            #    scaff_id += 1
#print scaffolds_new.keys()
pickle.dump(scaffolds_new,open(args.directory+'/scaffolds_iteration_'+str(int(args.iteration) -1)+'.p','w'))
update_bed(scaffolds_new)
ofile = open(args.directory+'/misasm_'+args.iteration+'.DONE','w')
ofile.close()
