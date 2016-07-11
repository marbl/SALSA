import networkx as nx
import operator
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import argparse

revcompl = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A'}[B] for B in x][::-1])

H = nx.Graph()
G = nx.Graph()
oriented = nx.Graph()

existing_nodes = set()
contigs = set()
edgemap = {}
contig_length = {}
nodes_to_edges = {}

id2seq = {}

recs = []
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c","--cleaned",help="cleaned assembly")
    args = parser.parse_args()

    input_handle = open(args.cleaned, "rU") 

    for record in SeqIO.parse(input_handle, "fasta"):
        id,seq = record.id, str(record.seq)
        #print id
        id2seq[id] = seq

    def break_cycle(nodes):
        nodeset = set()
        for node in nodes:
            nodeset.add(node.split(":")[0])
        nodeset = list(nodeset)
        weight = ""
        chosen_edge = ""
        if len(nodeset) == 2:
            u = nodeset[0]
            v = nodeset[1]
            if u in nodes_to_edges:
                edges = nodes_to_edges[u]
                sorted_edges = sorted(edges,key=operator.itemgetter(2),reverse=True)
                #print sorted_edges
                for each in sorted_edges:
                    if each[1].split(":")[0] == v:
                        weight = each[2]
                        chosen_edge = each

            if v in nodes_to_edges:
                edges = nodes_to_edges[v]
                sorted_edges = sorted(edges,key=operator.itemgetter(2),reverse=True)
                for each in sorted_edges:
                    if each[1].split(":")[0] == u:
                        if each[2] > weight:
                            weight = each[2]
                            chosen_edge = each
                        if each[2] < weight:
                            break

            if chosen_edge != "":
                start = chosen_edge[0]
                end = chosen_edge[1]
                path = []
                if start.split(":")[1] == 'B':
                    path.append(start.split(":")[0]+":E")
                    #path.append(start.split(":")[0]+":B")
                else:
                    path.append(start.split(":")[0]+":B")
                    #path.append(start.split(":")[0]+":E")

                path.append(start)
                path.append(end)

                if end.split(":")[1] == 'B':
                    path.append(end.split(":")[0]+":E")
                    #path.append(end.split(":")[0]+":B")
                else:
                    path.append(end.split(":")[0]+":B")
                    #path.append(end.split(":")[0]+":E")

                return path

    with open('contig_lengths_new','r') as f:
        lines = f.readlines()
        for line in lines:
            attrs = line.split()
            contig_length[attrs[0]] = long(attrs[-1])


    with open("new_links_sorted",'r') as f:
        for row in f:
            row = row.strip().split()
            v1, v2 = row[0:2]
            score = float(row[-1])
            count = float(row[4])
            c1 = v1.split(":")[0]
            c2 = v2.split(":")[0]
            contigs.add(c1)
            contigs.add(c2)

            if c1 not in nodes_to_edges:
                nodes_to_edges[c1] = []
            if c2 not in nodes_to_edges:
                nodes_to_edges[c2] = []
            nodes_to_edges[c1].append((v1,v2,float(row[4])))
            
            key =  c1+'$'+c2
            H.add_edge(c1,c2,weight=int(row[-1]))
            #oriented.add_edge(v1,v2,weight=int(row[-1]),count=float(row[3]))
            oriented.add_edge(v1,v2,weight=int(row[-1]))
            #print key
            if key not in edgemap:
                edgemap[key] =  int(row[-1])
            else:
                edgemap[key] +=  int(row[-1])

            key =  c2+'$'+c1    
            if key not in edgemap:
                edgemap[key] =  int(row[-1])
            else:
                edgemap[key] +=  int(row[-1])

            if v1 not in existing_nodes and v2 not in existing_nodes:
                if count < 25:
                    continue
                G.add_edge( v1, v2, score=score, t="x" )
                existing_nodes.add(v1)
                existing_nodes.add(v2)

    for ctg in list(contigs):
        G.add_edge( ctg+":B", ctg+":E", t="c" )

    g_idx = 1
    recs = []
    to_merge = set()

    backbone_paths = {}
    path_id = 1
    for subg in nx.connected_component_subgraphs(G):
        p0 = []

        for v in subg.nodes():
            if subg.degree(v) == 1:
                p0.append(v)

        if len(p0) != 2:
            path = break_cycle(subg.nodes())
            if path != None:
                #print path
                if len(path) == 2:
                    to_merge.add(path[0].split(':')[0])
                    continue
                backbone_paths[path_id] = path
                path_id += 1
            else:
                for each in subg.nodes():
                    to_merge.add(each.split(':')[0])
            continue
        else:
            path = nx.shortest_path(subg, p0[0], p0[1])
            if len(path) == 2:
                to_merge.add(path[0].split(':')[0])
                continue
            backbone_paths[path_id] = path
            path_id += 1
            curr_contig = ""

        g_idx += 1

    #now for each separate contig, find a maximum likely backbone path
    assignment = {}
    for each in to_merge:
        max_sum = -1
        max_path = -1
        for key in backbone_paths:
            path = backbone_paths[key]
            cur_sum = 0
            cnt = 0
            for node in path:
                if H.has_edge(each,node.split(':')[0]):
                    cur_sum += H[each][node.split(':')[0]]['weight']
                    cnt += 1
            if cnt != 0 and cur_sum > max_sum:
                max_sum = cur_sum
                max_path = key

        if max_sum != -1:
            assignment[each] = (max_path, max_sum,contig_length[each])
                

    #now that we have found the path, try putting contig at best position in the path

    count = len(assignment)

    path_to_contig = {}

    for each in assignment:
        key = assignment[each][0]
        if key not in path_to_contig:
            path_to_contig[key] = []
        path_to_contig[key].append((each, assignment[each][1],assignment[each][2]))

    for each in path_to_contig:
        contigs = path_to_contig[each]
        contigs_sorted = sorted(contigs, key = operator.itemgetter(1), reverse = True)
        path_to_contig[each] = contigs_sorted
        #print contigs_sorted


    ofile = open('ambigous_contigs','w')

    for path_id in path_to_contig:
        path = backbone_paths[path_id]
        temp_path = list(path)
        
        contigs = path_to_contig[path_id]
        contigs = [str(i[0]) for i in contigs]
        for contig in contigs:
            begin = contig + ":B"
            end = contig + ":E"
            total_max = -1
            orientation = ''
            pos = -1
            #check for positions in the middle of the path
            for i in range(1,len(path)-1,2):
                score_fow = -1
                score_rev = -1

                if oriented.has_edge(path[i],begin) and oriented.has_edge(end,path[i+1]):
                    score_fow = oriented[path[i]][begin]['weight'] + oriented[end][path[i+1]]['weight']

                if oriented.has_edge(path[i],end) and oriented.has_edge(begin,path[i+1]):
                    score_rev = oriented[path[i]][end]['weight'] + oriented[begin][path[i+1]]['weight']

                if score_fow >= score_rev:
                    if score_fow > total_max:
                        total_max = score_fow
                        orientation = 'fow'
                        pos = i
                else:
                    if score_rev > total_max:
                        total_max = score_rev
                        orientation = 'rev'
                        pos = i

                #check for start and end
                if oriented.has_edge(begin,path[0]):
                    score_fow = oriented[begin][path[0]]['weight']
                    if score_fow > total_max:
                        total_max = score_fow
                        orientation = 'fow'
                        pos = 0

                if oriented.has_edge(end,path[0]):
                    score_rev = oriented[end][path[0]]['weight']
                    if score_rev > total_max:
                        total_max = score_rev
                        orientation = 'rev'
                        pos = 0

                if oriented.has_edge(path[-1],begin):
                    score_fow = oriented[path[-1]][begin]['weight']
                    if score_fow > total_max:
                        total_max = score_fow
                        orientation = 'fow'
                        pos = len(path)

                if oriented.has_edge(path[-1],end):
                    score_rev = oriented[path[-1]][end]['weight']
                    if score_rev > total_max:
                        total_max = score_rev
                        orientation = 'rev'
                        pos = len(path)

            if total_max == -1:
                for i in range(2,len(path)-1,2):
                    if oriented.has_edge(end,path[i]):
                        if oriented[end][path[i]]['weight'] > total_max:
                            total_max = oriented[end][path[i]]['weight']
                            orientation = 'fow'
                            pos = i-1
                    if oriented.has_edge(begin,path[i]):
                        if oriented[begin][path[i]]['weight'] > total_max:
                            total_max = oriented[begin][path[i]]['weight']
                            orientation = 'rev'
                            pos = i-1

                    if oriented.has_edge(path[i+1],begin):
                        if oriented[path[i+1]][begin]['weight'] > total_max:
                            total_max = oriented[path[i+1]][begin]['weight']
                            orientation = 'fow'
                            pos = i+1

                    if oriented.has_edge(path[i+1],end):
                        if oriented[path[i+1]][end]['weight'] > total_max:
                            total_max = oriented[path[i+1]][end]['weight']
                            orientation = 'rev'
                            pos = i+1


            if total_max > 25:
                if orientation == 'fow':
                    if pos == 0:
                        path.insert(0,end)
                        path.insert(0,begin)
                    else:
                        if pos == len(path):
                            path.append(begin)
                            path.append(end)
                        else:
                            path.insert(pos+1,begin)
                            path.insert(pos+2,end)

                else:
                    if pos == 0:
                        path.insert(0,begin)
                        path.insert(0,end)
                    else:
                        if pos == len(path):
                            path.append(end)
                            path.append(begin)
                        else:
                            path.insert(pos+1,end)
                            path.insert(pos+2,begin)
            # else:
            #     '''
            #     The format is as follows
            #     small_contig    contig_to_the_left contig_to_the_right 
            #     if contig is to be place to the beginning or end, just one contig
            #     '''



        backbone_paths[path_id] = path

    for key in backbone_paths:
        if len(backbone_paths[key]) >= 4:
            print backbone_paths[key]

    # for key1 in backbone_paths:
    #     max_weight = 0
    #     max_path = ''
    #     for key2 in backbone_paths:
    #         if key1 != key2:
    #             path1 = backbone_paths[key1]
    #             path2 = backbone_paths[key2]
    #             weight = 0
    #             for contig1 in path1:
    #                 ctg1 = contig1.split(':')[0]
    #                 for contig2 in path2:
    #                     ctg2 = contig2.split(':')[0]
    #                     if H.has_edge(ctg1,ctg2):
    #                         weight += H[ctg1][ctg2]['weight']
    #             if weight > max_weight:
    #                 max_weight = weight
    #                 max_path = key2

    #     if max_path != '' and 1000 < max_weight < 4000:
    #         print backbone_paths[key1], backbone_paths[max_path], max_weight

    c_id = 1
    for key in backbone_paths:
        if len(backbone_paths[key]) >= 4:
            path = backbone_paths[key]
            curr_contig = ""
            for i in range(0,len(path)-1,2):
                curr = path[i]
                next = path[i+1]
                curr = curr.split(':')
                next = next.split(':')
                if curr[1] == 'B' and next[1] == 'E':
                    curr_contig += id2seq[curr[0]]
                if curr[1] == 'E' and next[1] == 'B':
                    curr_contig += revcompl(id2seq[curr[0]])
                if i != len(path) - 2:
                    for j in range(0,500):
                        curr_contig += 'N'
            rec = SeqRecord(Seq(curr_contig,generic_dna),id='scaffold_'+str(c_id))
            recs.append(rec)
            c_id += 1
    SeqIO.write(recs,'scaffolds.fasta','fasta')

if __name__ == '__main__':
    main()
