import networkx as nx
import operator
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import argparse

revcompl = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A'}[B] for B in x][::-1])



id2seq = {}

nodes_to_edges = {}

def main():


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

            #print chosen_edge
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

    parser = argparse.ArgumentParser()
    parser.add_argument("-c","--clean",help="cleaned assembly") 
    args = parser.parse_args()
    input_handle = open(args.clean, "rU") 

    for record in SeqIO.parse(input_handle, "fasta"):
        id,seq = record.id, str(record.seq)
        #print id
        id2seq[id] = seq


    G = nx.Graph()

    existing_nodes = set()
    contigs = set()
    with open("new_links_sorted") as f:
        for row in f:
            row = row.strip().split()
            v1, v2 = row[0:2]
            score = float(row[2])
            #count = float(row[3])
            c1 = v1.split(":")[0]
            c2 = v2.split(":")[0]
            contigs.add(c1)
            contigs.add(c2)
            if c1 not in nodes_to_edges:
                nodes_to_edges[c1] = []
            if c2 not in nodes_to_edges:
                nodes_to_edges[c2] = []
            nodes_to_edges[c1].append((v1,v2,float(row[4])))

            if v1 not in existing_nodes and v2 not in existing_nodes:
                #if score < 15 and count < 800:
                #    continue
                G.add_edge( v1, v2, score=score, t="x" )
                existing_nodes.add(v1)
                existing_nodes.add(v2)

    for ctg in list(contigs):
        G.add_edge( ctg+":B", ctg+":E", t="c" )

    g_idx = 1
    recs = []
    for subg in nx.connected_component_subgraphs(G):
        
        p0 = []
        curr_contig = ""
        for v in subg.nodes():
            if subg.degree(v) == 1:
                p0.append(v)
        if len(p0) != 2:
            path = break_cycle(subg.nodes())
            if path != None:
                print path
                #curr_contig = ""
                for i in range(0,len(path)-1,2):
                    curr = path[i]
                    next = path[i+1]
                    curr = curr.split(':')
                    next = next.split(':')
                    if curr[1] == 'B' and next[1] == 'E':
                        print 'adding ' + curr[0]
                        curr_contig += id2seq[curr[0]]
                    if curr[1] == 'E' and next[1] == 'B':
                        print 'adding rev ' + curr[0]
                        curr_contig += revcompl(id2seq[curr[0]])
                rec = SeqRecord(Seq(curr_contig,generic_dna),id=str(g_idx))
                recs.append(rec)
                print g_idx,subg.nodes()
            
        else:
            path = nx.shortest_path(subg, p0[0], p0[1])
            for each in path:
                print g_idx, each
            for i in range(0,len(path)-1,2):
                curr = path[i]
                next = path[i+1]
                curr = curr.split(':')
                next = next.split(':')
                if curr[1] == 'B' and next[1] == 'E':
                    curr_contig += id2seq[curr[0]]
                if curr[1] == 'E' and next[1] == 'B':
                    curr_contig += revcompl(id2seq[curr[0]])
            rec = SeqRecord(Seq(curr_contig,generic_dna),id=str(g_idx))
            recs.append(rec)
        
        g_idx += 1


    SeqIO.write(recs,'scaffolds.fasta','fasta')

        



    nx.write_gexf(G, "scaffolding.gexf")

if __name__ == '__main__':
    main()

