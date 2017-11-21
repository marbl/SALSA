import networkx as nx
import sys
import argparse

G = nx.Graph()
parser = argparse.ArgumentParser()
parser.add_argument("-i","--iteration", help="iteration number")
parser.add_argument("-d",'--directory',help="output directory")
args = parser.parse_args()

iteration = int(args.iteration)
def get_max_incident(start, end):
   max=0
   try:
      for e in G.successors(start):
         if (e != end and G[start][e]['weight'] > max):
            max = G[start][e]['weight']
      for e in G.predecessors(start):
         if (e != end and G[e][start]['weight'] > max):
            max = G[start][e]['weight']
   except:
      for e in G.neighbors(start):
         if (e != end and G[start][e]['weight'] > max):
            max = G[start][e]['weight']

   return max

#map for node to max edge weight for it
node2weight = {}

def get_max_weight(start,end):
   return max(get_max_incident(start, end), get_max_incident(end, start))

f = open(args.directory+'/contig_links_iteration_'+str(args.iteration),'r')
for line in f:
   line = line.strip().split()

   #LNK tig00001320_pilon tig00019085_pilon E:B 103.138095238 21659 True True True
   #LNK tig00006867_pilon tig00009112_pilon E:E 0.214285714286 1 True True True
#   if line[0] != "LNK":
#     continue

#   edgeType=line[3].split(":")
   G.add_node(line[0])
   G.add_node(line[1])
   if line[0] not in node2weight:
      node2weight[line[0]] = float(line[2])
   else:
      if float(line[2]) >= node2weight[line[0]]:
         node2weight[line[0]] = float(line[2])

   if line[1] not in node2weight:
      node2weight[line[1]] = float(line[2])
   else:
      if float(line[2]) >= node2weight[line[1]]:
         node2weight[line[1]] = float(line[2])

   G.add_edge(line[0], line[1], weight=float(line[2]),isNeighbor="False", isGood="False",links=int(line[3]))
f.close()

ofile = open(args.directory+'/contig_links_scaled_iteration_'+str(args.iteration),'w')
#print 'Graph Loaded Done'
f1 = open(args.directory+'/contig_links_iteration_'+str(args.iteration),'r')
for line in f1:
   #bestAlt=get_max_weight(u, v)
   attrs = line.strip().split()
   #print 'here'
   u = attrs[0]
   v = attrs[1]
   w = float(attrs[2])
   max_u = 0
   max_v = 0
   if u in node2weight:
      max_u = node2weight[u]
   if v in node2weight:
      max_v = node2weight[v]

   bestAlt = max(max_u,max_v)
   if bestAlt == w:
      #print 'here'
      bestAlt=get_max_weight(u, v)
   if bestAlt==0:
      bestAlt=1

   ofile.write(str(u)+'\t'+str(v)+'\t'+str(w)+'\t'+str(bestAlt)+'\t'+str(w/bestAlt)+'\t'+str(G[u][v]['links'])+'\t'+str(G[u][v]['isNeighbor'])+'\t'+str(G[u][v]['isGood'])+'\n')
   #print str(u)+'\t'+str(v)+'\t'+str(w)+'\t'+str(bestAlt)+'\t'+str(w/bestAlt)+'\t'+str("FALSE")+'\t'+str("FALSE")+'\t'+str("FALSE")
   #print "%s %s %f %f %f %d %s %s"%(u, v, d, bestAlt, d/bestAlt, G[u][v]['links'],G[u][v]['isNeighbor'], G[u][v]['isGood'])
ofile.close()
f1.close()