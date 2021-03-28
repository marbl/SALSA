#include <iostream>
#include <algorithm>
#include <map>
#include <string>
#include <cstring>
#include <fstream>
#include <cmath>
#include <getopt.h>
#include <vector>
#include <unordered_map>

#include <boost/config.hpp>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/undirected_graph.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

using namespace std;

typedef float Weight;
typedef boost::property<boost::edge_weight_t, Weight> WeightProperty;
typedef boost::property<boost::vertex_name_t, std::string> NameProperty;

typedef boost::adjacency_list < boost::listS, boost::vecS, boost::directedS,NameProperty, WeightProperty > Graph;

typedef boost::graph_traits < Graph >::vertex_descriptor Vertex;

typedef boost::property_map < Graph, boost::vertex_index_t >::type IndexMap;
typedef boost::property_map < Graph, boost::vertex_name_t >::type NameMap;

typedef boost::iterator_property_map < Vertex*, IndexMap, Vertex, Vertex& > PredecessorMap;
typedef boost::iterator_property_map < Weight*, IndexMap, Weight, Weight& > DistanceMap;

// typedef boost::property<boost::edge_weight_t, double> EdgeWeightProperty;
// typedef boost::adjacency_list<boost::listS, boost::vecS,boost::bidirectionalS,boost::no_property,EdgeWeightProperty> Graph;

// typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> Graph;

std::vector<std::string> split(std::string const & s, size_t count)
{
  vector<string> ret;

  for(unsigned int i = 0;i < s.length();i+=count)
  {
    ret.push_back(s.substr(i,count));
  }
    return ret;
}


char* getCharExpr(string s)
{
    char *a=new char[s.size()+1];
    a[s.size()]=0;
    memcpy(a,s.c_str(),s.size());
    return a;
}

void print_usage()
{
	printf("Usage: misasm -g gfa -l hic_links\n");
}

int main(int argc, char* argv[])
{
	int option = 0;
	string gfa="", hic_links="";
	while((option = getopt(argc,argv,"g:l:")) != -1)
	{
		switch(option)
		{
			case 'l':
			{
				//cout<<optarg<<endl;
				hic_links = optarg;
				break;
			}
			case 'g':
			{
				//cout<<optarg<<endl;
				gfa = optarg;
				break;
			}
			default:
			{
				//cout<<"here"<<endl;
				print_usage();
				exit(EXIT_FAILURE);
			}
		}
	}

	if(gfa == "" || hic_links == "")
	{
		print_usage();
		exit(EXIT_FAILURE);
	}

	unordered_map<string,int> unitig2id;
	unordered_map<int,string> id2unitig;
	int id = 0;
	Graph g;
	string line;
	ifstream edges(getCharExpr(gfa));
	while(getline(edges,line))
	{
		string u,v,y,z,a,b;
		double w = 1.0;
		istringstream iss(line);
		/*
		S       tig00000042     *       LN:i:30165
		*/
		if(line[0] == 'S')
		{
			iss >> u >> v >> y >> z;
			if(unitig2id.find(v+":FOW") == unitig2id.end())
			{
				unitig2id[v+":FOW"] = id;
				id2unitig[id] = v+":FOW";
				id++;
			}
			if(unitig2id.find(v+":REV") == unitig2id.end())
			{
				unitig2id[v+":REV"] = id;
				id2unitig[id] = v+":REV";
				id++;
			}
		}
		/*
		L       tig00000014     -       tig00000001     -       5I1M5I
		*/
		if(line[0] == 'L')
		{
			iss >> u >> v >> y >> z >> a >> b;
			string contig1 = v;
			string contig2 = z;
			string or1 = y;
			string or2 = a;
			string v1,v2,v3,v4;
			//cout<<or1<<"\t"<<or2<<endl;
			if(or1 == "+" && or2 == "+")
			{
				v1 = contig1 + ":FOW";
				v2 = contig2 + ":FOW";
				v3 = contig2 + ":REV";
				v4 = contig1 + ":REV";
			}
			if(or1 == "+" && or2 == "-")
			{
				v1 = contig1 + ":FOW";
				v2 = contig2 + ":REV";
				v3 = contig2 + ":FOW";
				v4 = contig1 + ":REV";
			}
			if(or1 == "-" && or2 == "+")
			{
				v1 = contig1 + ":REV";
				v2 = contig2 + ":FOW";
				v3 = contig2 + ":REV";
				v4 = contig1 + ":FOW";
			}
			if(or1 == "-" && or2 == "-")
			{
				v1 = contig1 + ":REV";
				v2 = contig2 + ":REV";
				v3 = contig2 + ":FOW";
				v4 = contig1 + ":FOW";
			}
			boost::add_edge(unitig2id[v1],unitig2id[v2],w,g);
			boost::add_edge(unitig2id[v3],unitig2id[v4],w,g);
		}
	}
	edges.close();
	cerr<<"Done loading GFA file"<<endl;
	cerr<<"Number of Nodes = "<<boost::num_vertices(g)<<endl;
	cerr<<"Number of Edges = "<<boost::num_edges(g)<<endl;


	/*
	Assume that HiC links are given with in the decreasing order of weights.
	Also, process only those shortest path queries for which no link is added before.
	This would reduce tons of computation time
	*/

    unordered_map<string,bool> added;
	ifstream hic(getCharExpr(hic_links));
	while(getline(hic,line))
	{
		string a,b,x,h;
		double c,d,e,f;
		istringstream iss(line);
		iss >> a >> b >> c >> d >> e >> f >> x >> h;
		if(e < 1)
		{
			cout<<line<<endl;
            continue;
		}
		string u = a.substr(0,a.length()-2);
		string v = b.substr(0,b.length()-2);
        if(unitig2id.find(u+":FOW") == unitig2id.end() && unitig2id.find(u+":REV") == unitig2id.end() && unitig2id.find(v+":FOW") == unitig2id.end() && unitig2id.find(v+":REV") == unitig2id.end())
        {
            cout<<line<<endl;
            continue;
        }
		// Create things for Dijkstra
		std::vector<Vertex> predecessorsFOW(boost::num_vertices(g)); // To store parents
		std::vector<Weight> distancesFOW(boost::num_vertices(g)); // To store distances

		IndexMap indexMap = boost::get(boost::vertex_index, g);
		PredecessorMap predecessorMapFOW(&predecessorsFOW[0], indexMap);
		DistanceMap distanceMapFOW(&distancesFOW[0], indexMap);

		boost::dijkstra_shortest_paths(g, unitig2id[u + ":FOW"], boost::distance_map(distanceMapFOW).predecessor_map(predecessorMapFOW));

		std::vector<Vertex> predecessorsREV(boost::num_vertices(g)); // To store parents
		std::vector<Weight> distancesREV(boost::num_vertices(g)); // To store distances

		PredecessorMap predecessorMapREV(&predecessorsREV[0], indexMap);
		DistanceMap distanceMapREV(&distancesREV[0], indexMap);

		boost::dijkstra_shortest_paths(g, unitig2id[u + ":REV"], boost::distance_map(distanceMapREV).predecessor_map(predecessorMapREV));

		long minpath = -1;
		string ori;
		double ratio = 0.1;
		vector<string> orientations;
		orientations.push_back("BB");
		orientations.push_back("BE");
		orientations.push_back("EB");
		orientations.push_back("EE");
		bool done = false;
		bool printed = false;
    bool inlimits = false;
    for(unsigned int i = 0;i < orientations.size();i++)
		{
            printed = false;
			string orientation = orientations[i];
			string v1,v2;
			if(orientation == "BB")
			{
				v1 = u + ":REV";
				v2 = v + ":FOW";
			}
			if(orientation == "BE")
			{
				v1 = u + ":REV";
				v2 = v + ":REV";
			}
			if(orientation == "EB")
			{
				v1 = u + ":FOW";
				v2 = v + ":FOW";
			}
			if(orientation == "EE")
			{
				v1 = u + ":FOW";
				v2 = v + ":REV";
			}
			// cout<<boost::edge(unitig2id[v1],unitig2id[v2],g).second<<endl;
			if(added.find(a) == added.end() && added.find(b) == added.end())
			{
				if(boost::edge(unitig2id[v1],unitig2id[v2],g).second)
				{
					//cout<<line<<endl;
					//cout<<"EDGE exists"<<endl;
                    added[u] = true;
					added[v] = true;
					cout<<u<<":"<<orientation[0]<<"\t"<<v<<":"<<orientation[1]<<"\t"<<c<<"\t"<<d<<"\t"<<e<<"\t"<<f<<"\t"<<x<<"\t"<<h<<endl;
					done = true;
					printed = true;
                    break;
				}
				else
				{
					//cout<<"here"<<endl;
					//string linkorient = string(a[a.length()-1]) + string(b[b.length()-1]);
					//extract contig name
				    //cout<<"here"<<endl;
					long path_len;
					if(orientation == "BB" || orientation == "BE")
					{
						path_len = distanceMapREV[unitig2id[v2]];
					    //cout<<"PATH LENGTH = "<<path_len<<endl;
                    }
					else
					{
						path_len = distanceMapFOW[unitig2id[v2]];
                        //cout<<"PATH LENGTH = "<<path_len<<endl;
					}
					if(path_len < INT_MAX && path_len > INT_MIN)
					{
                        //cout<<"Path Length = "<<path_len<<endl;
                        inlimits=true;
						if(minpath == -1)
						{
                            //cout<<"setting"<<endl;
                            //cout<<line<<endl;
                            //cout<<orientation<<endl;
							minpath = path_len;
							ori = orientation;
						}
						else
						{
							if(path_len < minpath)
							{
								ratio = 1.0*path_len/minpath;
								minpath = path_len;
								ori = orientation;
							}
						}
					}
				}
			}
		}
		//if(!done &&  minpath < INT_MAX && minpath > INT_MIN)
		//{
		// 	cout<<ratio<<"\t"<<minpath<<endl;
		//}
        if(!inlimits)
        {
            continue;
        }
		if(ratio > 0 && ratio < 1 && minpath < INT_MAX && minpath > INT_MIN && !done)
		{
			//cout<<ratio<<"\t"<<minpath<<endl;
			//cout<<"Where we found the paths in graph"<<endl;
            cout<<u+":"+ori[0]<<'\t'<<v+":"+ori[1]<<'\t'<<c<<"\t"<<d<<"\t"<<e<<"\t"<<f<<"\t"<<x<<"\t"<<h<<endl;
		    //cout << line << endl;
			added[a] = true;
            printed = true;
			added[b] = true;
		}
		else
		{
            //cout<<"here"<<endl;
			if(!printed)
            {
                cout<<line<<endl;
			    added[a] = true;
			    added[b] = true;
            }
		}
        //cout<<"======================="<<endl;
	}

	// BGL_FORALL_VERTICES(v, g, Graph)
 //  	{
 //  		cout<<id2unitig[1]<<"\t"<<id2unitig[v]<<"\t"<<distanceMap[v]<<endl;
 //  	}
	return 0;
}
