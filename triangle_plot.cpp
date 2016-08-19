#include <iostream>
#include <unordered_map>
#include <vector>
#include <fstream>
#include <string>
#include <cstring>
#include <algorithm>
#include <sstream>
#include <cmath>
#include <queue>
#include <set>

#include "cmdline.h"

using namespace std;

class BedRecord
{
public:
        string contig;
        long start;
        long end;
        string read;
        int flag;
        char strand;//+ forward - reverse
        BedRecord () {}
        BedRecord(string contig, long start, long end, string read, int flag, char strand);

};

BedRecord :: BedRecord(string contig, long start, long end, string read, int flag, char strand)
{
        this->contig = contig;
        this->start = start;
        this->end = end;
        this->read = read;
        this->flag = flag;
        this->strand = strand;
}

unordered_map<string,BedRecord> first_in_pair;
unordered_map<string, BedRecord> second_in_pair;
unordered_map<string,vector<pair<long,long> > > intervals;
unordered_map<string,long> contig_length;
unordered_map<string, vector<BedRecord> > contig_to_record;
unordered_map<string, double> contig_coverage;

pair<long,long> check(long from, long polong, vector<pair<long,long>> intervals)
{
	long ret = 0;
	long i;
	for(i = from; i < intervals.size();i++)
	{
		pair<long,long> curr = intervals[i];
		if(polong < curr.first)
			return make_pair(ret,i);
		else
		{
			if(polong >= curr.first && polong <= curr.second)
				ret++;
		}
	}
	return make_pair(ret,i);
}

char* getCharExpr(string s)
{
        char *a=new char[s.size()+1];
        a[s.size()]=0;
        memcpy(a,s.c_str(),s.size());
        return a;
}

int main(int argc, char *argv[])
{
	cmdline::parser p;
	p.add<string>("alignment", 'a', "bed file for alignment", true, "");
	p.add<string>("output", 'o', "coordinate output file", true, "");
	p.parse_check(argc, argv);	
	set<string> contigs;

	ifstream bedfile(getCharExpr(p.get<string>("alignment")));
	string line;
	string reference;
	while(getline(bedfile,line))
	{
		string contig, read;
		char strand;
		long start,end,flag;
		istringstream iss(line);
		if(!(iss >> contig >> start >> end >> read >> flag >> strand))
			break;
		// if(contig == "000008F")
		// {
			contigs.insert(contig);
			BedRecord rec(contig,start,end,read,flag,strand);
			if(contig_to_record.find(contig) == contig_to_record.end())
			{
				vector<BedRecord> recs;
				contig_to_record[contig] =recs;
			}
			else
			{
				contig_to_record[contig].push_back(rec);
			}
			if(read[read.length() -1 ] == '1')
	        {
	                first_in_pair[read.substr(0,read.length()-2)] = rec;
	        }
	        else
	        {
	                second_in_pair[read.substr(0,read.length()-2)] = rec;
	    	}
 		// }
		
	}

	for(unordered_map<string,BedRecord> :: iterator it = first_in_pair.begin(); it != first_in_pair.end(); ++it)
	{
		BedRecord rec1 = it->second;
		if(second_in_pair.find(it->first) != second_in_pair.end())
		{
			BedRecord rec2 = second_in_pair[it->first];
			if(rec1.contig == rec2.contig)
			{
				if(intervals.find(rec1.contig) == intervals.end())
				{
					vector<pair<long,long> > val;
					intervals[rec1.contig] = val;
				}
				long start,end;
				if(rec1.start <= rec2.start)
				{
					start = rec1.start;
					end = rec2.end;
				}
				else	
				{
					start = rec2.start;
					end = rec1.end;
				}
				//if(end - start <= 500000)
				intervals[rec1.contig].push_back(make_pair(start,end));
			}
		}
	}

	bedfile.close();
	cerr<<"bedfile loaded"<<endl;
	ofstream outputfile(getCharExpr(p.get<string>("output")));

	for(unordered_map<string,vector<pair<long,long> > > :: iterator it = intervals.begin(); it != intervals.end(); ++it)
	{
		string curr_contig = it->first;
		cerr<<curr_contig<<endl;
		outputfile<<curr_contig<<endl;
		vector<pair<long,long> > interval = it->second;
		long len = interval.size();
		for(long i = 0;i < len;i++)
		{
			pair<long,long> curr_interval = interval.at(i);
			double x = (curr_interval.first + curr_interval.second)/2.0;
			double y = abs(curr_interval.first - curr_interval.second)/2.0;
			outputfile<<x<<"\t"<<y<<endl;
			
		}
	}
	return 0;
}

