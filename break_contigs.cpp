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

int main()
{
	ifstream bedfile("alignment_unique.bed");
	string line;
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
	ofstream outputfile("alignment_unique_clean.bed");
	ofstream outputfile2("breakpoints");
	//ofstream outputfile1("length_to_coverage");
	/*
	for(unordered_map<string,vector<pair<long,long> > >  :: iterator it = intervals.begin(); it != intervals.end(); ++it)
	{
		vector<pair<long,long> > longs = it->second;
		sort(longs.begin(),longs.end());
		intervals[it->first] = longs;
	}*/
	ifstream lenfile("contig_lengths");
	while(getline(lenfile,line))
	{
		istringstream iss(line);
		string a;
		long n;
		iss >> a >> n;
		contig_length[a] = n;
	}
	
	//const char* args[] = {"000001F"};
	//const char* args[] = {"chr6"};
	vector<string> contigs;
	for(unordered_map<string,vector<pair<long,long> > > :: iterator it = intervals.begin(); it != intervals.end(); ++it)
	{
		contigs.push_back(it->first);
	}

	for(int i = 0;i < contigs.size();i++)
	{
		string curr_contig = contigs[i];
		//cerr<<curr_contig<<endl;
		long length = contig_length[curr_contig];
		vector<pair<long,long> > segments;
		if(contig_length[curr_contig] > 10000000)
		{
			//cout<<"Here"<<endl;
			long from  = 0;
			//ofstream ofile(getCharExpr("results/"+curr_contig+".coords.full"));
			//sort(intervals[curr_contig].begin(), intervals[curr_contig].end());
			vector<pair<long,long> > interval = intervals[curr_contig];
			vector<int> counts(length+1,0);


			//now load the intervals and store start and end points in priority queue
			// priority_queue<long, vector<long>, greater<long> > start_points;
			// priority_queue<long, vector<long>, greater<long> > end_points;

			long count = 0;

			const long sz = counts.size();
			
			long len = interval.size();
			for(long i = 0;i < len;i++)
			{
				pair<long,long> curr_interval = interval.at(i);
				// start_points.push(curr_interval.first);
				// end_points.push(curr_interval.second);
				counts.at(curr_interval.first) += 1;
				counts.at(curr_interval.second + 1) -= 1;
			}

			contig_coverage[curr_contig] = 0;
			for(long i = 1; i < sz; i++)
			{	
				//cout<<i<<endl;
				counts.at(i) += counts.at(i-1);
				contig_coverage[curr_contig] += counts.at(i);
			}

			//new idea
			//cout<<"before"<<endl;
			double threshold = contig_coverage[curr_contig] / 5 / contig_length[curr_contig];
			//cout<<"Threshold = "<<threshold<<endl;
			vector<int> threshold_array(sz-4000000);
			//cout<<"here"<<endl;
			for(long i = 0; i < sz-4000000;i++)
			{
				if(counts[i+2000000] <= threshold)
					threshold_array[i] = 1;
				else
					threshold_array[i] = -1;
			}
			//cout<<"arrray created"<<endl;
			long max_so_far = 0;
			long max_ending_here = 0;
			long start_index = 0;
			long max_start_index = 0;
			long max_end_index = 0;
			for (long i = 0; i < sz-4000000; i++)
		    {
		        max_ending_here = max_ending_here + threshold_array[i];
		        if (max_ending_here < 0)
		        {
		            max_ending_here = 0;
		            start_index = i + 1;
		        }
		        if (max_so_far < max_ending_here)
		        {
		            max_so_far = max_ending_here;
		            max_start_index = start_index;
		            max_end_index = i;
		        }
		    }
		    //cerr<<"sum = " << max_so_far<<endl;
		    long actual_start_index = max_start_index + 2000000;
			long actual_end_index = max_end_index + 2000000;
		    if(sz - actual_end_index > 2000000 && max_start_index != 0)
		    {
				
			//	cout<<actual_start_index<<"\t"<<actual_end_index<<endl;
				segments.push_back(make_pair(0,actual_start_index));
				segments.push_back(make_pair(actual_start_index+1,sz));
				
			}
			else
			{
				segments.push_back(make_pair(0,sz));
			}
			
		}
		else
		{
			segments.push_back(make_pair(0,length+1));
		}

		vector<BedRecord> curr_record = contig_to_record[curr_contig];
		for(int i = 0; i < curr_record.size();i++)
		{
			BedRecord rec = curr_record[i];
			long start = rec.start;
			long end = rec.end;
			for(int j = 0;j < segments.size();j++)
			{
				pair<long,long> seg = segments[j];
				if(start >= seg.first && end <= seg.second)
				{
					string new_id = curr_contig +"_" + to_string(seg.first) +"_" + to_string(seg.second);
					long new_start = start - seg.first ;
					long new_end =  end - seg.first;
					outputfile << new_id << "\t" << new_start << "\t" << new_end << "\t" << rec.read << "\t" << rec.flag << "\t" << rec.strand << endl;
				}
			}
		}
		outputfile2<<curr_contig<<endl;
		for(long i = 0;i < segments.size();i++)
		{
			outputfile2<<segments.at(i).first<<"\t"<<segments.at(i).second<<endl;
		}
	}
}

