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
#include <numeric>
//#include <thread>

#include "cmdline.h"

using namespace std;

//int nthread = std::thread::hardware_concurrency();

pair<long,long> check(long from, long polong, vector<pair<long,long>> intervals)
{
	long ret = 0;
	unsigned long i;
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

vector<int> slice(const vector<int>& v, long start=0, long end=-1) {
	//cout<<"IN SLICE"<<endl;
    long oldlen = v.size();
    long newlen;

    if (end == -1 or end >= oldlen){
        newlen = oldlen-start;
    } else {
        newlen = end-start;
    }

    vector<int> nv(newlen);

    for (long i=0; i<newlen; i++) {
        nv[i] = v[start+i];
    }
    return nv;
}

char* getCharExpr(string s)
{
        char *a=new char[s.size()+1];
        a[s.size()]=0;
        memcpy(a,s.c_str(),s.size());
        return a;
}

unordered_map<string,vector<pair<long,pair<long,long>>>> contig2breakpoints;
unordered_map<string,long> contig_length;
void load_breakpoints(string file)
{
	string line;
	ifstream bfile(getCharExpr(file));
	string contig;
	long pos;
	while(getline(bfile,line))
	{
		vector<pair<long,pair<long,long>>> breaks;
		istringstream iss(line);
		iss >> contig;
		vector<long>bps;
        while(iss >> pos)
		{
            bps.push_back(pos);

			//breaks.push_back(pos);
            if(pos - 2000000 >=0  && pos + 2000000 <= contig_length[contig])
            {
		        breaks.push_back(make_pair(pos,make_pair(pos - 2000000, pos + 2000000)));
            }
            //cout<<pos<<endl;
		}
       // for(int i = 0;i < bps.size();i++)
       // {
       //     if(i == 0 && i+1 < bps.size())
       //     {
       //         breaks.push_back(make_pair(bps[i],make_pair(0,bps[i+1])));
       //     }
       //     if(i > 0 && i + 1 < bps.size())
       //     {
       //         breaks.push_back(make_pair(bps[i],make_pair(bps[i-1],bps[i+1])));
       //     }
       //     if(i == bps.size()-1)
       //     {
       //         breaks.push_back(make_pair(bps[i],make_pair(bps[i-1],contig_length[contig])));
       //     }
       // }
		contig2breakpoints[contig] = breaks;
	}
	bfile.close();
}

double inspect(int pos, vector<int>& count)
{
	int low = 0;
	int total = 0;
	long sz = count.size();
    long start = sz/3;
    long end = 2*sz/3;
    long interval = 50000;
    for(long i = interval; i < sz;i+=interval)
	{
		if(pos - i >= start && pos + i < end)
		{
			total += 1;
			if(count[pos - i ] > count[pos] && count[pos + i] > count[pos])
			{
				low += 1;
			}
		}
		else
		{
			break;
		}
	}
    //cout<<"LOW = "<<low<<"\tTOTAL = "<<total<<endl;
	if(total == 0)
	{
		return 0;
	}
	else
	{
		return low*1.0/total;
	}
}

bool check_if_inside(string contig, long start, long end)
{
    if(contig2breakpoints.find(contig) == contig2breakpoints.end())
    {
        return false;
    }
    vector<pair<long,pair<long,long>>> positions = contig2breakpoints[contig];
    for(unsigned int i = 0;i < positions.size();i++)
    {
        if(start >= positions[i].second.first && end <= positions[i].second.second)
        {
            return true;
        }
    }
    return false;

}

int main(int argc, char *argv[])
{
	cmdline::parser p;
	p.add<string>("alignment", 'a', "bed file for alignment", true, "");
	//p.add<string>("outputdir", 'd', "coordinate output file", true, "");
	p.add<string>("breakpoints", 'b', "breakpoints", true, "");
  p.add<int>("min_size",'s',"Minimum mate pair separation for error findng",true,0);
	p.add<string>("contiglen", 'l', "length of contigs", true, "");
	p.add<int>("iteration",'i',"Iteration number",true,0);
  p.parse_check(argc, argv);
	string line;
	ifstream lenfile(getCharExpr(p.get<string>("contiglen")));
  //unused variable
	// int separation = p.get<int>("min_size");
	unordered_map<string,int> contig2cutoff;
    while(getline(lenfile,line))
	{
		istringstream iss(line);
		string a;
		long n;
		iss >> a >> n;
		contig_length[a] = n;
        contig2cutoff[a] = n/10;
	}
  //unused variable
	// int iteration = p.get<int>("iteration");
	lenfile.close();
	load_breakpoints(p.get<string>("breakpoints"));
	//cout<<"loaded breakpoints"<<endl;
	ifstream bedfile(getCharExpr(p.get<string>("alignment")));

	unordered_map<string,vector<int> > contig2coverage;
    string prev_line = "";
	string prev_contig="";
	long prev_start=-1;
	long prev_end=-1;
	string prev_read="";
	while(getline(bedfile,line))
	{
		string contig, read;
		// unused variable
		// char strand;
		// long flag
		long start,end;
		if(prev_line == "")
		{
			prev_line = line;
			continue;
		}
		istringstream iss(line);
		if(!(iss >> contig >> start >> end >> read))
			break;
		//if(contig != "scaffold_0")
		//{
        //    continue;
        //}
		if(contig_length[contig] < 1000000)
		{
			prev_contig = contig;
			prev_start = start;
			prev_end = end;
			prev_read = read;
			continue;
		}
		if(contig2coverage.find(contig) == contig2coverage.end())
		{
			vector<int> cov(contig_length[contig]+10,0);
			contig2coverage[contig] = cov;
		}
		//vector<int> cov = contig2coverage[contig];
		if(read.substr(0,read.length()-2) == prev_read.substr(0,prev_read.length()-2) && prev_contig == contig)
		{
			// cout<<"here"<<endl;
			if(prev_end <= start)
			{
				//cout<<"here"<<endl;
                if(check_if_inside(contig,prev_start,end+1))
                {
                    if(end - prev_start <= contig2cutoff[contig])
                    {
                        contig2coverage[contig].at(prev_start) += 1;
				        contig2coverage[contig].at(end+1) -= 1;
                    }
                }
				// cout<<contig2coverage[contig].at(prev_start)<<endl;
			}
			//contig2coverage[contig] = cov;
		}
		prev_contig = contig;
		prev_start = start;
		prev_end = end;
		prev_read = read;
	}
	bedfile.close();
	//cerr<<"bedfile loaded"<<endl;

    int total_breakpoints = 0;
    int suspicious_breakpoints = 0;

	for(unordered_map<string,vector<int> > :: iterator it = contig2coverage.begin(); it != contig2coverage.end();++it)
	{
		//cout<<"Testing contig " << it->first<<endl;
        string contig = it->first;
        if (contig2breakpoints.find(it->first) != contig2breakpoints.end())
        {
            cout<<contig;
            vector<pair<long,pair<long,long>>> positions = contig2breakpoints[it->first];

            for(unsigned int i = 0;i < positions.size();i++)
            {
                //cout<<"testing "<<i<<endl;
                long misasm_loc = positions[i].first;
                long start = positions[i].second.first;
                long end = positions[i].second.second;

								long start_pos = start;
                //unused variable
								//long end_pos = end;
                vector<int> cov = contig2coverage[it->first];
                vector<int> local_coverage;
                vector<int>tmpcov(cov);
                for(long j = start;j < end;j++)
                {
                    tmpcov[j] += tmpcov[j-1];
                    local_coverage.push_back(tmpcov[j]);;
                }
               // ofstream covfile(to_string(i)+"_coverage.txt");
                //for(long k = 0;k < local_coverage.size();k+=1000)
                //{
                //   covfile<<k<<"\t"<<local_coverage[k]<<endl;
                //}
                //covfile.close();
                contig2coverage[it->first] = cov;
                double average = accumulate(local_coverage.begin(),local_coverage.end(),0.0)/local_coverage.size();
                //cout<<"average coverage = "<<average<<endl;
                vector<long> misasm_pos;
               for(int div = 5; div <= 15; div++)
               {
                   double cutoff = average/div;
                   vector<int>delta;
                   //cout<<"Cutoff = " << cutoff<<endl;
                   for(unsigned long j = 0; j < local_coverage.size();j++)
                   {
                       if(local_coverage[j] < cutoff)
                       {
                           delta.push_back(5);
                       }
                       else
                       {
                           delta.push_back(-5);
                       }
                   }
                   long sz = delta.size();
                   //cout<<"Running Kadane algorithm"<<endl;
                   /*
                    * Now find maximum sum subarray of delta with Kadane's algorithm
                    */
                   int max_so_far = -100000, max_ending_here = 0;
                   long  start = sz/3, end = 2*sz/3,s = sz/3;
                   for(long j = sz/3 ;j < 2*sz/3;j++)
                   {
                       max_ending_here += delta.at(j);
                       if(max_so_far < max_ending_here)
                       {
                           max_so_far = max_ending_here;
                           start = s;
                           end = j;
                       }
                       if(max_ending_here < 0)
                       {
                           max_ending_here = 0;
                           s = j + 1;
                       }
                   }
                   //cout<<start<<"\t"<<end<<endl;
                   if(misasm_loc >= start + start_pos && misasm_loc <= end+start_pos)
                   {
                       misasm_pos.push_back((start+end)/2);
                   }
                   //cout<<"Possible Misassembly with cutoff "<<cutoff<<" in "<<contig<<" at "<<start<<"\t"<<end<<endl;
               }
               //cout<<"size = "<<"\t"<<misasm_pos.size()<<endl;
               if(misasm_pos.size() >= 8)
               {
                   suspicious_breakpoints += 1;
                   cout<<"\t"<<misasm_loc;
               }
               total_breakpoints += 1;
            }
            cout<<endl;
        }
    }
	cout<<"Total Joins = "<<total_breakpoints<<endl;
    cout<<"Suspicious Joins = "<<suspicious_breakpoints<<endl;
    cout<<"Percent Suspicious Joins = "<<suspicious_breakpoints*100.0/total_breakpoints<<endl;
	return 0;
}
