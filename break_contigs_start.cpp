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

unordered_map<string,vector<long>> contig2breakpoints;
unordered_map<string,long> contig_length;
void load_breakpoints(string file)
{
	string line;
	ifstream bfile(getCharExpr(file));
	string contig;
	long pos;
	while(getline(bfile,line))
	{
		vector<long> breaks;
		istringstream iss(line);
		iss >> contig;
		while(iss >> pos)
		{
			breaks.push_back(pos);
			//cout<<pos<<endl;
		}
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
void get_MAD_complete(string contig, vector<int>& count)
{
	if(contig2breakpoints.find(contig) == contig2breakpoints.end())
	{
		return;
	}
//	if(contig2breakpoints[contig].size() < 5)
//	{
//		cout<<"============================="<<endl;
//		cout<<contig<<endl;
//		for(int ind = 0; ind < contig2breakpoints[contig].size();ind++)
//		{
//			long pos = contig2breakpoints[contig][ind];
//			cout<<"POS = "<<pos<<" RATIO = "<<inspect(pos,count)<<endl;
//		}
//	}

//	cout<<"============================="<<endl;
	double median = 0;
	//unused variable
	//double mad = 0;
	long sz = count.size();
	vector<int> original = count;
	double sum = std::accumulate(original.begin(), original.end(), 0.0);
	double mean = sum / original.size();
	std::vector<double> diff(original.size());
	std::transform(original.begin(), original.end(), diff.begin(),
	std::bind2nd(std::minus<double>(), mean));
	//unused variables
  //double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
	//double stdev = std::sqrt(sq_sum / original.size());

	sort(count.begin(),count.end());

	// ints = curr_part.size();
	if(sz %2 == 0)
	{
		median = (count[sz/2] + count[sz/2-1])/2;
	}
	else
	{
		median = count[sz/2];
	}
	for(long i = 0;i < sz;i++)
	{
		count[i] = abs(count[i] - median);
	}
	sort(count.begin(),count.end());
	//mad = count[sz/2];

	vector<int>join_coverage;

//	cout<<"CONTIG "<<contig<<" MEAN = "<<mean<<" STD = "<<stdev<<" MEDIAN = "<<median<<" MAD = "<<mad<<endl;
	for(unsigned int ind = 0; ind < contig2breakpoints[contig].size();ind++)
	{
		long pos = contig2breakpoints[contig][ind];
		//cout<<"COVERAGE AT POSITION "<<pos <<" = "<<original[pos]<<endl;
		join_coverage.push_back(original[pos]);
	}
	for(unsigned int i = 1; i < join_coverage.size() - 1; i++)
	{
		if(join_coverage[i] < join_coverage[i-1] && join_coverage[i] < join_coverage[i+1])
		{
			cout<<contig<<"\t"<<contig2breakpoints[contig][i]<<endl;
		}
	}
}

void get_MAD(string contig, vector<int>& count)
{
	if(contig2breakpoints.find(contig) == contig2breakpoints.end())
	{
		return;
	}
	double median = 0;
	double mad = 0;
	long sz = count.size();
	vector<int> original = count;
	//unused variable
	// long check_index = 0;
	//cout<<"SIZE = "<<sz<<endl;

	for(unsigned int ind = 0; ind < contig2breakpoints[contig].size();ind++)
	{
		long bound = 5000000;
		long win_size;
		/*
		Check if window size is smaller than 5 mb
		*/
		long pos = contig2breakpoints[contig][ind];
		if(pos < bound/2)
		{
			win_size = pos;
			//cout<<"LOWER"<<endl;
		}
		else if(sz - pos < bound/2)
			{
				win_size = sz - pos -1;
				//cout<<"UPPER"<<endl;
			}
			else
			{
				//cout<<"MID"<<endl;
			}
		long start = pos - win_size;
		long end = pos + win_size;
		//cout<<"WINDOW Size = "<<start<<" - "<<end<<endl;
		vector<int> curr_part = slice(count,start,end);

		//COmpute mean and stdev
		double sum = std::accumulate(curr_part.begin(), curr_part.end(), 0.0);
		double mean = sum / curr_part.size();
		std::vector<double> diff(curr_part.size());
		std::transform(curr_part.begin(), curr_part.end(), diff.begin(),
        std::bind2nd(std::minus<double>(), mean));
		double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
		double stdev = std::sqrt(sq_sum / curr_part.size());

		//compute median and MAD
		sort(curr_part.begin(),curr_part.end());

		int curr_size = curr_part.size();
		if(curr_size  %2 == 0)
		{
			median = (curr_part[curr_size/2] + curr_part[curr_size/2-1])/2;
		}
		else
		{
			median = curr_part[curr_size/2];
		}
		for(long i = 0;i < curr_size;i++)
		{
			curr_part[i] = abs(curr_part[i] - median);
		}
		sort(curr_part.begin(),curr_part.end());
		mad = curr_part[curr_size/2];
		cout<<"MEDIAN = "<<median<<", MAD = "<<mad<<" MEAN = "<<mean<<" STDEV = "<<stdev<<endl;
		if(original[pos] < median - 3*mad)
		{
			cout<<"CHECKING SCAFFOLD "<<contig<<" AT POSITION "<<pos<<" YES; COV = "<<original[pos]<<endl;;
		}
		else
		{
			cout<<"CHECKING SCAFFOLD "<<contig<<" AT POSITION "<<pos<<" NO; COV = "<<original[pos]<<endl;;
		}

	}
}


int main(int argc, char *argv[])
{
	cmdline::parser p;
	p.add<string>("alignment", 'a', "bed file for alignment", true, "");
	//p.add<string>("outputdir", 'd', "coordinate output file", true, "");
	//p.add<string>("breakpoints", 'b', "breakpoints", true, "");
	p.add<string>("contiglen", 'l', "length of contigs", true, "");
	//p.add<int>("iteration",'i',"Iteration number",true,0);
    p.add<int>("min_size",'s',"Minimum mate pair separation for error findng",true,0);
    p.parse_check(argc, argv);
	string line;
	ifstream lenfile(getCharExpr(p.get<string>("contiglen")));
  // unused variabble
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
    //int iteration = p.get<int>("iteration");
	lenfile.close();
	//load_breakpoints(p.get<string>("breakpoints"));
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
		//unused variable
		// char strand;
		// long flag;
		long start,end;

		if(prev_line == "")
		{
			prev_line = line;
			continue;
		}
		istringstream iss(line);
		if(!(iss >> contig >> start >> end >> read))
			break;
		// if(contig == "000008F")
		// {
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
           long span_start = 0, span_end = 0;
           if(prev_end <= start)
           {
               span_start = prev_start;
               span_end = end;
           }
           else
           {
               if(prev_end >= start && prev_end <= end)
               {
                   span_start = prev_start;
                   span_end = end;
               }
               else
               {
                   if(end <= prev_start)
                   {
                       span_start = start;
                       span_end = prev_end;
                   }
                   else
                   {
                       if(end >= prev_start && end <= prev_end)
                       {
                           span_start = start;
                           span_end = prev_end;
                       }
                   }
               }
           }
           if(span_start >= 0 && span_end < contig_length[contig])
           {
                if(span_end - span_start <= contig2cutoff[contig])
                {
                    contig2coverage[contig].at(span_start) += 1;
                    contig2coverage[contig].at(span_end+1) -= 1;
                }

           } //contig2coverage[contig] = cov;
		}
		prev_contig = contig;
		prev_start = start;
		prev_end = end;
		prev_read = read;
	}
	bedfile.close();
	cerr<<"bedfile loaded"<<endl;

  // unused variables
  // int total_breakpoints = 0;
  // int suspicious_breakpoints = 0;

	for(unordered_map<string,vector<int> > :: iterator it = contig2coverage.begin(); it != contig2coverage.end();++it)
	{
		//cout<<"Testing contig " << it->first<<endl;
        vector<int> cov = contig2coverage[it->first];
		for(unsigned long i = 1; i< cov.size();i++)
		{
			cov[i] += cov[i-1];
			//cout<<(int)cov[i]<<endl;
		}
		contig2coverage[it->first] = cov;
		//get_MAD_complete(it->first,cov);
        string contig = it->first;

        /*
         *Delta array will store 1 if coverage is less than threshold , -1 otherwise
         */

        double average = accumulate(cov.begin(),cov.end(),0.0)/cov.size();
        vector<long> positions;
        for(int div = 5; div <= 15; div++)
        {
            double cutoff = average/div;
            vector<int>delta;
            //cout<<"Cutoff = " << cutoff<<endl;
            for(unsigned long i = 0; i < cov.size();i++)
            {
                if(cov[i] < cutoff)
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
            long  start = sz/100, end = sz/100,s = sz/100;
            for(long i = sz/100 ;i < 9*sz/10;i++)
            {
                max_ending_here += delta.at(i);
                if(max_so_far < max_ending_here)
                {
                    max_so_far = max_ending_here;
                    start = s;
                    end = i;
                }
                if(max_ending_here < 0)
                {
                    max_ending_here = 0;
                    s = i + 1;
                }
            }
            if(start >= sz/100 + 5000 && end <= 9*sz/10 - 5000)
            {
                positions.push_back((start+end)/2);
            }
            //cout<<"Possible Misassembly with cutoff "<<cutoff<<" in "<<contig<<" at "<<start<<"\t"<<end<<endl;
        }
        if(positions.size() > 0)
        {
            long consensus_pos = std::accumulate(positions.begin(),positions.end(),0.0)/positions.size();
            cout<<contig<<"\t"<<consensus_pos<<endl;
        }
    }
	return 0;
}
