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

#include "ThreadPool.h"
#include "cmdline.h"

using namespace std;

int nthread = std::thread::hardware_concurrency();
unordered_map<string,vector<int> > contig2coverage;

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

void get_MAD_complete(string contig)
{
	if(contig2breakpoints.find(contig) == contig2breakpoints.end())
	{
		return;
	}
	cout<<"============================="<<endl;
	double median = 0;
	double mad = 0;
	vector<int> count = contig2coverage[contig];
	long sz = count.size();
	vector<int> original = count;
	double sum = std::accumulate(original.begin(), original.end(), 0.0);
	double mean = sum / original.size();
	std::vector<double> diff(original.size());
	std::transform(original.begin(), original.end(), diff.begin(),
    std::bind2nd(std::minus<double>(), mean));
	double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
	double stdev = std::sqrt(sq_sum / original.size());

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
	mad = count[sz/2];

	vector<int>join_coverage;

	cout<<"CONTIG "<<contig<<" MEAN = "<<mean<<" STD = "<<stdev<<" MEDIAN = "<<median<<" MAD = "<<mad<<endl;
	for(int ind = 0; ind < contig2breakpoints[contig].size();ind++)
	{
		long pos = contig2breakpoints[contig][ind];
		cout<<"COVERAGE AT POSITION "<<pos <<" = "<<original[pos]<<endl;
		join_coverage.push_back(original[pos]);
	}
	for(int i = 1; i < join_coverage.size() - 1; i++)
	{
		if(join_coverage[i] < join_coverage[i-1] && join_coverage[i] < join_coverage[i+1])
		{
			cout<<"POSITION "<<contig2breakpoints[contig][i]<<" AMBIGUOUS"<<endl;
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
	long check_index = 0;
	//cout<<"SIZE = "<<sz<<endl;
	
	for(int ind = 0; ind < contig2breakpoints[contig].size();ind++)
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

	//First create a threadpool of nthreads
	ThreadPool POOL(nthread);

	cmdline::parser p;
	p.add<string>("alignment", 'a', "bed file for alignment", true, "");
	p.add<string>("outputdir", 'd', "coordinate output file", true, "");
	p.add<string>("breakpoints", 'b', "breakpoints", true, "");
	p.add<string>("contiglen", 'l', "length of contigs", true, "");
	p.parse_check(argc, argv);	
	string line;
	ifstream lenfile(getCharExpr(p.get<string>("contiglen")));
	while(getline(lenfile,line))
	{
		istringstream iss(line);
		string a;
		long n;
		iss >> a >> n;
		contig_length[a] = n;
	}

	lenfile.close();
	load_breakpoints(p.get<string>("breakpoints"));
	cout<<"loaded breakpoints"<<endl;
	ifstream bedfile(getCharExpr(p.get<string>("alignment")));
	
	
	string prev_line = "";
	string prev_contig="";
	long prev_start=-1;
	long prev_end=-1;
	string prev_read="";
	while(getline(bedfile,line))
	{
		string contig, read;
		char strand;
		long start,end,flag;
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
			// cout<<"here"<<endl;
			if(prev_end <= start)
			{
				//cout<<"here"<<endl;
				contig2coverage[contig].at(prev_start) += 1;
				contig2coverage[contig].at(end+1) -= 1;
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
	cerr<<"bedfile loaded"<<endl;

	std::vector< std::future<void> > results;


	for(unordered_map<string,vector<int> > :: iterator it = contig2coverage.begin(); it != contig2coverage.end();++it)
	{
		vector<int> cov = contig2coverage[it->first];
		for(long i = 1; i< cov.size();i++)
		{
			cov[i] += cov[i-1];
			//cout<<(int)cov[i]<<endl;
		}
		contig2coverage[it->first] = cov;
		//MULTITHREADED STUFF, CHECK CAREFULLY LATER
		results.emplace_back(POOL.enqueue(get_MAD_complete,it->first));
		for(auto && result: results)
      	{
      		result.get();
      	}
		//get_MAD_complete(it->first,cov);
		
		//cout<<stat.first<<"\t"<<stat.second<<endl;
	}
	
	// vector<int> cov = contig2coverage["scaffold_229"];
	// ofstream ofile1("scaffold_229.coverage");
	// for(long i = 0; i< cov.size();i += 10000)
	// {
	// 	//cout<<(int)cov[i]<<endl;
	// 	ofile1<<i<<"\t"<<cov[i]<<endl;
	// }
	// ofile1.close();

	// cov = contig2coverage["scaffold_30"];
	// ofstream ofile2("scaffold_30.coverage");
	// for(long i = 0; i< cov.size();i += 10000)
	// {
	// 	//cout<<(int)cov[i]<<endl;
	// 	ofile2<<i<<"\t"<<cov[i]<<endl;
	// }
	// ofile2.close();

	// cov = contig2coverage["scaffold_293"];
	// ofstream ofile3("scaffold_293.coverage");
	// for(long i = 0; i< cov.size();i += 10000)
	// {
	// 	//cout<<(int)cov[i]<<endl;
	// 	ofile3<<i<<"\t"<<cov[i]<<endl;
	// }
	// ofile3.close();

	// ofstream outputfile(getCharExpr(p.get<string>("outputdir")+"/alignment_unique_clean.bed"));
	// ofstream outputfile2(getCharExpr(p.get<string>("outputdir")+"/breakpoints"));
	//ofstream outputfile1("length_to_coverage");
	/*
	for(unordered_map<string,vector<pair<long,long> > >  :: iterator it = intervals.begin(); it != intervals.end(); ++it)
	{
		vector<pair<long,long> > longs = it->second;
		sort(longs.begin(),longs.end());
		intervals[it->first] = longs;
	}*/
	
	
	//const char* args[] = {"000001F"};
	//const char* args[] = {"chr6"};
	return 0;
}


