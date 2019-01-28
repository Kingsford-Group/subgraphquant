#ifndef __utils__
#define __utils__

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <string>
#include "boost/algorithm/string.hpp"


using namespace std;


static std::map<char,char> Nucleotide={{'A','T'},{'C','G'},{'G','C'},{'T','A'},{'R','Y'},{'Y','R'},{'S','W'},{'W','S'},{'K','M'},{'M','K'},{'B','V'},{'V','B'},{'D','H'},{'H','D'}, {'N','N'}, {'.','.'},{'-','-'}};


inline string ReverseComplement(string::iterator itbegin, string::iterator itend)
{
	string rcseq="";
	for(string::iterator it=itbegin; it!=itend; it++){
		rcseq += Nucleotide[toupper(*it)];
	}
	reverse(rcseq.begin(), rcseq.end());
	return rcseq;
};


inline string ReverseComplement(string::const_iterator itbegin, string::const_iterator itend)
{
	string rcseq="";
	for(string::const_iterator it=itbegin; it!=itend; it++){
		rcseq += Nucleotide[toupper(*it)];
	}
	reverse(rcseq.begin(), rcseq.end());
	return rcseq;
};


inline void WriteCorrection(string outputfile, const vector< vector<double> >& corrections, const vector<string>& transnames)
{
	ofstream output(outputfile, ios::out | ios::binary);
	int32_t numtrans=corrections.size();
	output.write(reinterpret_cast<char*>(&numtrans), sizeof(int32_t));
	for(int32_t i = 0; i < numtrans; i++){
		string name = transnames[i];
		vector<double> bias = corrections[i];
		int32_t namelen = name.size();
		int32_t seqlen=bias.size();
		output.write(reinterpret_cast<char*>(&namelen), sizeof(int32_t));
		output.write(reinterpret_cast<char*>(&seqlen), sizeof(int32_t));
		output.write(name.c_str(), namelen*sizeof(char));
		output.write(reinterpret_cast<char*>(bias.data()), seqlen*sizeof(double));
	}
	output.close();
};


inline void WriteMatrix(string outputfile, const vector< vector< tuple<int32_t,int32_t,double> > >& corrections, const vector<string>& transnames)
{
	ofstream output(outputfile, ios::out | ios::binary);
	int32_t numtrans=corrections.size();
	output.write(reinterpret_cast<char*>(&numtrans), sizeof(int32_t));
	for(int32_t i = 0; i < numtrans; i++){
		string name = transnames[i];
		vector< tuple<int32_t,int32_t,double> > bias = corrections[i];
		int32_t namelen = name.size();
		int32_t seqlen=bias.size();
		output.write(reinterpret_cast<char*>(&namelen), sizeof(int32_t));
		output.write(reinterpret_cast<char*>(&seqlen), sizeof(int32_t));
		output.write(name.c_str(), namelen*sizeof(char));
		for (int32_t j = 0; j < bias.size(); j++) {
			output.write(reinterpret_cast<char*>(&get<0>(bias[j])), sizeof(int32_t));
			output.write(reinterpret_cast<char*>(&get<1>(bias[j])), sizeof(int32_t));
			output.write(reinterpret_cast<char*>(&get<2>(bias[j])), sizeof(double));
		}
	}
	output.close();
};


inline void ReadSalmonTPM(string salmonquantfile, const map< string,int32_t >& TransIndex, vector<double>& expr, double minvalue = 1e-4)
{
	// clear variable
	expr.clear();
	expr.assign(TransIndex.size(), 0);
	// read quantification output
	ifstream input(salmonquantfile);
	string line;
	int32_t linecount = 0;
	while (getline(input, line)) {
		linecount ++;
		if (linecount == 1)
			continue;
		vector<string> strs;
		boost::split(strs, line, boost::is_any_of("\t"));
		map< string,int32_t >::const_iterator it = TransIndex.find(strs[0]);
		if( it == TransIndex.cend() )
			continue;
		expr[it->second] = stod(strs[3]);
	}
	for (int32_t i = 0; i < expr.size(); i++)
		if (expr[i] <= 0)
			expr[i] = minvalue;
	input.close();
};

#endif