#ifndef __SEQbias_T__
#define __SEQbias_T__

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <string>
#include "Eigen/Dense"
#include "boost/algorithm/string.hpp"
#include "boost/iostreams/filtering_stream.hpp"
#include "boost/iostreams/device/file.hpp"
#include "boost/iostreams/filter/gzip.hpp"

using namespace std;


void ReadSeqBias(string gzip_file, int32_t& contextLeft, int32_t& contextRight, vector<int32_t>& orders, 
	vector<int32_t>& shifts, vector<int32_t>& widths, Eigen::MatrixXd& Probs);


class SeqBiasModel_t{
public:
	int32_t contextLeft;
	int32_t contextRight;
	vector<int32_t> orders, shifts, widths;
	Eigen::MatrixXd Obs5Probs;
	Eigen::MatrixXd Obs3Probs;
	Eigen::MatrixXd Exp5Probs;
	Eigen::MatrixXd Exp3Probs;
public:
	SeqBiasModel_t(){};
	SeqBiasModel_t(string obs5_gzip_file, string obs3_gzip_file, string exp5_gzip_file, string exp3_gzip_file);
	
	SeqBiasModel_t(int32_t contextLeft, int32_t contextRight, vector<int32_t> orders, vector<int32_t> shifts, vector<int32_t> widths, 
		Eigen::MatrixXd Obs5Probs, Eigen::MatrixXd Obs3Probs, Eigen::MatrixXd Exp5Probs, Eigen::MatrixXd Exp3Probs):
		contextLeft(contextLeft), contextRight(contextRight), orders(orders), shifts(shifts), widths(widths),
		Obs5Probs(Obs5Probs), Obs3Probs(Obs5Probs), Exp5Probs(Exp5Probs), Exp3Probs(Exp3Probs) {};
};



#endif