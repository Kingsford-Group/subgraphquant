#ifndef __GCbias_T__
#define __GCbias_T__

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


void ReadGCmatrix(string gzip_file, int32_t& dtype, Eigen::MatrixXd& Counts);


class GCBiasModel_t{
public:
	int32_t GCbin;
	int32_t Condbin;
	// -3 === -2 === -1 === 0 === 1
	// -3 to -1 is contextLeft, 1 is contextRight
	int32_t contextLeft;
	int32_t contextRight;
	Eigen::MatrixXd GCBiasRatio;

public:
	GCBiasModel_t(){};
	GCBiasModel_t(string obs_gc_gzip_file, string exp_gc_gzip_file, double maxRatio=1000.0, int32_t contextLeft=3, int32_t contextRight=1);
};

#endif
