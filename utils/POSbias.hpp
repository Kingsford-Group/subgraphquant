#ifndef __POSbias_T__
#define __POSbias_T__

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
#include "spline.h"

using namespace std;

void ReadPosBiasFile(string gzip_file, uint32_t& numModels, vector<uint32_t>& lenBounds, vector< vector<double> >& masses);


class PosBiasModel_t{
private:
	vector<uint32_t> lenBounds;
	vector< vector<double> > obs5models;
	vector< vector<double> > obs3models;
	vector< vector<double> > exp5models;
	vector< vector<double> > exp3models;
	const vector<double> positionBins{{0.02, 0.04, 0.06, 0.08, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.92, 0.94, 0.96, 0.98, 1.0}};
public:
	uint32_t numModels;
	vector< tk::spline > obs5Splines;
	vector< tk::spline > obs3Splines;
	vector< tk::spline > exp5Splines;
	vector< tk::spline > exp3Splines;

	PosBiasModel_t(){};
	PosBiasModel_t(string obs5_gzip_file, string obs3_gzip_file, string exp5_gzip_file, string exp3_gzip_file);

	PosBiasModel_t(uint32_t numModels, vector<uint32_t> lenBounds, vector< vector<double> > obs5models, vector< vector<double> > obs3models, 
		vector< vector<double> > exp5models, vector< vector<double> > exp3models):
		numModels(numModels), lenBounds(lenBounds), obs5models(obs5models), obs3models(obs3models), exp5models(exp5models), exp3models(exp3models) 
		{
			assert((int32_t)obs5models.size()==numModels);
			assert((int32_t)obs3models.size()==numModels);
			assert((int32_t)exp5models.size()==numModels);
			assert((int32_t)exp3models.size()==numModels);
			ProcessSplines();
		};

	void ProcessSplines();

	int32_t lenClass(uint32_t len) const {
		vector<uint32_t>::const_iterator itbegin=lenBounds.begin();
		vector<uint32_t>::const_iterator itend=lenBounds.end();
		vector<uint32_t>::const_iterator ub=upper_bound(itbegin, itend, len);
		return min((int32_t)lenBounds.size()-1, (int32_t)distance(itbegin, ub));
	};

};



#endif