#ifndef __CombinedBias_T__
#define __CombinedBias_T__

#include <iostream>
#include <fstream>
#include <vector>
#include <tuple>
#include <map>
#include <algorithm>
#include <string>
#include <numeric>
#include <omp.h>
#include <mutex>
#include "Eigen/Dense"
#include "boost/algorithm/string.hpp"
#include "boost/iostreams/filtering_stream.hpp"
#include "boost/iostreams/device/file.hpp"
#include "boost/math/distributions/binomial.hpp"
#include "jellyfish/mer_dna.hpp"
#include "GCbias.hpp"
#include "SEQbias.hpp"
#include "POSbias.hpp"
#include "utils.hpp"


using namespace std;
using Mer = jellyfish::mer_dna_ns::mer_base_static<uint64_t, 4>;


void ReadFLD(string fld_gzip_file, vector<int32_t>& RawFLD, int32_t numbins=1001);


void FLDKDE(vector<int32_t>& RawFLD, vector<double>& FLD, int32_t kernel_n=10, double kernel_p=0.5);


void GetFLDbound(const vector<double>& FLD, int32_t& fldLow, int32_t& fldHigh);


vector<double> BiasEffLen_coverage(const string seq, const GCBiasModel_t& gcbias, const SeqBiasModel_t& seqbias, const vector<double>& FLD, int32_t fldLow, int32_t fldHigh);


vector<double> BiasEffLen_fragstart(const string seq, const GCBiasModel_t& gcbias, const SeqBiasModel_t& seqbias, const PosBiasModel_t& posbias, const vector<double>& FLD, int32_t fldLow, int32_t fldHigh);


vector< tuple<int32_t,int32_t,double> > BiasEffLen_matrix(const string seq, const GCBiasModel_t& gcbias, const SeqBiasModel_t& seqbias, const PosBiasModel_t& posbias, const vector<double>& FLD, int32_t fldLow, int32_t fldHigh);


void ProcessBias_coverage(const vector<string>& sequences, const GCBiasModel_t& gcbias, const SeqBiasModel_t& seqbias, const vector<double>& FLD,
	vector< vector<double> >& corrections, int32_t fldLow, int32_t fldHigh, int32_t threads=8);


void ProcessBias_fragstart(const vector<string>& sequences, const GCBiasModel_t& gcbias, const SeqBiasModel_t& seqbias, const PosBiasModel_t& posbias, const vector<double>& FLD,
	vector< vector<double> >& corrections, int32_t fldLow, int32_t fldHigh, int32_t threads=8);


void ProcessBias_matrix(const vector<string>& sequences, const GCBiasModel_t& gcbias, const SeqBiasModel_t& seqbias, const PosBiasModel_t& posbias, const vector<double>& FLD,
	vector< vector< tuple<int32_t,int32_t,double> > >& corrections, int32_t fldLow, int32_t fldHigh, int32_t threads=8);


void AdjustBias(vector< vector<double> >& corrections, const vector<string>& transnames,string salmonquantfile);

void AdjustBias(vector< vector< tuple<int32_t,int32_t,double> > >& corrections, const vector<string>& transnames,string salmonquantfile);

#endif
