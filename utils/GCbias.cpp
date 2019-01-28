#include "GCbias.hpp"


void ReadGCmatrix(string gzip_file, int32_t& dtype, Eigen::MatrixXd& Counts)
{
	// gzip file handle
	boost::iostreams::filtering_istream fpin;
	fpin.push(boost::iostreams::gzip_decompressor());
	fpin.push(boost::iostreams::file_source(gzip_file, std::ios_base::in | std::ios_base::binary));

	// get number of rows and columns of gc matrix
	fpin.read((char*)&dtype, sizeof(int32_t));
	typename Eigen::MatrixXd::Index rows, cols;
	fpin.read((char*)&rows, sizeof(typename Eigen::MatrixXd::Index));
	fpin.read((char*)&cols, sizeof(typename Eigen::MatrixXd::Index));

	vector<double> modelTotals;
	modelTotals.resize(rows, 0);
	fpin.read((char*)&modelTotals[0], rows*sizeof(double));

	Counts=Eigen::MatrixXd::Zero(rows, cols);
	fpin.read((char*)(Counts.data()), rows*cols*sizeof(typename Eigen::MatrixXd::Scalar));
};


GCBiasModel_t::GCBiasModel_t(string obs_gc_gzip_file, string exp_gc_gzip_file, double maxRatio, int32_t contextLeft, int32_t contextRight) :
	contextLeft(contextLeft), contextRight(contextRight)
{
	time_t CurrentTime;
	string CurrentTimeStr;
	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout<<"["<<CurrentTimeStr.substr(0, CurrentTimeStr.size()-1)<<"] "<<"Loading GC bias model."<<endl;

	// reading exp and ops GC bias matrices
	int32_t dtype_obs;
	Eigen::MatrixXd GC_obs;
	ReadGCmatrix(obs_gc_gzip_file, dtype_obs, GC_obs);

	int32_t dtype_exp;
	Eigen::MatrixXd GC_exp;
	ReadGCmatrix(exp_gc_gzip_file, dtype_exp, GC_exp);

	// GC matrix size
	GCbin=GC_obs.cols();
	Condbin=GC_obs.rows();

	// construct GC multiplier
	GCBiasRatio=Eigen::MatrixXd::Ones(Condbin, GCbin);
	if (dtype_obs != dtype_exp) {
		// there is one matrix not in log-scale, find it and turn into log-scale
		if (dtype_obs != 0) {
			for (int32_t i = 0; i < Condbin; i++) {
				for (int32_t j = 0; j < GCbin; j++) {
					GC_obs(i,j) = log(GC_obs(i,j) + 1e-10);
				}
			}
		}
		else{
			for (int32_t i = 0; i < Condbin; i++) {
				for (int32_t j = 0; j < GCbin; j++) {
					GC_exp(i,j) = log(GC_exp(i,j) + 1e-10);
				}
			}
		}
	}

	double minRatio = 1.0/maxRatio;
	// if log-scaled
	if (dtype_obs == 1) {
		for(int32_t i=0; i<Condbin; i++) {
			for(int32_t j=0; j<GCbin; j++){
				GCBiasRatio(i,j)=exp(GC_obs(i,j)-GC_exp(i,j));
				if(GCBiasRatio(i,j)<minRatio)
					GCBiasRatio(i,j)=minRatio;
				else if(GCBiasRatio(i,j)>maxRatio)
					GCBiasRatio(i,j)=maxRatio;
			}
		}
	}
	else{
		for(int32_t i=0; i<Condbin; i++) {
			for(int32_t j=0; j<GCbin; j++){
				GCBiasRatio(i,j)=GC_obs(i,j)/GC_exp(i,j);
				if(GCBiasRatio(i,j)<minRatio)
					GCBiasRatio(i,j)=minRatio;
				else if(GCBiasRatio(i,j)>maxRatio)
					GCBiasRatio(i,j)=maxRatio;
			}
		}
	}

	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout << "[" << CurrentTimeStr.substr(0, CurrentTimeStr.size()-1) << "] " << "Finish.\n";
};

