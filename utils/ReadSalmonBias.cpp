#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <omp.h>
#include <mutex>
#include "boost/algorithm/string.hpp"
#include "boost/iostreams/filtering_stream.hpp"
#include "boost/iostreams/device/file.hpp"
#include "boost/math/distributions/binomial.hpp"
#include "Eigen/Dense"
#include "jellyfish/mer_dna.hpp"
#include "spline.h"

using namespace std;
using Mer = jellyfish::mer_dna_ns::mer_base_static<uint64_t, 4>;

string ReverseComplement(string::iterator itbegin, string::iterator itend){
	string rcseq="";
	for(string::iterator it=itbegin; it!=itend; it++){
		if((*it)=='G' || (*it)=='g')
			rcseq+='C';
		else if((*it)=='C' || (*it)=='c')
			rcseq+='G';
		else if((*it)=='A' || (*it)=='a')
			rcseq+='T';
		else if((*it)=='T' || (*it)=='t')
			rcseq+='A';
	}
	reverse(rcseq.begin(), rcseq.end());
	return rcseq;
};

void ReadSeqBias(boost::iostreams::filtering_istream& fpin, int32_t& contextLength, int32_t& contextLeft, int32_t& contextRight, vector<int32_t>& orders, 
	vector<int32_t>& shifts, vector<int32_t>& widths, Eigen::MatrixXd& Probs){

	fpin.read((char*)&contextLength, sizeof(int32_t));
	fpin.read((char*)&contextLeft, sizeof(int32_t));
	fpin.read((char*)&contextRight, sizeof(int32_t));

	orders.resize(contextLength);
	shifts.resize(contextLength);
	widths.resize(contextLength);
	fpin.read((char*)&orders[0], contextLength*sizeof(int32_t));
	fpin.read((char*)&shifts[0], contextLength*sizeof(int32_t));
	fpin.read((char*)&widths[0], contextLength*sizeof(int32_t));

	typename Eigen::MatrixXd::Index prows, pcols;
	fpin.read((char*)&prows, sizeof(typename Eigen::MatrixXd::Index));
	fpin.read((char*)&pcols, sizeof(typename Eigen::MatrixXd::Index));

	Probs=Eigen::MatrixXd::Zero(prows, pcols);
	fpin.read((char*)(Probs.data()), prows*pcols*sizeof(typename Eigen::MatrixXd::Scalar));
}

void ReadGCBias(boost::iostreams::filtering_istream& fpin, int32_t& dtype, Eigen::MatrixXd& Counts){
	fpin.read((char*)&dtype, sizeof(int32_t));
	typename Eigen::MatrixXd::Index rows, cols;
	fpin.read((char*)&rows, sizeof(typename Eigen::MatrixXd::Index));
	fpin.read((char*)&cols, sizeof(typename Eigen::MatrixXd::Index));

	vector<double> modelTotals;
	modelTotals.resize(rows, 0);
	Counts=Eigen::MatrixXd::Zero(rows, cols);
	fpin.read((char*)&modelTotals[0], rows*sizeof(double));

	fpin.read((char*)(Counts.data()), rows*cols*sizeof(typename Eigen::MatrixXd::Scalar));
};

void ReadPosBias(boost::iostreams::filtering_istream& fpin, uint32_t& numModels, vector<uint32_t>& lenBounds, vector< vector<double> >& masses){
	fpin.read((char*)(&numModels), sizeof(uint32_t));
	lenBounds.resize(numModels);
	masses.clear();
	for(int i=0; i<numModels; i++)
		fpin.read((char*)(&lenBounds[i]), sizeof(uint32_t));
	for(int i=0; i<numModels; i++){
		uint32_t modelLen=0;
		fpin.read((char*)(&modelLen), sizeof(uint32_t));
		vector<double> tmpmodel(modelLen, 0);
		fpin.read((char*)tmpmodel.data(), modelLen*sizeof(double));
		masses.push_back(tmpmodel);
	}
};

void ReadFLD(boost::iostreams::filtering_istream& fpin, vector<int32_t>& FLD){
	FLD.resize(1001, -1);
	fpin.read((char*)&FLD[0], 1001*sizeof(int32_t));
};

void FLDKDE(vector<int32_t>& RawFLD, vector<double>& FLD, int32_t kernel_n=10, double kernel_p=0.5){
	// initialize fragment length distribution FLD
	FLD.assign(RawFLD.size(), 0);
	// calculate binomial kernel
	boost::math::binomial bino(kernel_n, kernel_p);
	vector<double> kernel(kernel_n+1, 0);
	for(int32_t i=0; i<kernel_n+1; i++)
		kernel[i]=pdf(bino, i);
	// calculate FLD based on kernel
	for(int32_t i=0; i<RawFLD.size(); i++){
		if(RawFLD[i]==0)
			continue;
		int32_t offset=max(0, i-kernel_n/2);
		while(offset<=i+kernel_n/2 && offset<RawFLD.size()){
			FLD[offset]+=RawFLD[i]*kernel[offset-i+kernel_n/2];
			offset++;
		}
	}
	double sum=0;
	for(int32_t i=0; i<RawFLD.size(); i++){
		FLD[i]+=1e-8;
		sum+=FLD[i];
	}
	for(int32_t i=0; i<RawFLD.size(); i++)
		FLD[i]/=sum;
};

class GCBiasModel_t{
public:
	int32_t GCbin, Condbin;
	// -3 === -2 === -1 === 0 === 1
	// -3 to -1 is contextLeft, 1 is contextRight
	int32_t contextLeft, contextRight;
	Eigen::MatrixXd GCBiasRatio;
public:
	GCBiasModel_t(){};
	GCBiasModel_t(int32_t dtype1, Eigen::MatrixXd GCObs, int32_t dtype2, Eigen::MatrixXd GCExp, double maxRatio, int32_t contextLeft, int32_t contextRight):
		contextLeft(contextLeft), contextRight(contextRight)
	{
		GCbin=GCObs.cols();
		Condbin=GCObs.rows();
		GCBiasRatio=Eigen::MatrixXd::Ones(Condbin, GCbin);
		if(dtype1!=dtype2){
			if(dtype1==1){
				for(int i=0; i<Condbin; i++)
					for(int j=0; j<GCbin; j++)
						GCObs(i,j)=exp(GCObs(i,j));
			}
			if(dtype2==1){
				for(int i=0; i<Condbin; i++)
					for(int j=0; j<GCbin; j++)
						GCExp(i,j)=exp(GCExp(i,j));
			}
			// assertion for normalized
			for(int i=0; i<Condbin; i++){
				double sumObs=0, sumExp=0;
				for(int j=0; j<GCbin; j++){
					sumObs+=GCObs(i,j);
					sumExp+=GCExp(i,j);
				}
				assert(fabs(sumObs-1)<1e-4);
				assert(fabs(sumExp-1)<1e-4);
			}
		}
		double minRatio = 1.0/maxRatio;
		if(dtype1==1){
			for(int i=0; i<Condbin; i++)
				for(int j=0; j<GCbin; j++){
					GCBiasRatio(i,j)=exp(GCObs(i,j)-GCExp(i,j));
					if(GCBiasRatio(i,j)<minRatio)
						GCBiasRatio(i,j)=minRatio;
					else if(GCBiasRatio(i,j)>maxRatio)
						GCBiasRatio(i,j)=maxRatio;
				}
		}
		else{
			for(int i=0; i<Condbin; i++)
				for(int j=0; j<GCbin; j++){
					GCBiasRatio(i,j)=GCObs(i,j)/GCExp(i,j);
					if(GCBiasRatio(i,j)<minRatio)
						GCBiasRatio(i,j)=minRatio;
					else if(GCBiasRatio(i,j)>maxRatio)
						GCBiasRatio(i,j)=maxRatio;
				}
		}
	};
};

class SeqBiasModel_t{
public:
	int32_t contextLeft, contextRight;
	vector<int32_t> orders, shifts, widths;
	Eigen::MatrixXd Obs5Probs, Obs3Probs, Exp5Probs, Exp3Probs;
public:
	SeqBiasModel_t(){};
	SeqBiasModel_t(int32_t contextLeft, int32_t contextRight, vector<int32_t> orders, vector<int32_t> shifts, vector<int32_t> widths, 
		Eigen::MatrixXd Obs5Probs, Eigen::MatrixXd Obs3Probs, Eigen::MatrixXd Exp5Probs, Eigen::MatrixXd Exp3Probs):
		contextLeft(contextLeft), contextRight(contextRight), orders(orders), shifts(shifts), widths(widths),
		Obs5Probs(Obs5Probs), Obs3Probs(Obs5Probs), Exp5Probs(Exp5Probs), Exp3Probs(Exp3Probs) {};
};

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
	PosBiasModel_t(uint32_t numModels, vector<uint32_t> lenBounds, vector< vector<double> > obs5models, vector< vector<double> > obs3models, 
		vector< vector<double> > exp5models, vector< vector<double> > exp3models):
		numModels(numModels), lenBounds(lenBounds), obs5models(obs5models), obs3models(obs3models), exp5models(exp5models), exp3models(exp3models) 
		{
			assert((int)obs5models.size()==numModels);
			assert((int)obs3models.size()==numModels);
			assert((int)exp5models.size()==numModels);
			assert((int)exp3models.size()==numModels);
			ProcessSplines();
		};
	void ProcessSplines(){
		for(int i=0; i<numModels; i++){
			vector<double> splineMass(obs5models[i].size()+2);
			vector<double> splineBins(obs5models[i].size()+2);
			// observed 5 model
			splineBins[0]=0;
			splineMass[0]=obs5models[i][0];
			for(int j=0; j<obs5models[i].size(); j++){
				splineBins[j+1]=positionBins[j]-0.01;
				splineMass[j+1]=obs5models[i][j]/(1+obs5models[i][0]+obs5models[i].back());
			}
			splineBins.back()=1;
			splineMass.back()=obs5models[i].back();
			tk::spline s_obs5;
			s_obs5.set_points(splineBins, splineMass);
			obs5Splines.push_back(s_obs5);
			// observed 3 model
			splineBins[0]=0;
			splineMass[0]=obs3models[i][0];
			for(int j=0; j<obs3models[i].size(); j++){
				splineBins[j+1]=positionBins[j]-0.01;
				splineMass[j+1]=obs3models[i][j]/(1+obs3models[i][0]+obs3models[i].back());
			}
			splineBins.back()=1;
			splineMass.back()=obs3models[i].back();
			tk::spline s_obs3;
			s_obs3.set_points(splineBins, splineMass);
			obs3Splines.push_back(s_obs3);
			// expected 5 model
			splineBins[0]=0;
			splineMass[0]=exp5models[i][0];
			for(int j=0; j<exp5models[i].size(); j++){
				splineBins[j+1]=positionBins[j]-0.01;
				splineMass[j+1]=exp5models[i][j]/(1+exp5models[i][0]+exp5models[i].back());
			}
			splineBins.back()=1;
			splineMass.back()=exp5models[i].back();
			tk::spline s_exp5;
			s_exp5.set_points(splineBins, splineMass);
			exp5Splines.push_back(s_exp5);
			// expected 3 model
			splineBins[0]=0;
			splineMass[0]=exp3models[i][0];
			for(int j=0; j<exp3models[i].size(); j++){
				splineBins[j+1]=positionBins[j]-0.01;
				splineMass[j+1]=exp3models[i][j]/(1+exp3models[i][0]+exp3models[i].back());
			}
			splineBins.back()=1;
			splineMass.back()=exp3models[i].back();
			tk::spline s_exp3;
			s_exp3.set_points(splineBins, splineMass);
			exp3Splines.push_back(s_exp3);
		}
	};
	int lenClass(uint32_t len) const{
		vector<uint32_t>::const_iterator itbegin=lenBounds.begin();
		vector<uint32_t>::const_iterator itend=lenBounds.end();
		vector<uint32_t>::const_iterator ub=upper_bound(itbegin, itend, len);
		return min((int)lenBounds.size()-1, (int)distance(itbegin, ub));
	};
};

void GetFLDbound(const vector<double>& FLD, int& fldLow, int& fldHigh){
	double quantileCutoffLow = 0.005;
	double FragCount=0;
	for(int i=0; i<FLD.size(); i++)
		FragCount+=FLD[i];
	bool lb=false, ub=false;
	double fldcummulative=0;
	for(int i=0; i<FLD.size(); i++){
		fldcummulative+=FLD[i];
		if(!lb && fldcummulative/FragCount > quantileCutoffLow){
			fldLow=i;
			lb=true;
		}
		if(!ub && fldcummulative/FragCount > 1-quantileCutoffLow){
			fldHigh=i-1;
			ub=true;
		}
		if(lb && ub)
			break;
	}
	cout<<"estimated fragment length lb = "<<fldLow<<"\t ub = "<<fldHigh<<endl;
};

vector<double> BiasCorrectTrans(string seq, const GCBiasModel_t& gcbias, const SeqBiasModel_t& seqbias, const PosBiasModel_t& posbias, const vector<double>& FLD, int fldLow, int fldHigh){
	vector<double> Correction(seq.size(), 0);

	// process GC condition
	vector<int> GCCondFW(seq.size(), 0);
	vector<int> GCCondRC(seq.size(), 0);
	vector<int> GCWindowFW(seq.size(), 0);
	vector<int> GCWindowRC(seq.size(), 0);
	for(int i=0; i<seq.size(); i++){
		int gcFWcount=0, gcRCcount=0;
		for(int j=i-gcbias.contextLeft; j<i+gcbias.contextRight+1; j++){
			if(j>=0 && j<seq.size())
				if(seq[j]=='G' || seq[j]=='g' || seq[j]=='C' || seq[j]=='c')
					gcFWcount++;
		}
		GCCondFW[i]=gcFWcount;
		GCWindowFW[i]=(min(i+gcbias.contextRight+1, (int)seq.size())-max(0, i-gcbias.contextLeft));
		for(int j=i+gcbias.contextLeft; j>i-gcbias.contextRight-1; j--){
			if(j>=0 && j<seq.size())
				if(seq[j]=='G' || seq[j]=='g' || seq[j]=='C' || seq[j]=='c')
					gcRCcount++;
		}
		GCCondRC[i]=gcRCcount;
		GCWindowRC[i]=(min((int)seq.size(), i+gcbias.contextLeft+1)-max(0, i-gcbias.contextRight));
	}
	// process GC content vector
	vector<int> GCRawCount(seq.size()+1, 0);
	int cummulative=0;
	for(int i=0; i<seq.size(); i++){
		if(seq[i]=='G' || seq[i]=='g' || seq[i]=='C' || seq[i]=='c'){
			cummulative++;
		}
		GCRawCount[i+1]=cummulative;
	}

	// process 5' seqbias positional ratio
	vector<double> Seq5Ratio(seq.size(), 1);
	vector<double> Seq3Ratio(seq.size(), 1);
	Mer mer;
	mer.k(seqbias.contextLeft+seqbias.contextRight+1);
	mer.from_chars(seq.c_str());
	for(int i=seqbias.contextLeft; i<seq.size()-seqbias.contextRight; i++){
		double obsvalue=0, expvalue=0;
		for(int j=0; j<seqbias.contextLeft+seqbias.contextRight+1; j++){
			int idx=mer.get_bits(seqbias.shifts[j], seqbias.widths[j]);
			obsvalue+=seqbias.Obs5Probs(idx, j);
			expvalue+=seqbias.Exp5Probs(idx, j);
		}
		Seq5Ratio[i]=exp(obsvalue-expvalue);
		mer.shift_left(seq[i+seqbias.contextRight+1]);
	}
	// process 3' seqbias positional ratio
	string rcseq=ReverseComplement(seq.begin(), seq.end());
	mer.from_chars(rcseq.c_str());
	for(int i=seqbias.contextLeft; i<rcseq.size()-seqbias.contextRight; i++){
		double obsvalue=0, expvalue=0;
		for(int j=0; j<seqbias.contextLeft+seqbias.contextRight+1; j++){
			int idx=mer.get_bits(seqbias.shifts[j], seqbias.widths[j]);
			obsvalue+=seqbias.Obs3Probs(idx, j);
			expvalue+=seqbias.Exp3Probs(idx, j);
		}
		Seq3Ratio[i]=exp(obsvalue-expvalue);
		mer.shift_left(rcseq[i+seqbias.contextRight+1]);
	}
	reverse(Seq3Ratio.begin(), Seq3Ratio.end());

	// process 5' and 3' posbias
	vector<double> Pos5Ratio(seq.size(), 1);
	vector<double> Pos3Ratio(seq.size(), 1);
	int li=posbias.lenClass((int)seq.size());
	for(int i=0; i<seq.size()-(seqbias.contextLeft+seqbias.contextRight+1); i++){
		double fracP=1.0*i/seq.size();
		double obs5=max(0.001, posbias.obs5Splines[li](fracP));
		double obs3=max(0.001, posbias.obs3Splines[li](fracP));
		double exp5=max(0.001, posbias.exp5Splines[li](fracP));
		double exp3=max(0.001, posbias.exp3Splines[li](fracP));
		Pos5Ratio[i]=obs5/exp5;
		Pos3Ratio[i]=obs3/exp3;
		// cout<<"POS\t"<<seq.size()<<"\t"<<li<<"\t"<<i<<"\t"<<obs3<<"\t"<<exp3<<endl;
	}

	// calculate correction
	double FragCount=0;
	int minFLDpossible=(seq.size()<FLD.size())?1:fldLow;
	int maxFLDpossible=(seq.size()<FLD.size())?seq.size():fldHigh;
	for(int i=0; i<maxFLDpossible; i++)
		FragCount+=FLD[i];
	for(int k=minFLDpossible; k<maxFLDpossible; k++){
		double flMassTotal=0;
		for(int i=0; i<seq.size()-k-1; i++){
			// seqBias
			double seq5factor=Seq5Ratio[i];
			double seq3factor=Seq3Ratio[i+k-1];
			// gcBias
			double gcfactor;
			double gccondition=1.0*(GCCondFW[i]+GCCondRC[i+k-1])/(GCWindowFW[i]+GCWindowRC[i+k-1]);
			int32_t gccondbin=(int)(gccondition*gcbias.Condbin);
			if(gccondbin==gcbias.Condbin)
				gccondbin--;
			double gcfraction=1.0*(GCRawCount[i+k]-GCRawCount[i])/k;
			int32_t gcbin=(int)(gcbias.GCbin*lrint(gcfraction*100)/100.0);
			if(gcbin==gcbias.GCbin)
				gcbin--;
			gcfactor=gcbias.GCBiasRatio(gccondbin, gcbin);
			// posBias
			double posfactor=Pos5Ratio[i]*Pos3Ratio[i+k-1];
			Correction[i]+=1.0*(FLD[k])*seq5factor*seq3factor*gcfactor*posfactor/FragCount;
			// Correction[i+k/2]+=1.0*(FLD[k])*seq5factor*seq3factor*gcfactor*posfactor/FragCount;
			// flMassTotal+=seq5factor*seq3factor*gcfactor;
			// if(k==110){
			// 	cout<<"SEQ\t"<<i<<"\t"<<(i+k-1)<<"\t"<<seq5factor<<"\t"<<seq3factor<<endl;
			// 	cout<<"GC\t"<<i<<"\t"<<(i+k-1)<<"\t"<<lrint(gcfraction*100)<<"\t"<<lrint(gccondition*100)<<"\t"<<gcfactor<<endl;
			// 	cout<<"POS\t"<<i<<"\t"<<(i+k-1)<<"\t"<<Pos5Ratio[i]<<"\t"<<Pos3Ratio[i+k-1]<<endl;
			// }
		}
		// cout<<k<<"\t"<<((FLD[k])/FragCount)<<"\t"<<flMassTotal<<endl;
	}
	return Correction;
};

vector<double> BiasCorrectTrans_wopos(string seq, const GCBiasModel_t& gcbias, const SeqBiasModel_t& seqbias, const vector<double>& FLD, int fldLow, int fldHigh){
	vector<double> Correction(seq.size(), 0);

	// process GC condition
	vector<int> GCCondFW(seq.size(), 0);
	vector<int> GCCondRC(seq.size(), 0);
	vector<int> GCWindowFW(seq.size(), 0);
	vector<int> GCWindowRC(seq.size(), 0);
	for(int i=0; i<seq.size(); i++){
		int gcFWcount=0, gcRCcount=0;
		for(int j=i-gcbias.contextLeft; j<i+gcbias.contextRight+1; j++){
			if(j>=0 && j<seq.size())
				if(seq[j]=='G' || seq[j]=='g' || seq[j]=='C' || seq[j]=='c')
					gcFWcount++;
		}
		GCCondFW[i]=gcFWcount;
		GCWindowFW[i]=(min(i+gcbias.contextRight+1, (int)seq.size())-max(0, i-gcbias.contextLeft));
		for(int j=i+gcbias.contextLeft; j>i-gcbias.contextRight-1; j--){
			if(j>=0 && j<seq.size())
				if(seq[j]=='G' || seq[j]=='g' || seq[j]=='C' || seq[j]=='c')
					gcRCcount++;
		}
		GCCondRC[i]=gcRCcount;
		GCWindowRC[i]=(min((int)seq.size(), i+gcbias.contextLeft+1)-max(0, i-gcbias.contextRight));
	}
	// process GC content vector
	vector<int> GCRawCount(seq.size()+1, 0);
	int cummulative=0;
	for(int i=0; i<seq.size(); i++){
		if(seq[i]=='G' || seq[i]=='g' || seq[i]=='C' || seq[i]=='c'){
			cummulative++;
		}
		GCRawCount[i+1]=cummulative;
	}

	// process 5' seqbias positional ratio
	vector<double> Seq5Ratio(seq.size(), 1);
	vector<double> Seq3Ratio(seq.size(), 1);
	Mer mer;
	mer.k(seqbias.contextLeft+seqbias.contextRight+1);
	mer.from_chars(seq.c_str());
	for(int i=seqbias.contextLeft; i<seq.size()-seqbias.contextRight; i++){
		double obsvalue=0, expvalue=0;
		for(int j=0; j<seqbias.contextLeft+seqbias.contextRight+1; j++){
			int idx=mer.get_bits(seqbias.shifts[j], seqbias.widths[j]);
			obsvalue+=seqbias.Obs5Probs(idx, j);
			expvalue+=seqbias.Exp5Probs(idx, j);
		}
		Seq5Ratio[i]=exp(obsvalue-expvalue);
		mer.shift_left(seq[i+seqbias.contextRight+1]);
	}
	// process 3' seqbias positional ratio
	string rcseq=ReverseComplement(seq.begin(), seq.end());
	mer.from_chars(rcseq.c_str());
	for(int i=seqbias.contextLeft; i<rcseq.size()-seqbias.contextRight; i++){
		double obsvalue=0, expvalue=0;
		for(int j=0; j<seqbias.contextLeft+seqbias.contextRight+1; j++){
			int idx=mer.get_bits(seqbias.shifts[j], seqbias.widths[j]);
			obsvalue+=seqbias.Obs3Probs(idx, j);
			expvalue+=seqbias.Exp3Probs(idx, j);
		}
		Seq3Ratio[i]=exp(obsvalue-expvalue);
		mer.shift_left(rcseq[i+seqbias.contextRight+1]);
	}
	reverse(Seq3Ratio.begin(), Seq3Ratio.end());

	// calculate correction
	double FragCount=0;
	int minFLDpossible=(seq.size()<FLD.size())?1:fldLow;
	int maxFLDpossible=(seq.size()<FLD.size())?seq.size():fldHigh;
	for(int i=0; i<maxFLDpossible; i++)
		FragCount+=FLD[i];
	for(int k=minFLDpossible; k<maxFLDpossible; k++){
		double flMassTotal=0;
		for(int i=0; i<seq.size()-k-1; i++){
			// seqBias
			double seq5factor=Seq5Ratio[i];
			double seq3factor=Seq3Ratio[i+k-1];
			// gcBias
			double gcfactor;
			double gccondition=1.0*(GCCondFW[i]+GCCondRC[i+k-1])/(GCWindowFW[i]+GCWindowRC[i+k-1]);
			int32_t gccondbin=(int)(gccondition*gcbias.Condbin);
			if(gccondbin==gcbias.Condbin)
				gccondbin--;
			double gcfraction=1.0*(GCRawCount[i+k]-GCRawCount[i])/k;
			int32_t gcbin=(int)(gcbias.GCbin*lrint(gcfraction*100)/100.0);
			if(gcbin==gcbias.GCbin)
				gcbin--;
			gcfactor=gcbias.GCBiasRatio(gccondbin, gcbin);
			Correction[i]+=1.0*(FLD[k])*seq5factor*seq3factor*gcfactor/FragCount;
		}
	}
	return Correction;
};

void ReadTransSequence(string filename, vector<string>& TransSequences, vector<string>& TransNames){
	ifstream input(filename);
	string line;
	string curname="";
	string curseq="";
	while(getline(input, line)){
		if(line[0]=='>'){
			if(curname!=""){
				TransSequences.push_back(curseq);
				TransNames.push_back(curname);
			}
			vector<string> strs;
			boost::split(strs, line, boost::is_any_of(" |"));
			curname=strs[0].substr(1);
			curseq="";
		}
		else
			curseq+=line;
	}
	TransSequences.push_back(curseq);
	TransNames.push_back(curname);
};

void ReadSalmonCov(string filename, map<string, double>& SalmonCov){
	SalmonCov.clear();
	ifstream input(filename);
	string line;
	int linecount=0;
	while(getline(input, line)){
		linecount++;
		if(linecount==1)
			continue;
		vector<string> strs;
		boost::split(strs, line, boost::is_any_of("\t"));
		double cov=stod(strs[4])/stod(strs[1]);
		SalmonCov[strs[0]]=cov;
	}
	input.close();
};

void CorrectNWrite(string filename, vector<string>& TransNames, vector<string>& TransSequences, map<string, double>& SalmonCov, GCBiasModel_t& gcbias, 
	SeqBiasModel_t& seqbias, PosBiasModel_t& posbias, vector<double>& FLD, int fldLow, int fldHigh, double threshold=1){

	boost::iostreams::filtering_ostream fp1;
	fp1.push(boost::iostreams::file_sink(filename, std::ios_base::out | std::ios_base::binary));
	ofstream fp2(filename+"names");
	for(int i=0; i<TransNames.size(); i++){
		map<string, double>::iterator it=SalmonCov.find(TransNames[i]);
		if(it!=SalmonCov.end() && it->second>threshold){
			vector<double> Correction=BiasCorrectTrans(TransSequences[i], gcbias, seqbias, posbias, FLD, fldLow, fldHigh);
			fp2<<TransNames[i]<<endl;
			int64_t length=Correction.size();
			fp1.write(reinterpret_cast<char*>(&length), sizeof(int64_t));
			fp1.write(reinterpret_cast<char*>(Correction.data()), length*sizeof(double));
		}
	}
	fp2.close();
};

void CountStartPos(string startposfile, map< string,vector<double> >& StartCount, int binsize){
	StartCount.clear();
	ifstream input(startposfile, ios::binary);
	int numtrans;
	input.read((char*)(&numtrans), sizeof(int));
	cout<<numtrans<<endl;
	for(int i=0; i<numtrans; i++){
		int namelen, vectorlen;
		input.read((char*)(&namelen), sizeof(int));
		input.read((char*)(&vectorlen), sizeof(int));

		char* buf=new char[namelen]();
		vector<int> poses(vectorlen, 0);
		vector<double> counts(vectorlen, 0);
		input.read(buf, namelen*sizeof(char));
		input.read((char*)poses.data(), vectorlen*sizeof(int));
		input.read((char*)counts.data(), vectorlen*sizeof(double));
		string name(buf, buf+namelen);

		vector<double> tmp((int)poses.back()/binsize+1, 0);
		for(int j=0; j<poses.size(); j++){
			int idx=poses[j]/binsize;
			tmp[idx]+=counts[j];
		}
		StartCount[name]=tmp;
	}
	input.close();
};

void WriteDeviation(string posfile, string singleposfile, string outfile, map< string,vector<double> >& Corrections, int binsize=50, int escapeedgebins=5){
	map< string,vector<double> > SingleCount;
	CountStartPos(singleposfile, SingleCount, binsize);

	map< string,vector<double> > AllCount;
	CountStartPos(posfile, AllCount, binsize);

	ofstream output(outfile);
	output<<"Name\tbinstart\tbinend\tcov\tcorrectedlen\tnstartstrue\tnsingleendmapped\n";

	for(map< string,vector<double> >::iterator it=Corrections.begin(); it!=Corrections.end(); it++){
		if(it->second.size()<=binsize*escapeedgebins*2)
			continue;
		map< string,vector<double> >::iterator itsingle=SingleCount.find(it->first);
		map< string,vector<double> >::iterator itall=AllCount.find(it->first);
		if(itall==AllCount.end())
			continue;

		double totalreads=0;
		double effectivelen=0;
		double cov=0;
		for(int i=0; i<itall->second.size(); i++)
			totalreads+=itall->second[i];
		for(int i=0; i<it->second.size(); i++)
			effectivelen+=it->second[i];
		cov=totalreads/effectivelen;

		// calculate corrected length for 50bp bin
		vector<double> BinCorrectedLen((int)ceil(1.0*it->second.size()/binsize)-2*escapeedgebins);
		for(int i=0; i<it->second.size(); i++)
			if(i>=escapeedgebins*binsize && i<binsize*(int)ceil(1.0*it->second.size()/binsize)-escapeedgebins*binsize){
				int idx=i/binsize-escapeedgebins;
				BinCorrectedLen[idx]+=it->second[i];
			}
		// write output
		for(int i=0; i<BinCorrectedLen.size(); i++){
			output<<it->first<<"\t"<<(escapeedgebins*binsize+i*binsize)<<"\t"<<((escapeedgebins+i+1)*binsize)<<"\t"<<cov<<"\t";
			if(itsingle!=SingleCount.end() && itsingle->second.size()>i+escapeedgebins)
				output<<BinCorrectedLen[i]<<"\t"<<(itall->second[i+escapeedgebins])<<"\t"<<(itsingle->second[i+escapeedgebins])<<endl;
			else
				output<<BinCorrectedLen[i]<<"\t"<<(itall->second[i+escapeedgebins])<<"\t0\t"<<endl;
		}
	}
	output.close();
};

void WriteRawCorrection(string outfile, map< string,vector<double> >& Corrections){
	ofstream fpout(outfile, ios::out | ios::binary);
	int numtrans=Corrections.size();
	fpout.write(reinterpret_cast<char*>(&numtrans), sizeof(int));
	for(map< string,vector<double> >::iterator it=Corrections.begin(); it!=Corrections.end(); it++){
		int namelen=it->first.size();
		int seqlen=it->second.size();
		string name=it->first;
		fpout.write(reinterpret_cast<char*>(&namelen), sizeof(int));
		fpout.write(reinterpret_cast<char*>(&seqlen), sizeof(int));
		fpout.write(name.c_str(), namelen*sizeof(char));
		fpout.write(reinterpret_cast<char*>(it->second.data()), seqlen*sizeof(double));
	}
	fpout.close();
};

void TestReading(string infile, map< string,vector<double> >& Corrections){
	ifstream fpin(infile, ios::binary);
	int numtrans;
	fpin.read((char*)&numtrans, sizeof(int));
	for(int i=0; i<numtrans; i++){
		int namelen;
		int seqlen;
		fpin.read((char*)&namelen, sizeof(int));
		fpin.read((char*)&seqlen, sizeof(int));

		char* buf=new char[namelen];
		string name;
		vector<double> mydata(seqlen, 0);
		fpin.read(buf, namelen*sizeof(char));
		fpin.read((char*)(mydata.data()), seqlen*sizeof(double));
		name=string(buf);

		vector<double> truedata=Corrections[name];
		double sum=0;
		assert(truedata.size()==mydata.size());
		for(int i=0; i<seqlen; i++)
			sum+=fabs(mydata[i]-truedata[i]);
		cout<<sum<<endl;
	}
};

int main(int argc, char* argv[]){
	if(argc==1){
		printf("readsalmonbias correction <BiasPath> <TransFasta> <SalmonQuant> <OutputFile> (fldLow fldHigh)\n");
		printf("readsalmonbias deviation <BiasPath> <TransFasta> <SalmonQuant> <PosFile> <SinglePosFile> <OutputFile> (fldLow fldHigh)\n");
	}
	else{
		string SeqBiasPath(argv[2]);
		string TransFasta(argv[3]);
		string SalmonQuant(argv[4]);
		string PosFile;
		string SinglePosFile;
		string OutputFile;
		int fldLow, fldHigh;

		if(string(argv[1])=="correction")
			OutputFile=string(argv[5]);
		else{
			PosFile=string(argv[5]);
			SinglePosFile=string(argv[6]);
			OutputFile=string(argv[7]);
		}

		// reading sequence bias
		SeqBiasModel_t seqbias;
		int32_t contextLength, contextLeft, contextRight;
		vector<int32_t> orders, shifts, widths;
		{
			boost::iostreams::filtering_istream fpin;
			fpin.push(boost::iostreams::file_source(SeqBiasPath+"/exp5_seq", std::ios_base::in | std::ios_base::binary));
			ReadSeqBias(fpin, contextLength, seqbias.contextLeft, seqbias.contextRight, seqbias.orders, seqbias.shifts, seqbias.widths, seqbias.Exp5Probs);
		}
		{
			boost::iostreams::filtering_istream fpin;
			fpin.push(boost::iostreams::file_source(SeqBiasPath+"/obs5_seq", std::ios_base::in | std::ios_base::binary));
			ReadSeqBias(fpin, contextLength, contextLeft, contextRight, orders, shifts, widths, seqbias.Obs5Probs);
			assert(contextLeft==seqbias.contextLeft);
			assert(contextRight=seqbias.contextRight);
			assert(orders.size()==seqbias.orders.size());
			for(int i=0; i<orders.size(); i++){
				assert(orders[i]==seqbias.orders[i]);
				assert(shifts[i]==seqbias.shifts[i]);
				assert(widths[i]==seqbias.widths[i]);
			}
		}
		{
			boost::iostreams::filtering_istream fpin;
			fpin.push(boost::iostreams::file_source(SeqBiasPath+"/exp3_seq", std::ios_base::in | std::ios_base::binary));
			ReadSeqBias(fpin, contextLength, contextLeft, contextRight, orders, shifts, widths, seqbias.Exp3Probs);
			assert(contextLeft==seqbias.contextLeft);
			assert(contextRight=seqbias.contextRight);
			assert(orders.size()==seqbias.orders.size());
			for(int i=0; i<orders.size(); i++){
				assert(orders[i]==seqbias.orders[i]);
				assert(shifts[i]==seqbias.shifts[i]);
				assert(widths[i]==seqbias.widths[i]);
			}
		}
		{
			boost::iostreams::filtering_istream fpin;
			fpin.push(boost::iostreams::file_source(SeqBiasPath+"/obs3_seq", std::ios_base::in | std::ios_base::binary));
			ReadSeqBias(fpin, contextLength, contextLeft, contextRight, orders, shifts, widths, seqbias.Obs3Probs);
			assert(contextLeft==seqbias.contextLeft);
			assert(contextRight=seqbias.contextRight);
			assert(orders.size()==seqbias.orders.size());
			for(int i=0; i<orders.size(); i++){
				assert(orders[i]==seqbias.orders[i]);
				assert(shifts[i]==seqbias.shifts[i]);
				assert(widths[i]==seqbias.widths[i]);
			}
		}

		// reading gc bias
		int32_t dtypeObs, dtypeExp;
		Eigen::MatrixXd GCCountsObs, GCCountsExp;
		{
			boost::iostreams::filtering_istream fpin;
			fpin.push(boost::iostreams::file_source(SeqBiasPath+"/exp_gc", std::ios_base::in | std::ios_base::binary));	
			ReadGCBias(fpin, dtypeExp, GCCountsExp);
		}
		{
			boost::iostreams::filtering_istream fpin;
			fpin.push(boost::iostreams::file_source(SeqBiasPath+"/obs_gc", std::ios_base::in | std::ios_base::binary));
			ReadGCBias(fpin, dtypeObs, GCCountsObs);
		}
		GCBiasModel_t gcbias(dtypeObs, GCCountsObs, dtypeExp, GCCountsExp, 1000.0, 3, 1);

		// reading positional bias model
		uint32_t numModels;
		vector<uint32_t> lenBounds;
		vector< vector<double> > obs5models, obs3models, exp5models, exp3models;
		{
			boost::iostreams::filtering_istream fpin;
			fpin.push(boost::iostreams::file_source(SeqBiasPath+"/obs5_pos", std::ios_base::in | std::ios_base::binary));	
			ReadPosBias(fpin, numModels, lenBounds, obs5models);
		}
		{
			boost::iostreams::filtering_istream fpin;
			fpin.push(boost::iostreams::file_source(SeqBiasPath+"/exp5_pos", std::ios_base::in | std::ios_base::binary));	
			ReadPosBias(fpin, numModels, lenBounds, exp5models);
		}
		{
			boost::iostreams::filtering_istream fpin;
			fpin.push(boost::iostreams::file_source(SeqBiasPath+"/obs3_pos", std::ios_base::in | std::ios_base::binary));	
			ReadPosBias(fpin, numModels, lenBounds, obs3models);
		}
		{
			boost::iostreams::filtering_istream fpin;
			fpin.push(boost::iostreams::file_source(SeqBiasPath+"/exp3_pos", std::ios_base::in | std::ios_base::binary));	
			ReadPosBias(fpin, numModels, lenBounds, exp3models);
		}
		cout<<"lenBounds"<<endl;
		for(int32_t i=0; i<lenBounds.size(); i++)
			cout<<lenBounds[i]<<" ";
		cout<<endl;
		cout<<"posBias"<<endl;
		for(int32_t i=0; i<numModels; i++){
			for(int32_t j=0; j<exp3models[i].size(); j++)
				cout<<(obs5models[i][j]/exp5models[i][j]*obs3models[i][j]/exp3models[i][j])<<" ";
			cout<<endl;
		}
		PosBiasModel_t posbias(numModels, lenBounds, obs5models, obs3models, exp5models, exp3models);

		// reading fragment length distribution
		vector<int32_t> RawFLD;
		vector<double> FLD;
		{
			boost::iostreams::filtering_istream fpin;
			fpin.push(boost::iostreams::file_source(SeqBiasPath+"/fld", std::ios_base::in | std::ios_base::binary));
			ReadFLD(fpin, RawFLD);
			FLDKDE(RawFLD, FLD);
		}

		if(string(argv[1])=="deviation" && argc==10){
			fldLow=atoi(argv[8]);
			fldHigh=atoi(argv[9]);
		}
		else if(string(argv[1])=="correction" && argc==8){
			fldLow=atoi(argv[6]);
			fldHigh=atoi(argv[7]);
		}
		else
			GetFLDbound(FLD, fldLow, fldHigh);


		vector<string> TransSequences;
		vector<string> TransNames;
		ReadTransSequence(TransFasta, TransSequences, TransNames);
		cout<<"num transcripts = "<<(TransSequences.size())<<endl;

		map<string, double> SalmonCov;
		ReadSalmonCov(SalmonQuant, SalmonCov);

		map<string, vector<double> > Corrections;
		mutex Corrections_mutex;
		omp_set_num_threads(8);
		#pragma omp parallel for
		for(int i=0; i<TransSequences.size(); i++){
			if(TransSequences[i].size()<seqbias.contextLeft+seqbias.contextRight+1){
				cout<<TransNames[i]<<"\t"<<(TransSequences[i].size())<<endl;
				continue;
			}
			vector<double> Correction=BiasCorrectTrans(TransSequences[i], gcbias, seqbias, posbias, FLD, fldLow, fldHigh);
			// vector<double> Correction=BiasCorrectTrans_wopos(TransSequences[i], gcbias, seqbias, FLD, fldLow, fldHigh);
			lock_guard<std::mutex> guard(Corrections_mutex);
			Corrections[TransNames[i]]=Correction;

		}
		cout<<"num corrected transcripts = "<<(Corrections.size())<<endl;

		if(string(argv[1])=="deviation")
			WriteDeviation(PosFile, SinglePosFile, OutputFile, Corrections, 50, 0);
		else
			WriteRawCorrection(OutputFile, Corrections);
	}
}