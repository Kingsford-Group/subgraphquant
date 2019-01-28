#include "CombinedBias.hpp"


void ReadFLD(string fld_gzip_file, vector<int32_t>& RawFLD, int32_t numbins)
{
	// gzip file handle
	boost::iostreams::filtering_istream fpin;
	fpin.push(boost::iostreams::gzip_decompressor());
	fpin.push(boost::iostreams::file_source(fld_gzip_file, std::ios_base::in | std::ios_base::binary));

	RawFLD.resize(numbins, -1);
	fpin.read((char*)&RawFLD[0], numbins*sizeof(int32_t));
};


void FLDKDE(vector<int32_t>& RawFLD, vector<double>& FLD, int32_t kernel_n, double kernel_p)
{
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


void GetFLDbound(const vector<double>& FLD, int32_t& fldLow, int32_t& fldHigh)
{
	time_t CurrentTime;
	string CurrentTimeStr;

	double quantileCutoffLow = 0.005;
	double FragCount=0;
	for(int32_t i=0; i<FLD.size(); i++)
		FragCount+=FLD[i];
	bool lb=false, ub=false;
	double fldcummulative=0;
	for(int32_t i=0; i<FLD.size(); i++){
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

	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout << "[" << CurrentTimeStr.substr(0, CurrentTimeStr.size()-1) << "] " << "Estimated fragment length lb ="<< fldLow <<" ub = "<< fldHigh <<"\n";
};


vector<double> BiasEffLen_coverage(const string seq, const GCBiasModel_t& gcbias, const SeqBiasModel_t& seqbias, 
	const vector<double>& FLD, int32_t fldLow, int32_t fldHigh)
{
	// process GC condition
	// count the number of G and C around the fragment start (FW) and end (RC) position within contextLeft and contextRight
	// use GCWindowFW and GCWindowRC to indicate the number if counted among how many nucleotides (window size)
	vector<int32_t> GCCondFW(seq.size(), 0);
	vector<int32_t> GCCondRC(seq.size(), 0);
	vector<int32_t> GCWindowFW(seq.size(), 0);
	vector<int32_t> GCWindowRC(seq.size(), 0);
	for(int32_t i = 0; i < seq.size(); i++){
		int32_t gcFWcount = 0, gcRCcount = 0;
		for(int32_t j = i-gcbias.contextLeft; j < i+gcbias.contextRight+1; j++){
			if(j >= 0 && j < seq.size())
				if(seq[j]=='G' || seq[j]=='g' || seq[j]=='C' || seq[j]=='c')
					gcFWcount++;
		}
		GCCondFW[i] = gcFWcount;
		GCWindowFW[i] = ( min(i+gcbias.contextRight+1, (int)seq.size()) - max(0, i-gcbias.contextLeft) );
		for(int32_t j = i+gcbias.contextLeft; j > i-gcbias.contextRight-1; j--){
			if(j >= 0 && j < seq.size())
				if(seq[j]=='G' || seq[j]=='g' || seq[j]=='C' || seq[j]=='c')
					gcRCcount++;
		}
		GCCondRC[i] = gcRCcount;
		GCWindowRC[i] = ( min((int)seq.size(), i+gcbias.contextLeft+1) - max(0, i-gcbias.contextRight) );
	}
	// process GC content vector
	// this is for calculating GC percentage of the entire fragment
	vector<int32_t> GCRawCount(seq.size()+1, 0);
	int32_t cummulative=0;
	for(int32_t i = 0; i < seq.size(); i++){
		if(seq[i]=='G' || seq[i]=='g' || seq[i]=='C' || seq[i]=='c'){
			cummulative++;
		}
		GCRawCount[i+1]=cummulative;
	}

	// SEQbias
	// process 5' seqbias positional ratio
	vector<double> Seq5Ratio(seq.size(), 1);
	vector<double> Seq3Ratio(seq.size(), 1);
	Mer mer;
	mer.k(seqbias.contextLeft+seqbias.contextRight+1);
	mer.from_chars(seq.c_str());
	for(int32_t i = seqbias.contextLeft; i < seq.size()-seqbias.contextRight; i++){
		double obsvalue=0, expvalue=0;
		for(int32_t j = 0; j < seqbias.contextLeft+seqbias.contextRight+1; j++){
			int32_t idx = mer.get_bits(seqbias.shifts[j], seqbias.widths[j]);
			obsvalue += seqbias.Obs5Probs(idx, j);
			expvalue += seqbias.Exp5Probs(idx, j);
		}
		Seq5Ratio[i] = exp(obsvalue-expvalue);
		mer.shift_left(seq[i+seqbias.contextRight+1]);
	}
	// process 3' seqbias positional ratio
	string rcseq = ReverseComplement(seq.cbegin(), seq.cend());
	mer.from_chars(rcseq.c_str());
	for(int32_t i = seqbias.contextLeft; i < rcseq.size()-seqbias.contextRight; i++){
		double obsvalue=0, expvalue=0;
		for(int32_t j = 0; j < seqbias.contextLeft+seqbias.contextRight+1; j++){
			int32_t idx=mer.get_bits(seqbias.shifts[j], seqbias.widths[j]);
			obsvalue+=seqbias.Obs3Probs(idx, j);
			expvalue+=seqbias.Exp3Probs(idx, j);
		}
		Seq3Ratio[i]=exp(obsvalue-expvalue);
		mer.shift_left(rcseq[i+seqbias.contextRight+1]);
	}
	reverse(Seq3Ratio.begin(), Seq3Ratio.end());

	// effective length calculated based on coverage probability
	// initialize correction
	vector<double> correction(seq.size(), 0);
	// calculate total number of possible fragments
	double FragCount=0;
	int32_t minFLDpossible=(seq.size()<FLD.size())?1:fldLow;
	int32_t maxFLDpossible=(seq.size()<FLD.size())?seq.size():fldHigh;
	for(int32_t i = 0; i < maxFLDpossible; i++)
		FragCount += FLD[i];
	// calculate the multiplier of generating a k length fragment at position i
	for(int32_t k = minFLDpossible; k < maxFLDpossible; k++){
		for(int32_t i = 0; i < seq.size()-k-1; i++){
			// seqBias
			double seq5factor=Seq5Ratio[i];
			double seq3factor=Seq3Ratio[i+k-1];
			// gcBias
			double gcfactor;
			double gccondition=1.0*(GCCondFW[i]+GCCondRC[i+k-1])/(GCWindowFW[i]+GCWindowRC[i+k-1]);
			int32_t gccondbin=(int32_t)(gccondition*gcbias.Condbin);
			if(gccondbin==gcbias.Condbin)
				gccondbin--;
			double gcfraction=1.0*(GCRawCount[i+k]-GCRawCount[i])/k;
			int32_t gcbin=(int32_t)(gcbias.GCbin*lrint(gcfraction*100)/100.0);
			if(gcbin==gcbias.GCbin)
				gcbin--;
			gcfactor=gcbias.GCBiasRatio(gccondbin, gcbin);

			double m = 1.0 * (FLD[k]) * seq5factor * seq3factor * gcfactor / FragCount / k;
			for (int32_t l = i; l < i + k; l++)
				correction[l] += m;
		}
	}

	return correction;
};


vector<double> BiasEffLen_fragstart(const string seq, const GCBiasModel_t& gcbias, const SeqBiasModel_t& seqbias, const PosBiasModel_t& posbias, 
	const vector<double>& FLD, int32_t fldLow, int32_t fldHigh)
{
	// process GC condition
	// count the number of G and C around the fragment start (FW) and end (RC) position within contextLeft and contextRight
	// use GCWindowFW and GCWindowRC to indicate the number if counted among how many nucleotides (window size)
	vector<int32_t> GCCondFW(seq.size(), 0);
	vector<int32_t> GCCondRC(seq.size(), 0);
	vector<int32_t> GCWindowFW(seq.size(), 0);
	vector<int32_t> GCWindowRC(seq.size(), 0);
	for(int32_t i = 0; i < seq.size(); i++){
		int32_t gcFWcount = 0, gcRCcount = 0;
		for(int32_t j = i-gcbias.contextLeft; j < i+gcbias.contextRight+1; j++){
			if(j >= 0 && j < seq.size())
				if(seq[j]=='G' || seq[j]=='g' || seq[j]=='C' || seq[j]=='c')
					gcFWcount++;
		}
		GCCondFW[i] = gcFWcount;
		GCWindowFW[i] = ( min(i+gcbias.contextRight+1, (int)seq.size()) - max(0, i-gcbias.contextLeft) );
		for(int32_t j = i+gcbias.contextLeft; j > i-gcbias.contextRight-1; j--){
			if(j >= 0 && j < seq.size())
				if(seq[j]=='G' || seq[j]=='g' || seq[j]=='C' || seq[j]=='c')
					gcRCcount++;
		}
		GCCondRC[i] = gcRCcount;
		GCWindowRC[i] = ( min((int)seq.size(), i+gcbias.contextLeft+1) - max(0, i-gcbias.contextRight) );
	}
	// process GC content vector
	// this is for calculating GC percentage of the entire fragment
	vector<int32_t> GCRawCount(seq.size()+1, 0);
	int32_t cummulative=0;
	for(int32_t i = 0; i < seq.size(); i++){
		if(seq[i]=='G' || seq[i]=='g' || seq[i]=='C' || seq[i]=='c'){
			cummulative++;
		}
		GCRawCount[i+1]=cummulative;
	}

	// SEQbias
	// process 5' seqbias positional ratio
	vector<double> Seq5Ratio(seq.size(), 1);
	vector<double> Seq3Ratio(seq.size(), 1);
	Mer mer;
	mer.k(seqbias.contextLeft+seqbias.contextRight+1);
	mer.from_chars(seq.c_str());
	for(int32_t i = seqbias.contextLeft; i < seq.size()-seqbias.contextRight; i++){
		double obsvalue=0, expvalue=0;
		for(int32_t j = 0; j < seqbias.contextLeft+seqbias.contextRight+1; j++){
			int32_t idx = mer.get_bits(seqbias.shifts[j], seqbias.widths[j]);
			obsvalue += seqbias.Obs5Probs(idx, j);
			expvalue += seqbias.Exp5Probs(idx, j);
		}
		Seq5Ratio[i] = exp(obsvalue-expvalue);
		mer.shift_left(seq[i+seqbias.contextRight+1]);
	}
	// process 3' seqbias positional ratio
	string rcseq = ReverseComplement(seq.cbegin(), seq.cend());
	mer.from_chars(rcseq.c_str());
	for(int32_t i = seqbias.contextLeft; i < rcseq.size()-seqbias.contextRight; i++){
		double obsvalue=0, expvalue=0;
		for(int32_t j = 0; j < seqbias.contextLeft+seqbias.contextRight+1; j++){
			int32_t idx=mer.get_bits(seqbias.shifts[j], seqbias.widths[j]);
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
	int32_t li = posbias.lenClass((int32_t)seq.size());
	for(int32_t i = 0; i < seq.size()-(seqbias.contextLeft+seqbias.contextRight+1); i++){
		double fracP = 1.0*i/seq.size();
		double obs5 = max(0.001, posbias.obs5Splines[li](fracP));
		double obs3 = max(0.001, posbias.obs3Splines[li](fracP));
		double exp5 = max(0.001, posbias.exp5Splines[li](fracP));
		double exp3 = max(0.001, posbias.exp3Splines[li](fracP));
		Pos5Ratio[i] = obs5/exp5;
		Pos3Ratio[i] = obs3/exp3;
	}

	// effective length calculated based on coverage probability
	// initialize correction
	vector<double> correction(seq.size(), 0);
	// calculate total number of possible fragments
	double FragCount=0;
	int32_t minFLDpossible=(seq.size()<FLD.size())?1:fldLow;
	int32_t maxFLDpossible=(seq.size()<FLD.size())?seq.size():fldHigh;
	for(int32_t i = 0; i < maxFLDpossible; i++)
		FragCount += FLD[i];
	// calculate the multiplier of generating a k length fragment at position i
	for(int32_t k = minFLDpossible; k < maxFLDpossible; k++){
		for(int32_t i = 0; i < seq.size()-k-1; i++){
			// seqBias
			double seq5factor = Seq5Ratio[i];
			double seq3factor = Seq3Ratio[i+k-1];
			// gcBias
			double gcfactor;
			double gccondition = 1.0*(GCCondFW[i]+GCCondRC[i+k-1])/(GCWindowFW[i]+GCWindowRC[i+k-1]);
			int32_t gccondbin = (int32_t)(gccondition*gcbias.Condbin);
			if(gccondbin == gcbias.Condbin)
				gccondbin--;
			double gcfraction = 1.0*(GCRawCount[i+k]-GCRawCount[i])/k;
			int32_t gcbin = (int32_t)(gcbias.GCbin*lrint(gcfraction*100)/100.0);
			if(gcbin == gcbias.GCbin)
				gcbin--;
			gcfactor = gcbias.GCBiasRatio(gccondbin, gcbin);
			// posBias
			double posfactor = Pos5Ratio[i]*Pos3Ratio[i+k-1];

			correction[i] += 1.0 * (FLD[k]) * seq5factor * seq3factor * gcfactor * posfactor / FragCount;
		}
	}

	return correction;
};


vector< tuple<int32_t,int32_t,double> > BiasEffLen_matrix(const string seq, const GCBiasModel_t& gcbias, const SeqBiasModel_t& seqbias, const PosBiasModel_t& posbias, 
	const vector<double>& FLD, int32_t fldLow, int32_t fldHigh)
{
	// process GC condition
	// count the number of G and C around the fragment start (FW) and end (RC) position within contextLeft and contextRight
	// use GCWindowFW and GCWindowRC to indicate the number if counted among how many nucleotides (window size)
	vector<int32_t> GCCondFW(seq.size(), 0);
	vector<int32_t> GCCondRC(seq.size(), 0);
	vector<int32_t> GCWindowFW(seq.size(), 0);
	vector<int32_t> GCWindowRC(seq.size(), 0);
	for(int32_t i = 0; i < seq.size(); i++){
		int32_t gcFWcount = 0, gcRCcount = 0;
		for(int32_t j = i-gcbias.contextLeft; j < i+gcbias.contextRight+1; j++){
			if(j >= 0 && j < seq.size())
				if(seq[j]=='G' || seq[j]=='g' || seq[j]=='C' || seq[j]=='c')
					gcFWcount++;
		}
		GCCondFW[i] = gcFWcount;
		GCWindowFW[i] = ( min(i+gcbias.contextRight+1, (int)seq.size()) - max(0, i-gcbias.contextLeft) );
		for(int32_t j = i+gcbias.contextLeft; j > i-gcbias.contextRight-1; j--){
			if(j >= 0 && j < seq.size())
				if(seq[j]=='G' || seq[j]=='g' || seq[j]=='C' || seq[j]=='c')
					gcRCcount++;
		}
		GCCondRC[i] = gcRCcount;
		GCWindowRC[i] = ( min((int)seq.size(), i+gcbias.contextLeft+1) - max(0, i-gcbias.contextRight) );
	}
	// process GC content vector
	// this is for calculating GC percentage of the entire fragment
	vector<int32_t> GCRawCount(seq.size()+1, 0);
	int32_t cummulative=0;
	for(int32_t i = 0; i < seq.size(); i++){
		if(seq[i]=='G' || seq[i]=='g' || seq[i]=='C' || seq[i]=='c'){
			cummulative++;
		}
		GCRawCount[i+1]=cummulative;
	}

	// SEQbias
	// process 5' seqbias positional ratio
	vector<double> Seq5Ratio(seq.size(), 1);
	vector<double> Seq3Ratio(seq.size(), 1);
	Mer mer;
	mer.k(seqbias.contextLeft+seqbias.contextRight+1);
	mer.from_chars(seq.c_str());
	for(int32_t i = seqbias.contextLeft; i < seq.size()-seqbias.contextRight; i++){
		double obsvalue=0, expvalue=0;
		for(int32_t j = 0; j < seqbias.contextLeft+seqbias.contextRight+1; j++){
			int32_t idx = mer.get_bits(seqbias.shifts[j], seqbias.widths[j]);
			obsvalue += seqbias.Obs5Probs(idx, j);
			expvalue += seqbias.Exp5Probs(idx, j);
		}
		Seq5Ratio[i] = exp(obsvalue-expvalue);
		mer.shift_left(seq[i+seqbias.contextRight+1]);
	}
	// process 3' seqbias positional ratio
	string rcseq = ReverseComplement(seq.cbegin(), seq.cend());
	mer.from_chars(rcseq.c_str());
	for(int32_t i = seqbias.contextLeft; i < rcseq.size()-seqbias.contextRight; i++){
		double obsvalue=0, expvalue=0;
		for(int32_t j = 0; j < seqbias.contextLeft+seqbias.contextRight+1; j++){
			int32_t idx=mer.get_bits(seqbias.shifts[j], seqbias.widths[j]);
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
	int32_t li = posbias.lenClass((int32_t)seq.size());
	for(int32_t i = 0; i < seq.size()-(seqbias.contextLeft+seqbias.contextRight+1); i++){
		double fracP = 1.0*i/seq.size();
		double obs5 = max(0.001, posbias.obs5Splines[li](fracP));
		double obs3 = max(0.001, posbias.obs3Splines[li](fracP));
		double exp5 = max(0.001, posbias.exp5Splines[li](fracP));
		double exp3 = max(0.001, posbias.exp3Splines[li](fracP));
		Pos5Ratio[i] = obs5/exp5;
		Pos3Ratio[i] = obs3/exp3;
	}

	// effective length calculated based on coverage probability
	vector< tuple<int32_t,int32_t,double> > correction;
	// calculate total number of possible fragments
	double FragCount=0;
	int32_t minFLDpossible=(seq.size()<FLD.size())?1:fldLow;
	int32_t maxFLDpossible=(seq.size()<FLD.size())?seq.size():fldHigh;
	for(int32_t i = 0; i < maxFLDpossible; i++)
		FragCount += FLD[i];
	// calculate the multiplier of generating a k length fragment at position i
	for(int32_t k = minFLDpossible; k < maxFLDpossible; k++){
		for(int32_t i = 0; i < seq.size()-k-1; i++){
			// seqBias
			double seq5factor = Seq5Ratio[i];
			double seq3factor = Seq3Ratio[i+k-1];
			// gcBias
			double gcfactor;
			double gccondition = 1.0*(GCCondFW[i]+GCCondRC[i+k-1])/(GCWindowFW[i]+GCWindowRC[i+k-1]);
			int32_t gccondbin = (int32_t)(gccondition*gcbias.Condbin);
			if(gccondbin == gcbias.Condbin)
				gccondbin--;
			double gcfraction = 1.0*(GCRawCount[i+k]-GCRawCount[i])/k;
			int32_t gcbin = (int32_t)(gcbias.GCbin*lrint(gcfraction*100)/100.0);
			if(gcbin == gcbias.GCbin)
				gcbin--;
			gcfactor = gcbias.GCBiasRatio(gccondbin, gcbin);
			// posBias
			double posfactor = Pos5Ratio[i]*Pos3Ratio[i+k-1];
			if ((FLD[k]) > 0 && seq5factor > 0 && seq3factor > 0 && gcfactor > 0 && posfactor > 0)
				correction.push_back( make_tuple(i, i+k, 1.0 * (FLD[k]) * seq5factor * seq3factor * gcfactor * posfactor / FragCount ) );
		}
	}

	return correction;
};


void ProcessBias_coverage(const vector<string>& sequences, const GCBiasModel_t& gcbias, const SeqBiasModel_t& seqbias, const vector<double>& FLD,
	vector< vector<double> >& corrections, int32_t fldLow, int32_t fldHigh, int32_t threads)
{
	time_t CurrentTime;
	string CurrentTimeStr;
	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout<<"["<<CurrentTimeStr.substr(0, CurrentTimeStr.size()-1)<<"] "<<"Calculating bias-corrected length."<<endl;

	// clear variable
	corrections.clear();
	for (int32_t i = 0; i < sequences.size(); i++) {
		vector<double> tmp(sequences[i].size(), 1);
		corrections.push_back(tmp);
	}
	corrections.shrink_to_fit();

	// calculate bias corrected length per position
	mutex corrections_mutex;
	omp_set_num_threads(threads);
	#pragma omp parallel for
	for(int32_t i = 0; i < sequences.size(); i++){
		if(sequences[i].size() < seqbias.contextLeft+seqbias.contextRight+1){
			cout <<"Not calculating effective length due to too-short length: "<< sequences[i] <<"\t"<< (sequences[i].size()) << endl;
			continue;
		}
		vector<double> tmp_correction = BiasEffLen_coverage(sequences[i], gcbias, seqbias, FLD, fldLow, fldHigh);
		lock_guard<std::mutex> guard(corrections_mutex);
		corrections[i] = tmp_correction;
	}

	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout << "[" << CurrentTimeStr.substr(0, CurrentTimeStr.size()-1) << "] " << "Finish.\n";
};


void ProcessBias_fragstart(const vector<string>& sequences, const GCBiasModel_t& gcbias, const SeqBiasModel_t& seqbias, const PosBiasModel_t& posbias, 
	const vector<double>& FLD, vector< vector<double> >& corrections, int32_t fldLow, int32_t fldHigh, int32_t threads)
{
	time_t CurrentTime;
	string CurrentTimeStr;
	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout<<"["<<CurrentTimeStr.substr(0, CurrentTimeStr.size()-1)<<"] "<<"Calculating bias-corrected length."<<endl;

	// clear variable
	corrections.clear();
	for (int32_t i = 0; i < sequences.size(); i++) {
		vector<double> tmp(sequences[i].size(), 1);
		corrections.push_back(tmp);
	}
	corrections.shrink_to_fit();

	// calculate bias corrected length per position
	mutex corrections_mutex;
	omp_set_num_threads(threads);
	#pragma omp parallel for
	for(int32_t i = 0; i < sequences.size(); i++){
		if(sequences[i].size() < seqbias.contextLeft+seqbias.contextRight+1){
			cout <<"Not calculating effective length due to too-short length: "<< sequences[i] <<"\t"<< (sequences[i].size()) << endl;
			continue;
		}
		vector<double> tmp_correction = BiasEffLen_fragstart(sequences[i], gcbias, seqbias, posbias, FLD, fldLow, fldHigh);
		lock_guard<std::mutex> guard(corrections_mutex);
		corrections[i] = tmp_correction;
	}

	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout << "[" << CurrentTimeStr.substr(0, CurrentTimeStr.size()-1) << "] " << "Finish.\n";
};


/*void ProcessBias_matrix(const vector<string>& sequences, const GCBiasModel_t& gcbias, const SeqBiasModel_t& seqbias, const PosBiasModel_t& posbias, 
	const vector<double>& FLD, vector< Eigen::MatrixXd >& corrections, int32_t fldLow, int32_t fldHigh, int32_t threads)
{
	time_t CurrentTime;
	string CurrentTimeStr;
	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout<<"["<<CurrentTimeStr.substr(0, CurrentTimeStr.size()-1)<<"] "<<"Calculating bias-corrected length."<<endl;

	// clear variable
	corrections.clear();
	for (int32_t i = 0; i < sequences.size(); i++) {
		cout << i <<"\t"<< sequences[i].size() << endl;
		Eigen::MatrixXd tmp = Eigen::MatrixXd::Zero(sequences[i].size(), sequences[i].size());
		corrections.push_back(tmp);
	}
	corrections.shrink_to_fit();

	// calculate bias corrected length per position
	mutex corrections_mutex;
	omp_set_num_threads(threads);
	#pragma omp parallel for
	for(int32_t i = 0; i < sequences.size(); i++){
		if(sequences[i].size() < seqbias.contextLeft+seqbias.contextRight+1){
			cout <<"Not calculating effective length due to too-short length: "<< sequences[i] <<"\t"<< (sequences[i].size()) << endl;
			continue;
		}
		Eigen::MatrixXd tmp_correction = BiasEffLen_matrix(sequences[i], gcbias, seqbias, posbias, FLD, fldLow, fldHigh);
		lock_guard<std::mutex> guard(corrections_mutex);
		corrections[i] = tmp_correction;
	}

	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout << "[" << CurrentTimeStr.substr(0, CurrentTimeStr.size()-1) << "] " << "Finish.\n";
};*/


void ProcessBias_matrix(const vector<string>& sequences, const GCBiasModel_t& gcbias, const SeqBiasModel_t& seqbias, const PosBiasModel_t& posbias, 
	const vector<double>& FLD, vector< vector< tuple<int32_t,int32_t,double> > >& corrections, int32_t fldLow, int32_t fldHigh, int32_t threads)
{
	time_t CurrentTime;
	string CurrentTimeStr;
	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout<<"["<<CurrentTimeStr.substr(0, CurrentTimeStr.size()-1)<<"] "<<"Calculating bias-corrected length."<<endl;

	// clear variable
	corrections.clear();
	for (int32_t i = 0; i < sequences.size(); i++) {
		vector< tuple<int32_t,int32_t,double> > tmp;
		corrections.push_back(tmp);
	}
	corrections.shrink_to_fit();

	// calculate bias corrected length per position
	mutex corrections_mutex;
	omp_set_num_threads(threads);
	#pragma omp parallel for
	for(int32_t i = 0; i < sequences.size(); i++){
		if(sequences[i].size() < seqbias.contextLeft+seqbias.contextRight+1){
			cout <<"Not calculating effective length due to too-short length: "<< sequences[i] <<"\t"<< (sequences[i].size()) << endl;
			continue;
		}
		vector< tuple<int32_t,int32_t,double> > tmp_correction = BiasEffLen_matrix(sequences[i], gcbias, seqbias, posbias, FLD, fldLow, fldHigh);
		lock_guard<std::mutex> guard(corrections_mutex);
		corrections[i] = tmp_correction;
	}

	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout << "[" << CurrentTimeStr.substr(0, CurrentTimeStr.size()-1) << "] " << "Finish.\n";
};


void AdjustBias(vector< vector<double> >& corrections, const vector<string>& transnames, string salmonquantfile)
{
	time_t CurrentTime;
	string CurrentTimeStr;
	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout<<"["<<CurrentTimeStr.substr(0, CurrentTimeStr.size()-1)<<"] "<<"Adjusting bias-corrected based on salmon output."<<endl;

	// read the corrected length from salmon quant file
	map<string,double> EffLen_salmon;
	ifstream input(salmonquantfile);
	string line;
	int32_t linecount = 0;
	while(getline(input, line)) {
		linecount ++;
		if (linecount == 1)
			continue;
		vector<string> strs;
		boost::split(strs, line, boost::is_any_of("\t"));
		EffLen_salmon[strs[0]] = stod(strs[2]);
	}
	input.close();
	// multiply the current corrections to make each transcript sum to the effective length from salmon
	int32_t num_skipped = 0;
	assert(corrections.size() == transnames.size());
	for (int32_t i = 0; i < corrections.size(); i++) {
		// retrieve tname, bias, and salmon efflen
		string tname = transnames[i];
		vector<double>& bias = corrections[i];
		map<string,double>::const_iterator itmap = EffLen_salmon.find(tname);
		if (itmap == EffLen_salmon.cend()) {
			num_skipped ++;
			continue;
		}
		double efflen = itmap->second;
		// calculate current sum
		double s = 0;
		for (int32_t j = 0; j < bias.size(); j++)
			s += bias[j];
		for (int32_t j = 0; j < bias.size(); j++)
			bias[j] *= efflen / s;
		// sanity check: the changed bias has been updated
		double new_s = 0;
		for (int32_t j = 0; j < bias.size(); j++)
			new_s += bias[j];
		if (fabs(new_s - efflen) > 1e-3)
			cout << i <<"\t"<< new_s <<"\t"<< (std::accumulate(bias.begin(), bias.end(), 0)) << endl;
		assert(fabs(new_s - efflen) < 1e-3);
	}

	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout << "[" << CurrentTimeStr.substr(0, CurrentTimeStr.size()-1) << "] " << "Finish. " << num_skipped << " are skipped.\n";
};


void AdjustBias(vector< vector< tuple<int32_t,int32_t,double> > >& corrections, const vector<string>& transnames, string salmonquantfile)
{
	time_t CurrentTime;
	string CurrentTimeStr;
	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout<<"["<<CurrentTimeStr.substr(0, CurrentTimeStr.size()-1)<<"] "<<"Adjusting bias-corrected based on salmon output."<<endl;

	// read the corrected length from salmon quant file
	map<string,double> EffLen_salmon;
	ifstream input(salmonquantfile);
	string line;
	int32_t linecount = 0;
	while(getline(input, line)) {
		linecount ++;
		if (linecount == 1)
			continue;
		vector<string> strs;
		boost::split(strs, line, boost::is_any_of("\t"));
		EffLen_salmon[strs[0]] = stod(strs[2]);
	}
	input.close();
	// multiply the current corrections to make each transcript sum to the effective length from salmon
	int32_t num_skipped = 0;
	assert(corrections.size() == transnames.size());
	for (int32_t i = 0; i < corrections.size(); i++) {
		// retrieve tname, bias, and salmon efflen
		string tname = transnames[i];
		vector< tuple<int32_t,int32_t,double> >& bias = corrections[i];
		map<string,double>::const_iterator itmap = EffLen_salmon.find(tname);
		if (itmap == EffLen_salmon.cend()) {
			num_skipped ++;
			continue;
		}
		double efflen = itmap->second;
		// calculate current sum
		double s = 0;
		for (int32_t j = 0; j < bias.size(); j++)
			s += get<2>(bias[j]);
		for (int32_t j = 0; j < bias.size(); j++)
			get<2>(bias[j]) *= efflen / s;
		// sanity check: the changed bias has been updated
		double new_s = 0;
		for (int32_t j = 0; j < corrections[i].size(); j++)
			new_s += get<2>(corrections[i][j]);
		if (fabs(new_s - efflen) > 1e-3)
			cout << i <<"\t"<< new_s << endl;
		assert(fabs(new_s - efflen) < 1e-3);
	}

	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout << "[" << CurrentTimeStr.substr(0, CurrentTimeStr.size()-1) << "] " << "Finish. " << num_skipped << " are skipped.\n";
};
