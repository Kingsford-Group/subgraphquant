#include "POSbias.hpp"


void ReadPosBiasFile(string gzip_file, uint32_t& numModels, vector<uint32_t>& lenBounds, vector< vector<double> >& masses)
{
	// file handle
	boost::iostreams::filtering_istream fpin;
	fpin.push(boost::iostreams::gzip_decompressor());
	fpin.push(boost::iostreams::file_source(gzip_file, std::ios_base::in | std::ios_base::binary));

	// read file content
	fpin.read((char*)(&numModels), sizeof(uint32_t));

	lenBounds.resize(numModels);
	for(int32_t i=0; i<numModels; i++)
		fpin.read((char*)(&lenBounds[i]), sizeof(uint32_t));

	masses.clear();
	for(int32_t i=0; i<numModels; i++){
		uint32_t modelLen=0;
		fpin.read((char*)(&modelLen), sizeof(uint32_t));
		vector<double> tmpmodel(modelLen, 0);
		fpin.read((char*)tmpmodel.data(), modelLen*sizeof(double));
		masses.push_back(tmpmodel);
	}
};


PosBiasModel_t::PosBiasModel_t(string obs5_gzip_file, string obs3_gzip_file, string exp5_gzip_file, string exp3_gzip_file)
{
	time_t CurrentTime;
	string CurrentTimeStr;
	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout<<"["<<CurrentTimeStr.substr(0, CurrentTimeStr.size()-1)<<"] "<<"Loading positional bias model."<<endl;

	// reading obs5
	ReadPosBiasFile(obs5_gzip_file, numModels, lenBounds, obs5models);

	// reading obs3
	uint32_t tmp_numModels;
	vector<uint32_t> tmp_lenBounds;
	ReadPosBiasFile(obs3_gzip_file, tmp_numModels, tmp_lenBounds, obs3models);
	assert(tmp_numModels == numModels);
	for (int32_t i = 0; i < tmp_lenBounds.size(); i++)
		assert(tmp_lenBounds[i] == lenBounds[i]);

	// reading exp5
	ReadPosBiasFile(exp5_gzip_file, tmp_numModels, tmp_lenBounds, exp5models);
	assert(tmp_numModels == numModels);
	for (int32_t i = 0; i < tmp_lenBounds.size(); i++)
		assert(tmp_lenBounds[i] == lenBounds[i]);

	// reading exp3
	ReadPosBiasFile(exp3_gzip_file, tmp_numModels, tmp_lenBounds, exp3models);
	assert(tmp_numModels == numModels);
	for (int32_t i = 0; i < tmp_lenBounds.size(); i++)
		assert(tmp_lenBounds[i] == lenBounds[i]);

	// interpolation between sampled points
	ProcessSplines();

	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout << "[" << CurrentTimeStr.substr(0, CurrentTimeStr.size()-1) << "] " << "Finish.\n";
};


void PosBiasModel_t::ProcessSplines(){
		for(int32_t i=0; i<numModels; i++){
			vector<double> splineMass(obs5models[i].size()+2);
			vector<double> splineBins(obs5models[i].size()+2);

			// observed 5 model
			splineBins[0]=0;
			splineMass[0]=obs5models[i][0];
			for(int32_t j=0; j<obs5models[i].size(); j++){
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
			for(int32_t j=0; j<obs3models[i].size(); j++){
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
			for(int32_t j=0; j<exp5models[i].size(); j++){
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
			for(int32_t j=0; j<exp3models[i].size(); j++){
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