#include "CombinedBias.hpp"

using namespace std;


void ReadTranscriptSequence(string filename, vector<string>& TransNames, vector<string>& TransSequences)
{
	// clear result variables
	TransNames.clear();
	TransSequences.clear();
	// read sequence fasta file
	ifstream input(filename);
	string line;
	string tname = "";
	string tseq = "";
	while(getline(input, line)) {
		if (line[0] == '>') {
			if (tname.size() > 0) {
				TransNames.push_back( tname );
				TransSequences.push_back( tseq );
			}
			vector<string> strs;
			boost::split(strs, line, boost::is_any_of("\t "));
			tname = strs[0].substr(1);
			tseq = "";
		}
		else
			tseq += line;
	}
	if (tname.size() > 0) {
		TransNames.push_back( tname );
		TransSequences.push_back( tseq );
	}
	input.close();
};


int32_t main(int32_t argc, char* argv[]) {
	string sequence_file(argv[1]);
	string salmonauxfolder(argv[2]);
	string output_file(argv[3]);
	int32_t is_matrix = stoi(argv[4]);
	string salmonquant = "";

	if (argc > 4)
		salmonquant = (string)argv[5];

	// read transcript sequences
	vector<string> TransNames;
	vector<string> TransSequences;
	ReadTranscriptSequence(sequence_file, TransNames, TransSequences);

	// // using a small set for testing
	// TransNames.resize(100);
	// TransSequences.resize(100);

	// reading bias profiles
	vector<int32_t> RawFLD;
	vector<double> FLD;
	int32_t fldLow, fldHigh;
	ReadFLD(salmonauxfolder+"/fld.gz", RawFLD);
	FLDKDE(RawFLD, FLD);
	GetFLDbound(FLD, fldLow, fldHigh);
	GCBiasModel_t gcbias(salmonauxfolder+"/obs_gc.gz", salmonauxfolder+"/exp_gc.gz", 1000.0, 3, 1);
	SeqBiasModel_t seqbias(salmonauxfolder+"/obs5_seq.gz", salmonauxfolder+"/obs3_seq.gz", salmonauxfolder+"/exp5_seq.gz", salmonauxfolder+"/exp3_seq.gz");
	PosBiasModel_t posbias(salmonauxfolder+"/obs5_pos.gz", salmonauxfolder+"/obs3_pos.gz", salmonauxfolder+"/exp5_pos.gz", salmonauxfolder+"/exp3_pos.gz");

	// perform bias correction
	if (is_matrix == 0) {
		vector< vector<double> > corrections;
		ProcessBias_fragstart(TransSequences, gcbias, seqbias, posbias, FLD, corrections, fldLow, fldHigh, 12);
		if (salmonquant != "")
			AdjustBias(corrections, TransNames, salmonquant);
		WriteCorrection(output_file, corrections, TransNames);
	}
	else {
		vector< vector< tuple<int32_t,int32_t,double> > > corrections;
		ProcessBias_matrix(TransSequences, gcbias, seqbias, posbias, FLD, corrections, fldLow, fldHigh, 12);
		if (salmonquant != "")
			AdjustBias(corrections, TransNames, salmonquant);
		WriteMatrix(output_file, corrections, TransNames);
	}
}