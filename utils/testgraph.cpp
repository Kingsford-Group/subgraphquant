#include "Transcript.hpp"
#include "GeneGraph.hpp"
#include "CombinedBias.hpp"

using namespace std;


int32_t main(int32_t argc, char* argv[]) {
	int32_t mode = atoi(argv[1]); // mode 0: without bias correction; 1: with bias correction
	string gtffile(argv[2]);
	string genomefasta(argv[3]);
	string salmonauxfolder(argv[4]);
	string salmonquantfile(argv[5]);
	string outprefix(argv[6]);

	vector<string> RefName;
	std::map<string,uint8_t> RefTable;
	BuildRefName(gtffile, RefName, RefTable);

	// read GTF and retrieve transcript sequence
	vector<Transcript_t> Transcripts;
	map< string, vector<int32_t> > GeneTransMap;
	map< string,int32_t > TransIndex;
	vector<string> sequences;
	ReadGTF (gtffile, Transcripts, RefTable);

	// // use subset of genes for testing
	// vector<string> caredgenes = {"ENSG00000197601.12", "ENSG00000084070.11", "ENSG00000242247.10", "ENSG00000166997.7", "ENSG00000171055.14"};
	// sort(caredgenes.begin(), caredgenes.end());
	// vector<Transcript_t> tmp_Transcripts;
	// for (Transcript_t& t : Transcripts) {
	// 	if (binary_search(caredgenes.begin(), caredgenes.end(), t.GeneID))
	// 		tmp_Transcripts.push_back(t);
	// }
	// Transcripts = tmp_Transcripts;

	Map_Gene_Trans (Transcripts, GeneTransMap);
	MapIndex(Transcripts, TransIndex);
	ReadTransSequence_genome(genomefasta, RefTable, Transcripts, sequences);

	vector<GeneGraph_t> GeneGraphs;
	ConstructAllGeneGraphs(Transcripts, GeneTransMap, GeneGraphs);
	WriteGraph(outprefix+"_graph_fragstart_beforebias.txt", GeneGraphs, RefName, Transcripts);

	if (mode == 1) {
		// read salmon quantification TPM
		vector<double> expr;
		ReadSalmonTPM(salmonquantfile, TransIndex, expr);

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
		vector< vector<double> > corrections;
		ProcessBias_fragstart(sequences, gcbias, seqbias, posbias, FLD, corrections, fldLow, fldHigh, 12);
		vector<string> transnames;
		for (const Transcript_t& t : Transcripts)
			transnames.push_back(t.TransID);
		AdjustBias(corrections, transnames, salmonquantfile);
		WriteCorrection(outprefix+"_corrections_fragstart.dat", corrections, transnames);

		UpdateAllBiasMultiplier(Transcripts, corrections, expr, GeneTransMap, GeneGraphs);

		WriteGraph(outprefix+"_graph_fragstart.txt", GeneGraphs, RefName, Transcripts);
	}
};