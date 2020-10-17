#include "Transcript.hpp"
#include "GeneGraph.hpp"
#include "CombinedBias.hpp"

using namespace std;


int32_t main(int32_t argc, char* argv[]) {
	string gtffile(argv[1]);
	string genomefasta(argv[2]);
	string salmonauxfolder(argv[3]);
	string salmonquantfile(argv[4]);
	string outprefix(argv[5]);

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
};