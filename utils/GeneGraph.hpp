#ifndef __GeneGraph_T__
#define __GeneGraph_T__

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <string>
#include "boost/algorithm/string.hpp"

using namespace std;
class Transcript_t;

// both segments and junctions are turned into Node_t
// for segments, Length is its original length, and BiasMultiplier is the original bias model
// for junctions, Length is 1, and BiasMultiplier is the average of the 2 adjacenct segments
class Node_t{
public:
	int32_t ID;
	uint8_t Chr;
	int32_t StartPos;
	int32_t EndPos;
	bool Strand;
	bool IsSegment;
	vector<int32_t> InEdges;
	vector<int32_t> OutEdges;

public:
	int32_t Length;
	// effective length of the node considering sequence biases and GC biases
	double BiasMultiplier;
	// the actual coverage of assigned reads. For each read, its contribution to the coverage is calculated by covered length / node length
	double Coverage;

public:
	Node_t(){};
	Node_t(uint8_t Chr, int32_t StartPos, int32_t EndPos, bool Strand, bool IsSegment): Chr(Chr), StartPos(StartPos), EndPos(EndPos), Strand(Strand), IsSegment(IsSegment) {BiasMultiplier = 1;};
	Node_t(int32_t ID, uint8_t Chr, int32_t StartPos, int32_t EndPos, bool Strand, bool IsSegment): ID(ID), Chr(Chr), StartPos(StartPos), EndPos(EndPos), Strand(Strand), IsSegment(IsSegment) {BiasMultiplier = 2;};

	bool operator < (const Node_t& rhs) const {
		if (Chr != rhs.Chr)
			return Chr < rhs.Chr;
		else if (StartPos != rhs.StartPos)
			return StartPos < rhs.StartPos;
		else
			return EndPos < rhs.EndPos;
	};

	void setID (int32_t id) {
		ID = id;
	};

	void CalculateLength() {
		Length = (IsSegment) ? (EndPos-StartPos) : 1;
		assert(Length > 0);
	};
};


// edge links original segments and junctions.
// it doesn't carry information of coverage or flow or length.
class Edge_t{
public:
	int32_t ID;
	int32_t Node_1;
	int32_t Node_2;
	vector<string> IncidentTranscripts;
	// estimated flow that satisfy flow constraints.
	double Flow;

public:
	Edge_t(){};
	Edge_t(int32_t Node_1, int32_t Node_2): Node_1(Node_1), Node_2(Node_2) {IncidentTranscripts.clear(); Flow = 0;};

	bool operator < (const Edge_t& rhs) const {
		if (Node_1 != rhs.Node_1)
			return Node_1 < rhs.Node_1;
		else
			return Node_2 < rhs.Node_2;
	};

	bool operator == (const Edge_t& rhs) const {
		return (Node_1 == rhs.Node_1 && Node_2 == rhs.Node_2);
	};

	void setID (int32_t id) {
		ID = id;
	};
};


class GeneGraph_t{
public:
	string GeneID;
	vector<Node_t> vNodes;
	vector<Edge_t> vEdges;

public:
	GeneGraph_t(){};
	GeneGraph_t(string GeneID, vector<Node_t> vNodes, vector<Edge_t> vEdges): GeneID(GeneID), vNodes(vNodes), vEdges(vEdges) {};
	GeneGraph_t(const vector<Transcript_t>& transcripts);

	void UpdateNodeEdgeReference();
	void CalculateBiasMultiplier(const vector<Transcript_t>& transcripts, const vector< vector<double> >& biascorrections, const vector<double>& expr);

	void RemoveZeroProbNodes();

	string PrintContent(const vector<string>& RefName, const vector<Transcript_t>& transcripts) const;
};


void ConstructAllGeneGraphs(const vector<Transcript_t>& transcripts, const map< string, vector<int32_t> >& GeneTransMap, vector<GeneGraph_t>& GeneGraphs);


void UpdateAllBiasMultiplier(const vector<Transcript_t>& transcripts, const vector< vector<double> >& corrections, const vector<double>& expr,
	const map< string, vector<int32_t> >& GeneTransMap, vector<GeneGraph_t>& GeneGraphs);


void WriteGraph(string outputfile, const vector<GeneGraph_t>& GeneGraphs, const vector<string>& RefName, const vector<Transcript_t>& transcripts);

#endif