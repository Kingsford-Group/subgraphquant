#ifndef __Transcript_H__
#define __Transcript_H__

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <string>
#include "boost/algorithm/string.hpp"
#include "utils.hpp"

using namespace std;


class Exon_t {
public:
	int32_t StartPos;
	int32_t EndPos;

public:
	Exon_t () {};
	Exon_t (int32_t StartPos, int32_t EndPos): StartPos(StartPos), EndPos(EndPos) {};

	bool operator < (const Exon_t& rhs) const {
		if (StartPos != rhs.StartPos)
			return StartPos < rhs.StartPos;
		else
			return EndPos < rhs.EndPos;
	};

	bool operator == (const Exon_t& rhs) const {
		return (StartPos == rhs.StartPos) && (EndPos == rhs.EndPos);
	};
};


class Transcript_t {
public:
	uint8_t Chr;
	int32_t StartPos;
	int32_t EndPos;
	bool Strand;
	string TransID;
	string GeneID;
	vector<Exon_t> vExons;

public:
	Transcript_t () {};
	Transcript_t (uint8_t Chr, int32_t StartPos, int32_t EndPos, bool Strand, string TransID, string GeneID):
		Chr(Chr), StartPos(StartPos), EndPos(EndPos), Strand(Strand), TransID(TransID), GeneID(GeneID) {};
	Transcript_t (string gtfline, const map<string,uint8_t>& ChrTable);

	bool operator < (const Transcript_t& rhs) const {
		if (Chr != rhs.Chr)
			return Chr < rhs.Chr;
		else if (StartPos != rhs.StartPos)
			return StartPos < rhs.StartPos;
		else
			return EndPos < rhs.EndPos;
	};

	static bool CompareTransID (const Transcript_t& lhs, const Transcript_t& rhs) {
		return lhs.TransID < rhs.TransID;
	};

	bool IsWithinExon (pair<int32_t, int32_t> region);
};


void BuildRefName(string gtffile, vector<string>& RefName, std::map<string,uint8_t> & RefTable);


string GetFeature (string line, string key);


int32_t ReadGTF (string gtffile, vector<Transcript_t>& Transcripts, const map<string, uint8_t>& ChrTable);


int32_t Map_Gene_Trans (const vector<Transcript_t>& Transcripts, map< string, vector<int32_t> >& GeneTransMap);


int32_t MapIndex (const vector<Transcript_t>& Transcripts, map< string,int32_t >& TransIndex);


void ReadTransSequence_transcript(string transcriptfasta, const vector<Transcript_t>& Transcripts, vector<string>& sequences);


void ReadTransSequence_genome(string genomefasta, const map<string, uint8_t>& ChrTable, const vector<Transcript_t>& Transcripts, vector<string>& sequences);


#endif
