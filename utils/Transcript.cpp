#include "Transcript.hpp"

Transcript_t::Transcript_t (string gtfline, const map<string,uint8_t>& ChrTable)
{
	vector<string> strs;
	boost::split(strs, gtfline, boost::is_any_of("\t"));
	assert(strs[2] == "transcript");

	// position info
	map<string, uint8_t>::const_iterator it = ChrTable.find(strs[0]);
	assert(it != ChrTable.cend());
	Chr = it->second;
	StartPos = stoi(strs[3]) - 1;
	EndPos = stoi(strs[4]);
	Strand = (strs[6] == "+");

	// TransID and GeneID
	TransID = GetFeature(gtfline, "transcript_id");
	GeneID = GetFeature(gtfline, "gene_id");

	vExons.clear();
};


bool Transcript_t::IsWithinExon (pair<int32_t, int32_t> region)
{
	// check the large range of transcript, if region doesn't fall into transcript, just return false
	if (StartPos >= region.first || EndPos <= region.second)
		return false;
	// go through each exon, check whether the region is totally within an exon
	for (vector<Exon_t>::const_iterator it = vExons.cbegin(); it != vExons.cend(); it++) {
		if (it->StartPos <= region.first && it->EndPos >= region.second)
			return true;
	}

	return false;
};


void BuildRefName(string gtffile, vector<string>& RefName, std::map<string,uint8_t> & RefTable)
{
	time_t CurrentTime;
	string CurrentTimeStr;
	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout<<"["<<CurrentTimeStr.substr(0, CurrentTimeStr.size()-1)<<"] "<<"Constructing reference table from GTF."<<endl;

	// clear variables
	RefName.clear();
	RefTable.clear();
	// construct ref info from GTF file
	ifstream input(gtffile);
	string line;
	while (getline(input,line)) {
		if (line[0] == '#')
			continue;
		vector<string> strs;
		boost::split(strs, line, boost::is_any_of("\t"));
		if (RefName.size() == 0 || RefName.back() != strs[0])
			RefName.push_back(strs[0]);
	}
	input.close();
	// sort and get unique RefName
	sort(RefName.begin(), RefName.end());
	RefName.resize( distance(RefName.begin(), unique(RefName.begin(),RefName.end())) );
	// construct RefTable
	for (uint8_t i = 0; i < RefName.size(); i++)
		RefTable[RefName[i]] = i;

	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout << "[" << CurrentTimeStr.substr(0, CurrentTimeStr.size()-1) << "] " << "Finish. Number of chromosomes = " << RefTable.size() <<"\n";
};


string GetFeature (string line, string key)
{
	int32_t s = line.find(key);
	if (s == string::npos)
		return "";
	int32_t t = line.find(";", s+1);
	if (t == string::npos)
		return "";

	return line.substr(s + (int32_t)key.size() + 2, t - s - (int32_t)key.size() - 3);
};


int32_t ReadGTF (string gtffile, vector<Transcript_t>& Transcripts, const map<string, uint8_t>& ChrTable)
{
	time_t CurrentTime;
	string CurrentTimeStr;
	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout<<"["<<CurrentTimeStr.substr(0, CurrentTimeStr.size()-1)<<"] "<<"Reading annotation GTF."<<endl;

	// clear variable
	Transcripts.clear();
	// in case there is only exon records, no transcript lines, create a transcript for each exon, and then merge
	vector<Transcript_t> trans_for_exon;

	ifstream input(gtffile);
	string line;

	while (getline(input, line)) {
		if (line[0] == '#')
			continue;
		vector<string> strs;
		boost::split(strs, line, boost::is_any_of("\t"));
		map<string, uint8_t>::const_iterator it = ChrTable.find(strs[0]);
		if (it == ChrTable.cend())
			continue;
		if (strs[2] == "transcript") {
			Transcript_t t(line, ChrTable);
			Transcripts.push_back( t );
		}
		else if (strs[2] == "exon") {
			string transid = GetFeature(line, "transcript_id");
			string geneid = GetFeature(line, "gene_id");
			Exon_t e( stoi(strs[3]) - 1, stoi(strs[4]) );
			if (Transcripts.size() > 0 && Transcripts.back().TransID == transid)
				Transcripts.back().vExons.push_back( e );
			else{
				map<string, uint8_t>::const_iterator it = ChrTable.find(strs[0]);
				assert(it != ChrTable.cend());
				uint8_t Chr = it->second;
				Transcript_t t(Chr, e.StartPos, e.EndPos, (strs[6] == "+"), transid, geneid);
				t.vExons.push_back( e );
				trans_for_exon.push_back( t );
			}
		}
	}

	input.close();

	// for the extra exons, sort based on transcript id, and merge
	sort(trans_for_exon.begin(), trans_for_exon.end(), Transcript_t::CompareTransID);

	// merge
	int32_t s = 0;
	int32_t t = s+1;
	while (s < trans_for_exon.size()) {
		for (t = s+1; t < trans_for_exon.size() && trans_for_exon[t].TransID == trans_for_exon[s].TransID; t++) {
			assert( trans_for_exon[t].vExons.size() == 1 );
			trans_for_exon[s].vExons.push_back( trans_for_exon[t].vExons[0] );
		}
		// sort the exons in order, and adjust the transcript StartPos and EndPos
		assert(trans_for_exon[s].vExons.size() > 0);
		sort(trans_for_exon[s].vExons.begin(), trans_for_exon[s].vExons.end());
		trans_for_exon[s].StartPos = trans_for_exon[s].vExons.front().StartPos;
		trans_for_exon[s].EndPos = trans_for_exon[s].vExons.back().EndPos;
		Transcripts.push_back( trans_for_exon[s] );
	}

	// sort exon based in order / reverse order based on strand
	for (vector<Transcript_t>::iterator it = Transcripts.begin(); it != Transcripts.end(); it++) {
		sort((it->vExons).begin(), (it->vExons).end());
		assert( it->StartPos <= (it->vExons).front().StartPos && it->EndPos >=(it->vExons).back().EndPos );
		if (!(it->Strand))
			reverse((it->vExons).begin(), (it->vExons).end());
		(it->vExons).reserve( (it->vExons).size() );
	}

	// sort all transcripts in coordinate order
	sort(Transcripts.begin(), Transcripts.end());

	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout << "[" << CurrentTimeStr.substr(0, CurrentTimeStr.size()-1) << "] " << "Finish. Number of transcripts = " << Transcripts.size() <<"\n";

	return 0;
};


int32_t Map_Gene_Trans (const vector<Transcript_t>& Transcripts, map< string, vector<int32_t> >& GeneTransMap)
{
	// clear variable
	GeneTransMap.clear();

	// for each gene, find the corresponding transcripts, store the indexes in the map
	for (int32_t i = 0; i < Transcripts.size(); i++) {
		const string& geneid = Transcripts[i].GeneID;
		map< string, vector<int32_t> >::iterator it = GeneTransMap.find(geneid);
		if (it != GeneTransMap.end())
			(it->second).push_back(i);
		else {
			vector<int32_t> tmp;
			tmp.push_back(i);
			GeneTransMap[geneid] = tmp;
		}
	}

	return 0;
};


int32_t MapIndex (const vector<Transcript_t>& Transcripts, map< string,int32_t >& TransIndex)
{
	TransIndex.clear();
	for (int32_t i = 0; i < Transcripts.size(); i++)
		TransIndex[Transcripts[i].TransID] = i;
	return 0;
};


void ReadTransSequence_transcript(string transcriptfasta, const vector<Transcript_t>& Transcripts, vector<string>& sequences)
{
	time_t CurrentTime;
	string CurrentTimeStr;
	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout<<"["<<CurrentTimeStr.substr(0, CurrentTimeStr.size()-1)<<"] "<<"Reading transcript sequence."<<endl;

	// clear variables, require sequences in the same order as Transcripts
	sequences.clear();
	sequences.assign(Transcripts.size(), "");
	// mapping transcript ID to its index in Transcripts
	map<string,int32_t> TransIndex;
	for (int32_t i = 0; i < Transcripts.size(); i++)
		TransIndex[Transcripts[i].TransID] = i;
	// read transcript fasta file
	ifstream input(transcriptfasta);
	string line;
	string transname = "";
	string transseq = "";
	while(getline(input, line)) {
		if (line[0] == '>') {
			// adding the transcript sequence to record
			if (transname != "") {
				assert(transseq != "");
				map<string,int32_t>::const_iterator itidx = TransIndex.find(transname);
				if (itidx != TransIndex.cend())
					sequences[itidx->second] = transseq;
			}
			// get new transcript name
			vector<string> strs;
			boost::split(strs, line, boost::is_any_of("\t |"));
			transname = strs[0].substr(1);
			transseq = "";
		}
		else 
			transseq += line;
	}
	// adding last record
	if (transname != "") {
		assert(transseq != "");
		map<string,int32_t>::const_iterator itidx = TransIndex.find(transname);
		if (itidx != TransIndex.cend())
			sequences[itidx->second] = transseq;
	}
	// sanity checks
	// sanity check about the sequence length
	for (int32_t i = 0; i < sequences.size(); i++) {
		int32_t len_seq = sequences[i].size();
		int32_t len_gtf = 0;
		for (const Exon_t& e : Transcripts[i].vExons)
			len_gtf += e.EndPos - e.StartPos;
		if (len_seq != len_gtf)
			cout << "Transcript length in FASTA does not match its length in gene annotation: " << Transcripts[i].TransID <<" "<< len_seq <<"!="<< len_gtf << endl;
		assert(len_seq == len_gtf);
	}
	input.close();

	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout << "[" << CurrentTimeStr.substr(0, CurrentTimeStr.size()-1) << "] " << "Finish.\n";
};


void ReadTransSequence_genome(string genomefasta, const map<string, uint8_t>& ChrTable, const vector<Transcript_t>& Transcripts, 
	vector<string>& sequences)
{
	time_t CurrentTime;
	string CurrentTimeStr;
	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout<<"["<<CurrentTimeStr.substr(0, CurrentTimeStr.size()-1)<<"] "<<"Retrieving transcript sequence from genome."<<endl;

	// clear variables, require sequences in the same order as Transcripts
	sequences.clear();
	sequences.assign(Transcripts.size(), "");

	// read genome sequence
	vector<string> g(ChrTable.size(), "");
	ifstream input(genomefasta);
	string line;
	string chrname = "";
	string chrseq = "";
	while (getline(input, line)) {
		if (line[0] == '>') {
			if (chrname != "") {
				map<string,uint8_t>::const_iterator itchr = ChrTable.find(chrname);
				if (itchr != ChrTable.cend())
					g[itchr->second] = chrseq;
			}
			// read new chrname and chrseq
			vector<string> strs;
			boost::split(strs, line, boost::is_any_of("\t |"));
			chrname = strs[0].substr(1);
			chrseq = "";
		}
		else
			chrseq += line;
	}
	// store the last chr
	if (chrname != "") {
		map<string,uint8_t>::const_iterator itchr = ChrTable.find(chrname);
		if (itchr != ChrTable.cend())
			g[itchr->second] = chrseq;
	}
	input.close();

	// substrings of genome are the transcripts
	for (int32_t i = 0; i < Transcripts.size(); i++) {
		const Transcript_t& t = Transcripts[i];
		int32_t idxchr = t.Chr;
		string transseq = "";
		if (t.StartPos < 0 || t.EndPos > g[idxchr].size())
			cout << "watch here\n";
		assert(t.StartPos >= 0 && t.EndPos <= g[idxchr].size());
		if (t.Strand) {
			for (const Exon_t& e : t.vExons)
				transseq += g[idxchr].substr(e.StartPos, e.EndPos-e.StartPos);
		}
		else {
			for (vector<Exon_t>::const_reverse_iterator it = t.vExons.rbegin(); it != t.vExons.rend(); it++)
				transseq += g[idxchr].substr(it->StartPos, it->EndPos - it->StartPos);
			transseq = ReverseComplement(transseq.begin(), transseq.end());
		}
		sequences[i] = transseq;
	}

	// sanity check on length
	for (int32_t i = 0; i < Transcripts.size(); i++) {
		const Transcript_t& t = Transcripts[i];
		int32_t length = 0;
		for (const Exon_t& e : t.vExons)
			length += e.EndPos - e.StartPos;
		assert(length == sequences[i].size());
	}

	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout << "[" << CurrentTimeStr.substr(0, CurrentTimeStr.size()-1) << "] " << "Finish.\n";
};

