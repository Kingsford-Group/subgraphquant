#include "Transcript.hpp"
#include "GeneGraph.hpp"
#include "CombinedBias.hpp"
#include "ProgressBar.hpp"
#include <iomanip>
#include <omp.h>
#include <mutex>

using namespace std;


bool smaller_nodelist(const vector<int32_t>& a, const vector<int32_t>& b)
{
	for (int32_t i = 0; i < a.size() && i < b.size(); i++) {
		if (a[i] != b[i])
			return a[i] < b[i];
	}
	return a.size() < b.size();
};


bool equal_nodelist(const vector<int32_t>& a, const vector<int32_t>& b)
{
	if (a.size() != b.size())
		return false;
	for (int32_t i = 0; i < a.size(); i++) {
		if (a[i] != b[i])
			return false;
	}

	return true;
};


void ReadPrefixPaths(string pathfile, map< string, vector< vector<int32_t> > >& eq_nodes)
{
	eq_nodes.clear();

	ifstream input(pathfile);
	string line;
	while (getline(input, line)) {
		if (line[0] == '#')
			continue;
		vector<string> strs;
		boost::split(strs, line, boost::is_any_of("\t"));
		assert( strs.size() == 2 );
		vector<string> strnode_list;
		boost::split(strnode_list, strs[1], boost::is_any_of(","));
		// type convert
		string g_id = strs[0];
		vector<int32_t> nodes;
		for (string & s : strnode_list)
			nodes.push_back( stoi(s) );
		// add to the map
		map< string, vector< vector<int32_t> > >::iterator it = eq_nodes.find(g_id);
		if (it == eq_nodes.end()) {
			vector< vector<int32_t> > tmp;
			tmp.push_back( nodes );
			eq_nodes[g_id] = tmp;
		}
		else
			(it->second).push_back( nodes );
	}
	input.close();
};


map< string, vector< vector<int32_t> > > Remove_source_sink(const map< string, vector< vector<int32_t> > >& eq_nodes, const vector<GeneGraph_t>& GeneGraphs)
{
	map< string, vector< vector<int32_t> > > out_eq_nodes;

	for (const GeneGraph_t& g : GeneGraphs) {
		// find the index of the source and sink node
		int32_t s = g.vNodes[0].ID;
		int32_t t = g.vNodes.back().ID;

		map< string, vector< vector<int32_t> > >::const_iterator it = eq_nodes.find(g.GeneID);
		vector< vector<int32_t> > tmp_lists;

		if (it != eq_nodes.cend()) {
			for (const vector<int32_t>& nl : (it->second)) {
				tmp_lists.push_back( nl );
				vector<int32_t>& nl2 = tmp_lists.back();
				if (nl2[0] == s)
					nl2.erase( nl2.begin() );
				if (nl2.back() == t)
					nl2.erase( nl2.begin() + nl2.size() - 1 );
			}

			assert( tmp_lists.size() == (it->second).size() );
			out_eq_nodes[g.GeneID] = tmp_lists;
		}
	}
	return out_eq_nodes;
};


void ReadSimpleEdges(string simpleedgefile, map< string, vector< vector<int32_t> > >& se_nodes)
{
	se_nodes.clear();

	ifstream input(simpleedgefile);
	string line;
	while (getline(input, line)) {
		if (line[0] == '#')
			continue;
		vector<string> strs;
		boost::split(strs, line, boost::is_any_of("\t"));

		vector<string> strnode_list;
		boost::split(strnode_list, strs[1], boost::is_any_of(","));
		// type convert
		string g_id = strs[0];
		vector<int32_t> nodes;
		for (string & s : strnode_list)
			nodes.push_back( stoi(s) );
		// add to the map
		map< string, vector< vector<int32_t> > >::iterator it = se_nodes.find(g_id);
		if (it == se_nodes.end()) {
			vector< vector<int32_t> > tmp;
			tmp.push_back( nodes );
			se_nodes[g_id] = tmp;
		}
		else
			(it->second).push_back( nodes );
	}
	input.close();
};


// the following function is deprecated
void ReadEq_GeneNodes(string eqfile, map< string, vector< vector<int32_t> > >& eq_nodes)
{
	time_t CurrentTime;
	string CurrentTimeStr;
	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout<<"["<<CurrentTimeStr.substr(0, CurrentTimeStr.size()-1)<<"] "<<"Reading paths supported by reads in splice graphs."<<endl;

	eq_nodes.clear();
	// read file
	ifstream input(eqfile);
	string line;
	while (getline(input, line)) {
		if (line[0] == '#')
			continue;
		// string split based on \t, the elements are GeneIDs, Aligned_Nodes, Aligned_Edges, ReadCount, Weights
		vector<string> strs;
		boost::split(strs, line, boost::is_any_of("\t"));
		// get the list of genes
		vector<string> gene_list;
		strs[0] = strs[0].substr(1, strs[0].size() - 2);
		boost::split(gene_list, strs[0], boost::is_any_of(","));
		// get the string of nodes corresponding to each gene
		vector<string> strnode_list;
		strs[1] = strs[1].substr(1, strs[1].size() - 2);
		boost::split(strnode_list, strs[1], boost::is_any_of(")"));
		for (int32_t i = 0; i < strnode_list.size(); i++) {
			while (strnode_list[i].front() == '(' || strnode_list[i].front() == ',')
				strnode_list[i] = strnode_list[i].substr(1, strnode_list[i].size() - 1);
			while (strnode_list[i].back() == ')' || strnode_list[i].back() == ',')
				strnode_list[i] = strnode_list[i].substr(0, strnode_list[i].size() - 1);
			if (strnode_list[i] == "") {
				strnode_list.erase(strnode_list.begin() + i);
				i--;
			}
		}
		assert( gene_list.size() == strnode_list.size() );
		// store the int value of node list 
		for (int32_t i = 0; i < gene_list.size(); i++) {
			vector<string> nodes;
			// while (strnode_list[i][0] == '(' || strnode_list[i][0] == ',')
			// 	strnode_list[i] = strnode_list[i].substr(1, strnode_list[i].size() - 1);
			// while (strnode_list[i].back() == ')' || strnode_list[i].back() == ',')
			// 	strnode_list[i] = strnode_list[i].substr(0, strnode_list[i].size() - 1);
			boost::split(nodes, strnode_list[i], boost::is_any_of(","));
			// assert that each eq_class read must cover at least one node
			assert( nodes.size() > 0 );
			vector<int32_t> nodes_indices;
			for (int32_t j = 0; j < nodes.size(); j++)
				nodes_indices.push_back( stoi(nodes[j]) );
			map< string, vector< vector<int32_t> > >::iterator itmap = eq_nodes.find(gene_list[i]);
			if (itmap == eq_nodes.end()) {
				vector< vector<int32_t> > tmp;
				tmp.push_back( nodes_indices );
				eq_nodes[ gene_list[i] ] = tmp;
			}
			else
				(itmap->second).push_back( nodes_indices );
		}
	}
	input.close();

	// remove duplicated node lists
	int32_t count_paths = 0;
	for (map< string, vector< vector<int32_t> > >::iterator itmap = eq_nodes.begin(); itmap != eq_nodes.end(); itmap++) {
		sort((itmap->second).begin(), (itmap->second).end(), smaller_nodelist);
		vector< vector<int32_t> >::iterator ituniq = unique((itmap->second).begin(), (itmap->second).end(), equal_nodelist);
		(itmap->second).resize( distance((itmap->second).begin(),ituniq) );
		count_paths += (itmap->second).size();
	}

	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout << "[" << CurrentTimeStr.substr(0, CurrentTimeStr.size()-1) << "] " << "Finish reading " << count_paths <<" from "<< eq_nodes.size() <<" genes.\n";
};


bool is_path_exist(const vector<int32_t>& nl_query, const vector<int32_t>& nl_ref, const GeneGraph_t& g, pair<int32_t, int32_t>& region)
{
	/*
		Synopsis: take a query path (nodelist) and a reference path (nodelist), return whether the query is a subpath of the reference. If so, find the start/end position (in basepair)
		Input: 
			nl_query: node list of the query. Node list is represented by the index of nodes.
			nl_ref: node list of the reference. Node list is represented by the index of nodes.
			g: Gene graph that encodes the corresponding nodes and edges.
			region: when the query is a subpath of the reference, return the start and end basepair position of the subpath (reference path coordinate)
		Mark:
			In this version, the subpath does not require to be continuous, but can be composed of several blocks (in correponsding topology order). The region is the start of the first block and the end of the last.
	*/
	int32_t idx_query = 0;
	region = make_pair(0, 0);
	int32_t covered_length = 0;
	for (int32_t i = 0; i < nl_ref.size(); i++) {
		if (idx_query < nl_query.size() && nl_ref[i] == nl_query[idx_query]) {
			if (idx_query == 0)
				region.first = covered_length;
			if(idx_query + 1 == nl_query.size())
				region.second = covered_length + g.vNodes[ nl_ref[i] ].EndPos - g.vNodes[ nl_ref[i] ].StartPos;
			idx_query ++;
		}
		covered_length += g.vNodes[ nl_ref[i] ].EndPos - g.vNodes[ nl_ref[i] ].StartPos;
	}

	if (idx_query >= nl_query.size())
		return true;
	else
		return false;
};


vector<double> per_gene_path_bias2(const GeneGraph_t& g, const vector< vector<int32_t> >& path_nodelist, const map<string, string>& trans_seqs, const map<string,double>& trans_expr,
	const map<string,double>& trans_salmon_efflen, const GCBiasModel_t& gcbias, const SeqBiasModel_t& seqbias, const PosBiasModel_t& posbias, const vector<double>& FLD, int32_t fldLow, int32_t fldHigh)
{
	// a map from transcript name to the node list
	map< string, vector<int32_t> > trans_nodelist;
	for (vector<Edge_t>::const_iterator itedge = g.vEdges.cbegin(); itedge != g.vEdges.cend(); itedge++) {
		for (vector<string>::const_iterator ittrans = (itedge->IncidentTranscripts).cbegin(); ittrans != (itedge->IncidentTranscripts).cend(); ittrans++) {
			map< string,vector<int32_t> >::iterator it = trans_nodelist.find( *ittrans );
			if (it == trans_nodelist.end()) {
				vector<int32_t> tmp{itedge->Node_1};
				trans_nodelist[ *ittrans ] = tmp;
			}
			else
				(it->second).push_back( itedge->Node_1 );
		}
	}
	// add the t node to the node list of each transcript
	for (map< string,vector<int32_t> >::iterator it = trans_nodelist.begin(); it != trans_nodelist.end(); it++) {
		(it->second).push_back( g.vNodes.back().ID );
	}

	// result vector
	vector<double> path_efflen(path_nodelist.size(), 0);

	// matrix form of effective length of all positions in each transcript
	map< string, vector< tuple<int32_t,int32_t,double> > > trans_matrix_efflen;

	for (map< string,vector<int32_t> >::iterator it = trans_nodelist.begin(); it != trans_nodelist.end(); it++) {
		const string& t = it->first;
		const string& seq = trans_seqs.at(t);
		const double& efflen = trans_salmon_efflen.at(t);
		if (seq.size() < seqbias.contextLeft+seqbias.contextRight+1)
			continue;
		vector< tuple<int32_t,int32_t,double> > tmp = BiasEffLen_matrix(trans_seqs.at(t), gcbias, seqbias, posbias, FLD, fldLow, fldHigh);
		AdjustBias_single(tmp, t, efflen);

		trans_matrix_efflen[t] = tmp;
	}

	// for each path, find the transcript that it is a substring of
	for (int32_t i = 0; i < path_nodelist.size(); i++) {
		const vector<int32_t>& nl = path_nodelist[i];
		vector<string> related_trans;
		vector<double> related_expr;
		vector<double> related_bias;
		for (map< string,vector<int32_t> >::const_iterator it = trans_nodelist.cbegin(); it != trans_nodelist.cend(); it++) {
			// check whether the path node list falls into transcript node list
			const string& t = it->first;
			const string& seq = trans_seqs.at(t);
			if (seq.size() < seqbias.contextLeft+seqbias.contextRight+1)
				continue;
			pair<int32_t, int32_t> region;
			if ( is_path_exist(nl, it->second, g, region) ) {
				related_trans.push_back( it->first );
				related_expr.push_back( trans_expr.at(it->first) );
				// calculate first and last node length in the query list
				int32_t nodelen_first = g.vNodes[nl.front()].EndPos - g.vNodes[nl.front()].StartPos;
				int32_t nodelen_last = g.vNodes[nl.back()].EndPos - g.vNodes[nl.back()].StartPos;
				// calculate bias
				const vector< tuple<int32_t,int32_t,double> >& matrix_bias = trans_matrix_efflen.at( t );
				double this_bias = 0.0;
				for (vector< tuple<int32_t,int32_t,double> >::const_iterator ittuple = matrix_bias.cbegin(); ittuple != matrix_bias.cend(); ittuple++) {
					if (get<0>(*ittuple) >= region.first && get<0>(*ittuple) < region.first + nodelen_first && get<1>(*ittuple) < region.second && get<1>(*ittuple) >= region.second - nodelen_last)
						this_bias += get<2>(*ittuple);
				}
				related_bias.push_back( this_bias );
			}
		}

		if (related_trans.size() == 0) {
			path_efflen[i] = 0;
			continue;
		}

		assert( related_trans.size() != 0 );
		assert( related_trans.size() == related_expr.size() );
		assert( related_trans.size() == related_bias.size() );
		// weight average of bias weighted by salmon expression (if not all zero), or simple average if the expressions of related transcripts are all zero.
		bool zero_expr = true;
		for (const double& e : related_expr) {
			if (fabs(e) > 1e-4)
				zero_expr = false;
		}
		if (zero_expr) {
			path_efflen[i] = std::accumulate(related_bias.begin(), related_bias.end(), 0.0);
			path_efflen[i] /= related_bias.size();
		}
		else {
			vector<double> ele_mul(related_bias.size(), 0.0);
			std::transform(related_bias.begin(), related_bias.end(), related_expr.begin(), ele_mul.begin(), std::multiplies<double>() );
			path_efflen[i] = std::accumulate(ele_mul.begin(), ele_mul.end(), 0.0) / std::accumulate(related_expr.begin(), related_expr.end(), 0.0);
		}
	}

	return path_efflen;
};


vector<double> per_gene_path_bias3(const GeneGraph_t& g, const vector< vector<int32_t> >& eq_nodelist, 
	const vector< vector<int32_t> >& se_nodelist, const map<string, string>& trans_seqs, const map<string,double>& trans_expr, const map<string,double>& trans_salmon_efflen, 
	const GCBiasModel_t& gcbias, const SeqBiasModel_t& seqbias, const PosBiasModel_t& posbias, const vector<double>& FLD, int32_t fldLow, int32_t fldHigh)
{
	// a map from transcript name to the node list
	map< string, vector<int32_t> > trans_nodelist;
	for (vector<Edge_t>::const_iterator itedge = g.vEdges.cbegin(); itedge != g.vEdges.cend(); itedge++) {
		for (vector<string>::const_iterator ittrans = (itedge->IncidentTranscripts).cbegin(); ittrans != (itedge->IncidentTranscripts).cend(); ittrans++) {
			map< string,vector<int32_t> >::iterator it = trans_nodelist.find( *ittrans );
			if (it == trans_nodelist.end()) {
				vector<int32_t> tmp{itedge->Node_1};
				trans_nodelist[ *ittrans ] = tmp;
			}
			else
				(it->second).push_back( itedge->Node_1 );
		}
	}
	// add the t node to the node list of each transcript
	for (map< string,vector<int32_t> >::iterator it = trans_nodelist.begin(); it != trans_nodelist.end(); it++) {
		(it->second).push_back( g.vNodes.back().ID );
	}

	// result vector
	vector<double> eq_nodes_efflen(eq_nodelist.size(), 0);
	vector<double> se_nodes_efflen(se_nodelist.size(), 0);

	// matrix form of effective length of all positions in each transcript
	map< string, vector< tuple<int32_t,int32_t,double> > > trans_matrix_efflen;

	// position-wise bias correction/efflen for reference transcripts
	for (map< string,vector<int32_t> >::iterator it = trans_nodelist.begin(); it != trans_nodelist.end(); it++) {
		const string& t = it->first;
		const string& seq = trans_seqs.at(t);
		const double& efflen = trans_salmon_efflen.at(t);
		if (seq.size() < seqbias.contextLeft+seqbias.contextRight+1)
			continue;
		vector< tuple<int32_t,int32_t,double> > tmp = BiasEffLen_matrix(trans_seqs.at(t), gcbias, seqbias, posbias, FLD, fldLow, fldHigh);
		AdjustBias_single(tmp, t, efflen);

		trans_matrix_efflen[t] = tmp;
	}

	// for each path, find the transcript that it is a substring of
	for (int32_t i = 0; i < eq_nodelist.size(); i++) {
		const vector<int32_t>& nl = eq_nodelist[i];
		vector<string> related_trans;
		vector<double> related_expr;
		vector<double> related_bias;
		for (map< string,vector<int32_t> >::const_iterator it = trans_nodelist.cbegin(); it != trans_nodelist.cend(); it++) {
			const string& t = it->first;
			const string& seq = trans_seqs.at(t);
			if (seq.size() < seqbias.contextLeft+seqbias.contextRight+1)
				continue;
			pair<int32_t, int32_t> region;
			// check whether the path node list falls into transcript node list
			if ( is_path_exist(nl, it->second, g, region) ) {
				related_trans.push_back( it->first );
				related_expr.push_back( max(1.0, trans_expr.at(it->first)) );
				// calculate first and last node length in the query list
				int32_t nodelen_first = g.vNodes[nl.front()].EndPos - g.vNodes[nl.front()].StartPos;
				int32_t nodelen_last = g.vNodes[nl.back()].EndPos - g.vNodes[nl.back()].StartPos;
				// calculate bias
				const vector< tuple<int32_t,int32_t,double> >& matrix_bias = trans_matrix_efflen.at( t );
				double this_bias = 0.0;
				for (vector< tuple<int32_t,int32_t,double> >::const_iterator ittuple = matrix_bias.cbegin(); ittuple != matrix_bias.cend(); ittuple++) {
					if (get<0>(*ittuple) >= region.first && get<0>(*ittuple) < region.first + nodelen_first && get<1>(*ittuple) < region.second && get<1>(*ittuple) >= region.second - nodelen_last)
						this_bias += get<2>(*ittuple);
				}
				// test: hard clip on extremely small efflen for single node
				if (nl.size() == 1 && nodelen_first > fldLow && (nodelen_first - fldLow) * 0.01 > this_bias ) {
					// cout << "hard clipping bias for " << g.GeneID << " node " << nl[0] << " old bias = " << this_bias << " new bias = " << (nodelen_first - fldLow) * 0.01 << endl;
					this_bias = (nodelen_first - fldLow) * 0.01;
				}
				related_bias.push_back( this_bias );
			}
		}

		if (related_trans.size() == 0) {
			eq_nodes_efflen[i] = 0;
			continue;
		}

		assert( related_trans.size() != 0 );
		assert( related_trans.size() == related_expr.size() );
		assert( related_trans.size() == related_bias.size() );
		// weight average of bias weighted by salmon expression (if not all zero), or simple average if the expressions of related transcripts are all zero.
		bool zero_expr = true;
		for (const double& e : related_expr) {
			if (fabs(e) > 1e-4)
				zero_expr = false;
		}
		if (zero_expr) {
			eq_nodes_efflen[i] = std::accumulate(related_bias.begin(), related_bias.end(), 0.0);
			eq_nodes_efflen[i] /= related_bias.size();
		}
		else {
			vector<double> ele_mul(related_bias.size(), 0.0);
			std::transform(related_bias.begin(), related_bias.end(), related_expr.begin(), ele_mul.begin(), std::multiplies<double>() );
			eq_nodes_efflen[i] = std::accumulate(ele_mul.begin(), ele_mul.end(), 0.0) / std::accumulate(related_expr.begin(), related_expr.end(), 0.0);
		}
	}

	// for each path, find the transcript that it is a substring of
	for (int32_t i = 0; i < se_nodelist.size(); i++) {
		const vector<int32_t>& nl = se_nodelist[i];
		vector<string> related_trans;
		vector<double> related_expr;
		vector<double> related_bias;
		for (map< string,vector<int32_t> >::const_iterator it = trans_nodelist.cbegin(); it != trans_nodelist.cend(); it++) {
			const string& t = it->first;
			const string& seq = trans_seqs.at(t);
			if (seq.size() < seqbias.contextLeft+seqbias.contextRight+1)
				continue;
			pair<int32_t, int32_t> region;
			// check whether the path node list falls into transcript node list
			if ( is_path_exist(nl, it->second, g, region) ) {
				related_trans.push_back( it->first );
				related_expr.push_back( trans_expr.at(it->first) );
				// calculate first and last node length in the query list
				int32_t nodelen_first = g.vNodes[nl.front()].EndPos - g.vNodes[nl.front()].StartPos;
				int32_t nodelen_last = g.vNodes[nl.back()].EndPos - g.vNodes[nl.back()].StartPos;
				// calculate bias
				const vector< tuple<int32_t,int32_t,double> >& matrix_bias = trans_matrix_efflen.at( t );
				double this_bias = 0.0;
				for (vector< tuple<int32_t,int32_t,double> >::const_iterator ittuple = matrix_bias.cbegin(); ittuple != matrix_bias.cend(); ittuple++) {
					if (get<0>(*ittuple) >= region.first && get<0>(*ittuple) < region.first + nodelen_first && get<1>(*ittuple) < region.second && get<1>(*ittuple) >= region.second - nodelen_last)
						this_bias += get<2>(*ittuple);
				}
				related_bias.push_back( this_bias );
			}
		}

		if (related_trans.size() == 0) {
			se_nodes_efflen[i] = 0;
			continue;
		}

		assert( related_trans.size() != 0 );
		assert( related_trans.size() == related_expr.size() );
		assert( related_trans.size() == related_bias.size() );
		// weight average of bias weighted by salmon expression (if not all zero), or simple average if the expressions of related transcripts are all zero.
		bool zero_expr = true;
		for (const double& e : related_expr) {
			if (fabs(e) > 1e-4)
				zero_expr = false;
		}
		if (zero_expr) {
			se_nodes_efflen[i] = std::accumulate(related_bias.begin(), related_bias.end(), 0.0);
			se_nodes_efflen[i] /= related_bias.size();
		}
		else {
			vector<double> ele_mul(related_bias.size(), 0.0);
			std::transform(related_bias.begin(), related_bias.end(), related_expr.begin(), ele_mul.begin(), std::multiplies<double>() );
			se_nodes_efflen[i] = std::accumulate(ele_mul.begin(), ele_mul.end(), 0.0) / std::accumulate(related_expr.begin(), related_expr.end(), 0.0);
		}
	}

	// return: concatenation of eq_nodes effective length and simple_edge_nodes effective length
	eq_nodes_efflen.insert(eq_nodes_efflen.end(), se_nodes_efflen.begin(), se_nodes_efflen.end());
	return eq_nodes_efflen;
};


void process_hyperedge_efflen(map< string,vector<double> >& result_eq_efflen, map< string,vector<double> >& result_se_efflen, 
	const vector<GeneGraph_t>& GeneGraphs, const map< string, vector< vector<int32_t> > >& eq_nodes, const map< string, vector< vector<int32_t> > >& se_nodes, 
	const map< string,int32_t >& TransIndex, const vector<string>& sequences, const vector<double>& expr, const vector<double>& salmon_efflen, 
	const GCBiasModel_t& gcbias, const SeqBiasModel_t& seqbias, const PosBiasModel_t& posbias, const vector<double>& FLD, int32_t fldLow, int32_t fldHigh, int32_t threads=4)
{
	time_t CurrentTime;
	string CurrentTimeStr;
	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout<<"["<<CurrentTimeStr.substr(0, CurrentTimeStr.size()-1)<<"] "<<"Calculating bias-corrected length."<<endl;

	// result
	result_eq_efflen.clear();
	result_se_efflen.clear();

	ProgressBar progressBar(GeneGraphs.size(), 70);

	// loop over each gene
	mutex corrections_mutex;
	omp_set_num_threads(8);
	#pragma omp parallel for
	// for (int32_t i = 0; i < 100; i++) {
	for (int32_t i = 0; i < GeneGraphs.size(); i++) {
		const GeneGraph_t& g = GeneGraphs[i];
		map<string, string> trans_seqs;
		map<string, double> trans_expr;
		map<string, double> trans_salmon_efflen;
		for (const Edge_t& e : g.vEdges) {
			for (const string& t : e.IncidentTranscripts) {
				map<string, int32_t>::const_iterator ittrans = TransIndex.find(t);
				assert( ittrans != TransIndex.cend() );
				trans_seqs[t] = sequences[ ittrans->second ];
				trans_expr[t] = expr[ ittrans->second ];
				trans_salmon_efflen[t] = salmon_efflen[ ittrans->second ];
			}
		}
		vector<double> eq_nodes_efflen;
		vector<double> se_nodes_efflen;
		map< string, vector< vector<int32_t> > >::const_iterator it_eq = eq_nodes.find(g.GeneID);
		map< string, vector< vector<int32_t> > >::const_iterator it_se = se_nodes.find(g.GeneID);
		if (it_eq != eq_nodes.cend()) {
			if (it_se == se_nodes.cend())
				cout << g.GeneID << endl;
			assert(it_se != se_nodes.cend());
			// efflen = per_gene_path_bias2(g, it->second, trans_seqs, trans_expr, trans_salmon_efflen, gcbias, seqbias, posbias, FLD, fldLow, fldHigh);
			eq_nodes_efflen = per_gene_path_bias3(g, it_eq->second, it_se->second, trans_seqs, trans_expr, trans_salmon_efflen, gcbias, seqbias, posbias, FLD, fldLow, fldHigh);
			int32_t n_lists_eq = (it_eq->second).size();
			int32_t n_lists_se = (it_se->second).size();
			se_nodes_efflen.insert(se_nodes_efflen.end(), eq_nodes_efflen.begin() + n_lists_eq, eq_nodes_efflen.end());
			assert(eq_nodes_efflen.size() == n_lists_eq + n_lists_se);
			assert(se_nodes_efflen.size() == n_lists_se);
			eq_nodes_efflen.resize(n_lists_eq);
		}
		
		lock_guard<std::mutex> guard(corrections_mutex);
		result_eq_efflen[g.GeneID] = eq_nodes_efflen;
		result_se_efflen[g.GeneID] = se_nodes_efflen;
		if (it_eq != eq_nodes.cend() && eq_nodes_efflen.size() != (it_eq->second).size())
			cout << "watch here\n";
		if (it_eq != eq_nodes.cend() && se_nodes_efflen.size() != (it_se->second).size())
			cout << "watch here.\n";
		++progressBar;
		progressBar.display();
	}

	progressBar.done();
	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout << "[" << CurrentTimeStr.substr(0, CurrentTimeStr.size()-1) << "] " << "Finish.\n";
};


bool sanity_check(const vector<GeneGraph_t>& GeneGraphs, const map< string, vector< vector<int32_t> > >& eq_nodes)
{
	for (const GeneGraph_t& g : GeneGraphs) {
		// a map from transcript name to the node list
		map< string, vector<int32_t> > trans_nodelist;
		for (vector<Edge_t>::const_iterator itedge = g.vEdges.cbegin(); itedge != g.vEdges.cend(); itedge++) {
			for (vector<string>::const_iterator ittrans = (itedge->IncidentTranscripts).cbegin(); ittrans != (itedge->IncidentTranscripts).cend(); ittrans++) {
				map< string,vector<int32_t> >::iterator it = trans_nodelist.find( *ittrans );
				if (it == trans_nodelist.end()) {
					vector<int32_t> tmp{itedge->Node_1};
					trans_nodelist[ *ittrans ] = tmp;
				}
				else
					(it->second).push_back( itedge->Node_1 );
			}
		}
		// add the t node to the node list of each transcript
		for (map< string,vector<int32_t> >::iterator it = trans_nodelist.begin(); it != trans_nodelist.end(); it++) {
			(it->second).push_back( g.vNodes.back().ID );
		}

		// find paths
		map< string, vector< vector<int32_t> > >::const_iterator it_nodes = eq_nodes.find(g.GeneID);
		if (it_nodes == eq_nodes.cend())
			continue;
		const vector< vector<int32_t> >& path_nodelist = it_nodes->second;

		for (int32_t i = 0; i < path_nodelist.size(); i++) {
			bool path_exist = false;
			const vector<int32_t>& nl = path_nodelist[i];
			for (map< string,vector<int32_t> >::const_iterator it = trans_nodelist.cbegin(); it != trans_nodelist.cend(); it++) {
				// check whether the path node list falls into transcript node list
				pair<int32_t, int32_t> region;
				if ( is_path_exist(nl, it->second, g, region) ) 
					path_exist = true;
			}
			if (! path_exist) {
				cout << "watch here\n";
				return false;
			}
		}
	}

	return true;
};


void write_hyperedge_efflen(string outputfile, const map< string,vector<double> >& hyperedge_efflen, const map< string, vector< vector<int32_t> > >& eq_nodes)
{
	ofstream output(outputfile, ios::out);
	output << "# gene\tnode_list\tefflen\n";

	for (map< string, vector< vector<int32_t> > >::const_iterator it = eq_nodes.cbegin(); it != eq_nodes.cend(); it++) {
		// find the hyperedge effective length
		map< string,vector<double> >::const_iterator it_len = hyperedge_efflen.find( it->first );
		if (it_len == hyperedge_efflen.cend() || (it_len->second).size() == 0)
			continue;
		// assert( it_len != hyperedge_efflen.cend() );
		if ( (it_len->second).size() != (it->second).size() )
			cout << "watch here\t" << (it->first) <<"\n";
		assert( (it_len->second).size() == (it->second).size() );

		for (int32_t i = 0; i < (it->second).size(); i++) {
			const vector<int32_t>& nl = (it->second)[i];
			const double& this_efflen = (it_len->second)[i];

			// get the string of nodelist
			string strnode_list = "";
			for (const int32_t& n : nl)
				strnode_list += to_string(n) + ",";
			if (strnode_list.back() == ',')
				strnode_list = strnode_list.substr(0, strnode_list.size() - 1);
			// write to file
			output << it->first <<"\t"<< strnode_list <<"\t"<< std::setprecision(9) << this_efflen << endl;
		}
	}

	output.close();
};


void write_path_efflen_inorder(string outputfile, const map< string,vector<double> >& result_efflen, const map< string, vector< vector<int32_t> > >& eq_nodes, 
	const map< string, vector< vector<int32_t> > >& original_eq_nodes, const vector<GeneGraph_t>& GeneGraphs)
{
	// map from gene name to index in GeneGraphs
	map< string,int32_t > name_index;
	for (int32_t i = 0; i < GeneGraphs.size(); i++)
		name_index[ GeneGraphs[i].GeneID ] = i;

	ofstream output(outputfile, ios::out);
	output << "# gene\tnode_list\tefflen\n";

	int32_t count_duplicate = 0;

	for (map< string, vector< vector<int32_t> > >::const_iterator it = eq_nodes.cbegin(); it != eq_nodes.cend(); it++) {

		// find the original node list corresponding to this gene that contain the source and sink
		map< string, vector< vector<int32_t> > >::const_iterator it_ori = original_eq_nodes.find( it->first );
		assert( it_ori != original_eq_nodes.cend() );

		// find the path effective length
		map< string,vector<double> >::const_iterator it_len = result_efflen.find( it->first );
		assert( it_len != result_efflen.cend() );
		assert( (it->second).size() == (it_ori->second).size() );
		assert( (it_len->second).size() == (it->second).size() );

		// find the corresponding gene and source and sink
		const GeneGraph_t& g = GeneGraphs[ name_index.at(it->first) ];
		int32_t s = g.vNodes[0].ID;
		int32_t t = g.vNodes.back().ID;

		for (const vector<int32_t>& ori_nl : (it_ori->second)) {
			vector<int32_t> nl(ori_nl.cbegin(), ori_nl.cend());
			// get the version without source/sink
			if (nl[0] == s)
				nl.erase( nl.begin() );
			if (nl.back() == t)
				nl.erase( nl.begin() + nl.size() - 1 );
			// find in the current node list
			int32_t idx = -1;
			for (idx = 0; idx < (it->second).size(); idx++) {
				if (equal_nodelist(nl, (it->second)[idx]))
					break;
			}
			assert( idx >= 0 && idx < (it->second).size() );
			// find the effective length
			const double& this_efflen = (it_len->second)[idx];
			// a special case: if nl ends with t, but after removing t, nl is an duplicated entry, then don't count the effective length again
			bool is_duplicated = false;
			if (ori_nl.back() == t) {
				// check whether nl is duplicated
				idx++;
				for (; idx < (it->second).size(); idx++)
					if (equal_nodelist(nl, (it->second)[idx])) {
						is_duplicated = true;
						count_duplicate++;
						break;
					}
			}
			// get the string of nodelist
			string strnode_list = "";
			for (const int32_t& n : ori_nl)
				strnode_list += to_string(n) + ",";
			if (strnode_list.back() == ',')
				strnode_list = strnode_list.substr(0, strnode_list.size() - 1);
			// write to file
			if (!is_duplicated)
				output << it->first <<"\t"<< strnode_list <<"\t"<< std::setprecision(9) << this_efflen << endl;
			else
				output << it->first <<"\t"<< strnode_list <<"\t"<< std::setprecision(9) << 0 << endl;
		}
	}
	output.close();

	cout << "Find " << count_duplicate << " entries after removing T.\n";
};


int32_t main(int32_t argc, char* argv[])
{
	if (argc == 1) {
		printf("pathbias <GTF file> <Genome Fasta> <Salmon quantification> <graph prefix>\n");
		exit(-1);
	}

	string gtffile = argv[1];
	string genomefasta = argv[2];
	string salmonquantfile = argv[3];
	string prefix = argv[4];

	string graphfile = prefix + "_graph_fragstart_beforebias.txt";
	string simpleedgefile = prefix + "_simpleedge_to_path.txt";
	string output_efflen_hyperedge = prefix + "_hyperedge_efflen.txt";
	string output_efflen_simpleedge = prefix + "_simpleedge_efflen.txt";
	string output_efflen_prefix = prefix + "_path_efflen.txt";
	size_t t = prefix.find_last_of("/");
	string eqfile = prefix.substr(0, t) + "/eq_classes.txt";

	size_t p_folder = salmonquantfile.find_last_of("/");
	string salmonauxfolder = salmonquantfile.substr(0, p_folder) + "/aux_info/";

	vector<string> RefName;
	std::map<string,uint8_t> RefTable;
	BuildRefName(gtffile, RefName, RefTable);

	// read GTF and retrieve transcript sequence
	vector<Transcript_t> Transcripts;
	ReadGTF (gtffile, Transcripts, RefTable);
	vector<string> sequences;
	ReadTransSequence_genome(genomefasta, RefTable, Transcripts, sequences);

	// gene to transcript map and transcript to its index
	map< string, vector<int32_t> > GeneTransMap;
	map< string,int32_t > TransIndex;
	Map_Gene_Trans (Transcripts, GeneTransMap);
	MapIndex(Transcripts, TransIndex);

	// read salmon quantification TPM
	vector<double> expr;
	ReadSalmonTPM(salmonquantfile, TransIndex, expr);
	vector<double> salmon_efflen;
	ReadSalmonEfflen(salmonquantfile, TransIndex, salmon_efflen);

	// read gene graphs
	vector<GeneGraph_t> GeneGraphs;
	ReadGraph(graphfile, GeneGraphs, RefTable);

	map< string, vector< vector<int32_t> > > eq_nodes;

	// sanity check on the equivalent class file
	ReadEq_GeneNodes(eqfile, eq_nodes);
	bool all_path_exist = sanity_check(GeneGraphs, eq_nodes);
	assert(all_path_exist);

	// read node lists corresponding to simple edges (single node or pair of connected nodes
	map< string, vector< vector<int32_t> > > se_nodes;
	ReadSimpleEdges(simpleedgefile, se_nodes);

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


	map< string,vector<double> > hyperedge_efflen;
	map< string,vector<double> > simpleedge_efflen;
	process_hyperedge_efflen(hyperedge_efflen, simpleedge_efflen, GeneGraphs, eq_nodes, se_nodes, TransIndex, sequences, expr, 
		salmon_efflen, gcbias, seqbias, posbias, FLD, fldLow, fldHigh);

	write_hyperedge_efflen(output_efflen_hyperedge, hyperedge_efflen, eq_nodes);
	write_hyperedge_efflen(output_efflen_simpleedge, simpleedge_efflen, se_nodes);

	// ReadPrefixPaths(pathfile, prefix_paths);

	// write_path_efflen_inorder(output_efflen_prefix, result_efflen, eq_nodes, original_eq_nodes, GeneGraphs);
};
