#include "GeneGraph.hpp"
#include "Transcript.hpp"

GeneGraph_t::GeneGraph_t(const vector<Transcript_t>& transcripts)
{
	string gene_id = transcripts[0].GeneID;
	uint8_t chr = transcripts[0].Chr;
	bool strand = transcripts[0].Strand;
	// make sure all transcripts in the vector correspond to the same gene
	for (int32_t i = 1; i < transcripts.size(); i++) {
		assert(transcripts[i].GeneID == gene_id);
		assert(transcripts[i].Chr == chr);
		assert(transcripts[i].Strand == strand);
	}
	// record gene ID
	GeneID = gene_id;
	// construct non-overlapping unique intervals from exons
	// intervals can only be split, not merged
	vector<Exon_t> nonoverlap_exons;
	vector< pair<int32_t,bool> > boundaries; // pair<position, is_end>
	for (const Transcript_t& t : transcripts) {
		for (const Exon_t& e : t.vExons) {
			boundaries.push_back( make_pair(e.StartPos, false) );
			boundaries.push_back( make_pair(e.EndPos, true) );
		}
	}
	sort(boundaries.begin(), boundaries.end(), 
		[](pair<int32_t,bool> a, pair<int32_t,bool> b){if (a.first != b.first) return a.first < b.first; else return a.second < b.second;} );
	int32_t counter = 1;
	int32_t last_position = boundaries[0].first;
	assert(boundaries[0].second == false);
	for (int32_t i = 1; i < boundaries.size(); i++) {
		if (boundaries[i].second == false) {
			if (counter > 0 && last_position != boundaries[i].first) {
				Exon_t e(last_position, boundaries[i].first);
				nonoverlap_exons.push_back(e);
			}
			last_position = boundaries[i].first;
			counter ++;
		}
		else {
			assert(counter > 0);
			if (last_position != boundaries[i].first) {
				Exon_t e(last_position, boundaries[i].first);
				nonoverlap_exons.push_back(e);
			}
			last_position = boundaries[i].first;
			counter --;
		}
	}
	// construct graph node from non-overlapping unique intervals
	for (const Exon_t& e : nonoverlap_exons) {
		Node_t node(chr, e.StartPos, e.EndPos, strand, true);
		vNodes.push_back(node);
	}
	// adding source and sink node
	{
		Node_t node(chr, -1, -1, strand, false);
		vNodes.push_back(node);
	}
	{
		Node_t node(chr, std::numeric_limits<int32_t>::max(), std::numeric_limits<int32_t>::max(), strand, false);
		vNodes.push_back(node);
	}
	sort(vNodes.begin(), vNodes.end());
	if (!strand)
		reverse(vNodes.begin(), vNodes.end());
	for (int32_t i = 0; i < vNodes.size(); i++) {
		vNodes[i].setID(i);
		vNodes[i].CalculateLength();
	}
	// adding edges that correspond to known transcripts
	for (int32_t i = 0; i < transcripts.size(); i++) {
		const Transcript_t& t = transcripts[i];
		int32_t idx_t = 0;
		int32_t idx_node = 0;
		vector<int32_t> connected_node_indices;
		connected_node_indices.push_back(0);
		while (idx_node < vNodes.size() && idx_t < t.vExons.size()) {
			// vNodes must be either totally not overlap with exons in t, or totally within
			assert(vNodes[idx_node].EndPos <= t.vExons[idx_t].StartPos || t.vExons[idx_t].EndPos <= vNodes[idx_node].StartPos 
				|| (vNodes[idx_node].StartPos >= t.vExons[idx_t].StartPos && vNodes[idx_node].EndPos <= t.vExons[idx_t].EndPos));
			if (strand && vNodes[idx_node].EndPos <= t.vExons[idx_t].StartPos)
				idx_node ++;
			else if (!strand && vNodes[idx_node].StartPos >= t.vExons[idx_t].EndPos)
				idx_node ++;
			else if (strand && t.vExons[idx_t].EndPos <= vNodes[idx_node].StartPos)
				idx_t ++;
			else if (!strand && t.vExons[idx_t].StartPos >= vNodes[idx_node].EndPos)
				idx_t ++;
			else {
				connected_node_indices.push_back(idx_node);
				idx_node ++;
			}
		}
		connected_node_indices.push_back(vNodes.size() - 1);
		// connecting nodes in order
		assert(is_sorted(connected_node_indices.begin(), connected_node_indices.end()));
		for (int32_t j = 1; j < connected_node_indices.size(); j++) {
			Edge_t edge(connected_node_indices[j-1], connected_node_indices[j]);
			edge.IncidentTranscripts.push_back(t.TransID);
			vEdges.push_back(edge);
		}
	}
	sort(vEdges.begin(), vEdges.end());
	// merge same edges
	vector<Edge_t> newvEdges;
	for (const Edge_t& e : vEdges) {
		if (newvEdges.size() == 0 || !(e == newvEdges.back()))
			newvEdges.push_back(e);
		else
			newvEdges.back().IncidentTranscripts.insert( newvEdges.back().IncidentTranscripts.end(), e.IncidentTranscripts.cbegin(), e.IncidentTranscripts.cend() );
	}
	vEdges = newvEdges;
	vEdges.shrink_to_fit();
	for (int32_t i = 0; i < vEdges.size(); i++) {
		vector<string>& incitrans = vEdges[i].IncidentTranscripts;
		sort(incitrans.begin(), incitrans.end());
		incitrans.resize( distance(incitrans.begin(), unique(incitrans.begin(), incitrans.end())) );
		vEdges[i].setID(i);
	}
	UpdateNodeEdgeReference();
};


void GeneGraph_t::UpdateNodeEdgeReference()
{
	for (int32_t i = 0; i < vNodes.size(); i++) {
		vNodes[i].InEdges.clear();
		vNodes[i].OutEdges.clear();
	}
	for (const Edge_t& e : vEdges) {
		vNodes[e.Node_1].OutEdges.push_back(e.ID);
		vNodes[e.Node_2].InEdges.push_back(e.ID);
	}
};

void GeneGraph_t::CalculateBiasMultiplier(const vector<Transcript_t>& transcripts, const vector< vector<double> >& biascorrections, const vector<double>& expr)
{
	// initialize BiasMultiplier to 0
	for (Node_t& n : vNodes)
		n.BiasMultiplier = 0;
	// calculate the total expression
	vector<double> expr_node(vNodes.size(), 0);
	for (int32_t i = 0; i < transcripts.size(); i++) {
		const Transcript_t& t = transcripts[i];
		const vector<double>& bias = biascorrections[i];
		assert(vNodes[0].Strand == transcripts[0].Strand);
		bool strand = vNodes[0].Strand;
		int32_t covered_length = 0;
		int32_t idx_t = 0;
		int32_t idx_node = 0;
		while (idx_node < vNodes.size() && idx_t < t.vExons.size()) {
			// vNodes must be either totally not overlap with exons in t, or totally within
			assert(vNodes[idx_node].EndPos <= t.vExons[idx_t].StartPos || t.vExons[idx_t].EndPos <= vNodes[idx_node].StartPos 
				|| (vNodes[idx_node].StartPos >= t.vExons[idx_t].StartPos && vNodes[idx_node].EndPos <= t.vExons[idx_t].EndPos));
			if (strand && vNodes[idx_node].EndPos <= t.vExons[idx_t].StartPos)
				idx_node ++;
			else if (!strand && vNodes[idx_node].StartPos >= t.vExons[idx_t].EndPos)
				idx_node ++;
			else if (strand && t.vExons[idx_t].EndPos <= vNodes[idx_node].StartPos) {
				covered_length += t.vExons[idx_t].EndPos - t.vExons[idx_t].StartPos;
				idx_t ++;
			}
			else if (!strand && t.vExons[idx_t].StartPos >= vNodes[idx_node].EndPos) {
				covered_length += t.vExons[idx_t].EndPos - t.vExons[idx_t].StartPos;
				idx_t ++;
			}
			else {
				if (covered_length + vNodes[idx_node].EndPos - vNodes[idx_node].StartPos > bias.size())
					cout << "watch here\n";
				assert(covered_length + vNodes[idx_node].EndPos - vNodes[idx_node].StartPos <= bias.size());
				if (strand) {
					for (int32_t pos = vNodes[idx_node].StartPos - t.vExons[idx_t].StartPos; pos < vNodes[idx_node].EndPos - t.vExons[idx_t].StartPos; pos++)
						vNodes[idx_node].BiasMultiplier += bias[covered_length + pos] * expr[i];
					expr_node[idx_node] += expr[i];
				}
				else {
					for (int32_t pos = t.vExons[idx_t].EndPos - vNodes[idx_node].EndPos; pos < t.vExons[idx_t].EndPos - vNodes[idx_node].StartPos; pos++)
						vNodes[idx_node].BiasMultiplier += bias[covered_length + pos] * expr[i];
					expr_node[idx_node] += expr[i];
				}
				idx_node ++;
			}
		}
	}
	// convert the effective length to an actual multiplier
	for (int32_t i = 1; i+1 < vNodes.size(); i++) {
		assert(vNodes[i].BiasMultiplier >= 0);
		vNodes[i].BiasMultiplier /= (expr_node[i] * vNodes[i].Length);
	}
	// add BiasMultiplier for the source and sink
	vNodes[0].BiasMultiplier = 0;
	vNodes.back().BiasMultiplier = 0;
};


// The following code reconnect by looking at the EXISTING TRANSCRIPT, which is NOT consistent with what we want.
void GeneGraph_t::RemoveZeroProbNodes()
{
	// assume that BiasMultiplier has already been calculated, check on this assumption
	double sum_bias_multiplier = 0;
	for (const Node_t& n : vNodes) {
		assert(n.BiasMultiplier >= 0);
		sum_bias_multiplier += n.BiasMultiplier;
	}
	assert(sum_bias_multiplier > 0);
	// search for zero-BiasMultiplier nodes, and reconnect the in-nodes to out-nodes by looking at the transcripts
	// first create a new vector of nodes to record the subset of nodes
	// at the same time, create new edges with respect to the old node indexes, and merge duplicated edges as usual
	vector<Node_t> newvNodes;
	vector<int32_t> ID_removed;
	for (const Node_t& n : vNodes) {
		if (n.StartPos == -1 || n.StartPos == std::numeric_limits<int32_t>::max() || n.BiasMultiplier > 0)
			newvNodes.push_back(n);
		else {
			// record the node ID to be removed
			ID_removed.push_back( n.ID );
		}
	}
	// if no nodes are removed, then don't need to care about all the following reconnecting / node index mapping
	if (ID_removed.size() == 0)
		return;
	// reconnect disjoint nodes after removal: collect paths corresponding to incident transcripts, reconnect two sides of the removal ones
	map< string,vector<int32_t> > TransPaths;
	for (const Edge_t& e : vEdges) {
		for (int32_t i = 0; i < e.IncidentTranscripts.size(); i++) {
			string s = e.IncidentTranscripts[i];
			map< string,vector<int32_t> >::iterator itpath = TransPaths.find(s);
			if (itpath != TransPaths.end())
				(itpath->second).push_back( e.ID );
			else {
				vector<int32_t> tmp = {e.ID};
				TransPaths[s] = tmp;
			}
		}
	}
	for (map< string,vector<int32_t> >::const_iterator itpath = TransPaths.cbegin(); itpath != TransPaths.cend(); itpath ++) {
		int32_t idx_pathelement = 0;
		int32_t in_node = -1;
		while (idx_pathelement < itpath->second.size()) {
			int32_t idx_edge = itpath->second[idx_pathelement];
			// make assertion that the previous edge Node_2 == this edge Node_1
			assert(idx_pathelement == 0 || vEdges[itpath->second[idx_pathelement-1]].Node_2 == vEdges[idx_edge].Node_1);
			if (binary_search(ID_removed.begin(), ID_removed.end(), vEdges[idx_edge].Node_2) && !binary_search(ID_removed.begin(), ID_removed.end(), vEdges[idx_edge].Node_1)) {
				assert(in_node == -1);
				in_node = vEdges[idx_edge].Node_1;
			}
			else if (binary_search(ID_removed.begin(), ID_removed.end(), vEdges[idx_edge].Node_1) && !binary_search(ID_removed.begin(), ID_removed.end(), vEdges[idx_edge].Node_2)) {
				assert(in_node >= 0 and in_node < vNodes.size());
				int32_t out_node = vEdges[idx_edge].Node_2;
				// add new connecting edges
				Edge_t tmp(in_node, out_node);
				tmp.IncidentTranscripts.clear();
				tmp.IncidentTranscripts.push_back( itpath->first );
				vEdges.push_back(tmp);
				// reset in_node
				in_node = -1;
			}
			idx_pathelement ++;
		}
	}
	// update nodes, update and merge edges
	{
		vNodes = newvNodes;
		vNodes.shrink_to_fit();
		sort(vEdges.begin(), vEdges.end());
		vector<Edge_t> newvEdges;
		for (const Edge_t& e : vEdges) {
			if (newvEdges.size() == 0 || !(e == newvEdges.back()))
				newvEdges.push_back(e);
			else
				newvEdges.back().IncidentTranscripts.insert( newvEdges.back().IncidentTranscripts.end(), e.IncidentTranscripts.cbegin(), e.IncidentTranscripts.cend() );
		}
		for (Edge_t& e : newvEdges) {
			vector<string>& incitrans = e.IncidentTranscripts;
			sort(incitrans.begin(), incitrans.end());
			incitrans.resize( distance(incitrans.begin(), unique(incitrans.begin(),incitrans.end())) );
		}
		vEdges = newvEdges;
	}

	// remove edges related to the removed node, and merge edges
	sort(vEdges.begin(), vEdges.end());
	vector<Edge_t> newvEdges;
	for (const Edge_t& e : vEdges) {
		if (binary_search(ID_removed.begin(), ID_removed.end(), e.Node_1) || binary_search(ID_removed.begin(), ID_removed.end(), e.Node_2))
			continue;
		if (newvEdges.size() == 0 || !(e == newvEdges.back()))
			newvEdges.push_back(e);
		else {
			vector<string>& incitrans = newvEdges.back().IncidentTranscripts;
			incitrans.insert(incitrans.end(), e.IncidentTranscripts.begin(), e.IncidentTranscripts.end());
		}
	}
	for (int32_t i = 0; i < newvEdges.size(); i++)
		newvEdges[i].setID(i);
	vEdges = newvEdges;
	vEdges.shrink_to_fit();
	// mapping the old node index to the new node index
	map<int32_t,int32_t> NodeIDMap;
	for (int32_t i = 0; i < vNodes.size(); i++)
		NodeIDMap[vNodes[i].ID] = i;
	// updating the node index of the edges
	for (Edge_t& e : vEdges) {
		map<int32_t,int32_t>::const_iterator it1 = NodeIDMap.find(e.Node_1);
		map<int32_t,int32_t>::const_iterator it2 = NodeIDMap.find(e.Node_2);
		assert(it1 != NodeIDMap.cend() && it2 != NodeIDMap.cend());
		e.Node_1 = it1->second;
		e.Node_2 = it2->second;
	}
	// update node ID
	for (int32_t i = 0; i < vNodes.size(); i++)
		vNodes[i].setID(i);
	// update the edge pointer in nodes
	UpdateNodeEdgeReference();

	// sanity check: the transcript paths should be proper
	TransPaths.clear();
	for (const Edge_t& e : vEdges) {
		for (int32_t i = 0; i < e.IncidentTranscripts.size(); i++) {
			string s = e.IncidentTranscripts[i];
			map< string,vector<int32_t> >::iterator itpath = TransPaths.find(s);
			if (itpath != TransPaths.end())
				(itpath->second).push_back( e.ID );
			else {
				vector<int32_t> tmp = {e.ID};
				TransPaths[s] = tmp;
			}
		}
	}
	for (map< string,vector<int32_t> >::const_iterator itpath = TransPaths.cbegin(); itpath != TransPaths.cend(); itpath ++) {
		for (int32_t idx_pathelement = 1; idx_pathelement < itpath->second.size(); idx_pathelement++)
			assert(vEdges[itpath->second[idx_pathelement-1]].Node_2 == vEdges[itpath->second[idx_pathelement]].Node_1);
	}
};


string GeneGraph_t::PrintContent(const vector<string>& RefName, const vector<Transcript_t>& transcripts) const
{
	string s = "";
	for (const Node_t& n : vNodes) {
		s += "node\t" + to_string(n.ID) +"\t"+ RefName[n.Chr] +"\t"+ to_string(n.StartPos) +"\t"+ to_string(n.EndPos) +"\t"+ to_string(n.Strand) +"\t";
		s += to_string(n.Coverage) +"\t"+ to_string(n.BiasMultiplier) +"\t";
		// InEdges
		string s_inedge = "";
		for (const int32_t& inedge : n.InEdges)
			s_inedge += to_string(inedge) +",";
		s_inedge = s_inedge.substr(0, s_inedge.size() - 1);
		s += "(" + s_inedge + ")\t";
		// OutEdges
		string s_outedge = "";
		for (const int32_t& outedge : n.OutEdges)
			s_outedge += to_string(outedge) +",";
		s_outedge = s_outedge.substr(0, s_outedge.size() - 1);
		s += "(" + s_outedge + ")\n";
	}
	for (const Edge_t& e : vEdges) {
		s += "edge\t" + to_string(e.ID) +"\t"+ to_string(e.Node_1) +"\t"+ to_string(e.Node_2) +"\t(";
		string s_incitrans = "";
		for (int32_t i = 0; i < e.IncidentTranscripts.size(); i++)
			s_incitrans += e.IncidentTranscripts[i] +",";
		s_incitrans = s_incitrans.substr(0, s_incitrans.size() - 1);
		s += s_incitrans +")\n";
	}
	return s;
};


void ConstructAllGeneGraphs(const vector<Transcript_t>& transcripts, const map< string, vector<int32_t> >& GeneTransMap, vector<GeneGraph_t>& GeneGraphs)
{
	time_t CurrentTime;
	string CurrentTimeStr;
	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout<<"["<<CurrentTimeStr.substr(0, CurrentTimeStr.size()-1)<<"] "<<"Constructing gene graph."<<endl;

	GeneGraphs.clear();
	for (map< string, vector<int32_t> >::const_iterator it = GeneTransMap.cbegin(); it != GeneTransMap.cend(); it++) {
		vector<Transcript_t> tmp_transcripts;
		for (const int32_t& i : it->second)
			tmp_transcripts.push_back(transcripts[i]);
		GeneGraph_t gg(tmp_transcripts);
		GeneGraphs.push_back(gg);
	}

	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout << "[" << CurrentTimeStr.substr(0, CurrentTimeStr.size()-1) << "] " << "Finish. Number of graphs = " << GeneGraphs.size() <<"\n";
};


void UpdateAllBiasMultiplier(const vector<Transcript_t>& transcripts, const vector< vector<double> >& corrections, const vector<double>& expr,
	const map< string, vector<int32_t> >& GeneTransMap, vector<GeneGraph_t>& GeneGraphs)
{
	time_t CurrentTime;
	string CurrentTimeStr;
	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout<<"["<<CurrentTimeStr.substr(0, CurrentTimeStr.size()-1)<<"] "<<"Updating bias correction for gene graph."<<endl;

	for (vector<GeneGraph_t>::iterator it = GeneGraphs.begin(); it != GeneGraphs.end(); it++) {
		// find corresponding transcript index
		map< string,vector<int32_t> >::const_iterator ittrans = GeneTransMap.find(it->GeneID);
		assert(ittrans != GeneTransMap.cend());
		// collect corresponding transcripts
		vector<Transcript_t> tmp_transcripts;
		for (const int32_t& i : ittrans->second)
			tmp_transcripts.push_back( transcripts[i] );
		// collect the bias correction vector of corresponding transcripts
		vector< vector<double> > tmp_corrections;
		for (const int32_t& i : ittrans->second)
			tmp_corrections.push_back( corrections[i] );
		// collect corresponding salmon quantification TPM
		vector<double> tmp_expr;
		for (const int32_t& i : ittrans->second)
			tmp_expr.push_back( expr[i] );
		it->CalculateBiasMultiplier( tmp_transcripts, tmp_corrections, tmp_expr );
	}

	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout << "[" << CurrentTimeStr.substr(0, CurrentTimeStr.size()-1) << "] " << "Finish.\n";
};


void WriteGraph(string outputfile, const vector<GeneGraph_t>& GeneGraphs, const vector<string>& RefName, const vector<Transcript_t>& transcripts)
{
	ofstream output(outputfile, ios::out);
	output << "# graph\tGeneID\tNumNodes\tNumEdges\n";
	output << "# node\tID\tChr\tStartPos\tEndPos\tStrand\tCoverage\tBiasMultiplier\tInEdges\tOutEdges\n";
	output << "# edge\tID\tNode_1\tNode_2\tIncidentTranscripts\n";
	for (const GeneGraph_t& g : GeneGraphs) {
		output <<"graph\t"<< g.GeneID <<"\t"<< (g.vNodes.size()) <<"\t"<< (g.vEdges.size()) << endl;
		string s = g.PrintContent(RefName, transcripts);
		output << s;
	}
	output.close();
};

