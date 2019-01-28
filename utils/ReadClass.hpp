#ifndef __ReadClass_T__
#define __ReadClass_T__

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <string>
#include "htslib/sam.h"


class GeneGraph_t;


// a SingleAlignment_t means one alignment in multi-mapping. It contains all aligned pieces of the whole fragment.
// the alignment positions of first and second mate are merged: the unique nucleotides over the aligned pieces of the whole fragment
class SingleAlignment_t
{
public:
	int32_t GraphID;
	vector<int32_t> NodeIDs;
	vector< pair<int32_t,int32_t> > AlignedPieces; // aligned pieces should be in genome coordinate

public:
	SingleAlignment_t(){};

	bool Smaller_node (const SingleAlignment_t& rhs) const {
		if (GraphID != rhs.GraphID)
			return GraphID < rhs.GraphID;
		else {
			// compare each node
			for (int32_t i = 0; i < min(NodeIDs.size(), rhs.NodeIDs.size()); i++) {
				if (NodeIDs[i] != rhs.NodeIDs[i])
					return NodeIDs[i] < rhs.NodeIDs[i];
			}
			// if one has more nodes than the other
			if (NodeIDs.size() != rhs.NodeIDs.size())
				return NodeIDs.size() < rhs.NodeIDs.size();
			return false;
		}
	};

	bool Equal_node (const SingleAlignment_t& rhs) const {
		if (GraphID != rhs.GraphID)
			return false;
		else {
			// compare each node
			for (int32_t i = 0; i < min(NodeIDs.size(), rhs.NodeIDs.size()); i++) {
				if (NodeIDs[i] != rhs.NodeIDs[i])
					return false;
			}
			// if one has more nodes than the other
			if (NodeIDs.size() != rhs.NodeIDs.size())
				return false;
			return true;
		}
	};

	inline void SortPieces(bool graph_strand){
		assert(NodeIDs.size() == AlignedPieces.size());
		// index of nodes / aligned pieces
		vector<int32_t> indexes(NodeIDs.size(), 0);
		std::iota(indexes.begin(), indexes.end(), 0);
		sort(indexes.begin(), indexes.end(), [&AlignedPieces](int32_t a, int32_t b)
			{if (AlignedPieces[a].first != AlignedPieces[b].first) return AlignedPieces[a].first < AlignedPieces[b].first; 
				else return AlignedPieces[a].second < AlignedPieces[b].second;} );
		if (!graph_strand)
			reverse(indexes.begin(), reverse.end());
		// sort NodeIDs and AlignedPieces
		vector<int32_t> newNodeIDs;
		vector< pair<int32_t,int32_t> > newAlignedPieces;
		for (int32_t i = 0; i < indexes.size(); i++) {
			newNodeIDs.push_back( NodeIDs[indexes[i]] );
			newAlignedPieces.push_back( newAlignedPieces[indexes[i]] );
		}
		NodeIDs = newNodeIDs;
		AlignedPieces = newAlignedPieces;
		// after sorting according to aligned pieces, the nodes should be also in order
		assert(is_sorted(NodeIDs.begin(), NodeIDs.end()));
	};
};


class MultiAlignment_t
{
public:
	string ReadName;
	vector<SingleAlignment_t> Alignments; 

public:
	MultiAlignment_t(){};
	MultiAlignment_t()
};


class FragmentEquivalence_t
{};


#endif