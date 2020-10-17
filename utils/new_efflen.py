import sys
from trie_conversion import PrefixTrie, RType, load_data
import struct
from collections import defaultdict
import time
from TranscriptClass import ReadGTF, Map_Gene_Trans
import pickle

def module_desc():
    '''
    This module handles effective length calculation and related matters.
    The way effective length is calculated in this module is consistent with
    GS paper Sect. 3.4.
    '''
    pass

class Tick:
    def __init__(self, name = ""):
        self._start = time.perf_counter()
        self._name = name
        self._last = self._start

    def __enter__(self):
        return self

    def tick(self, *desc):
        print(" |{:.2f}(+{:.2f})s {} - Tick".format(time.perf_counter() - self._start,
            time.perf_counter() - self._last,
            self._name), *desc)
        self._last = time.perf_counter()

    def __call__(self, *desc):
        return self.tick(*desc)


    def __exit__(self, type, value, traceback):
        print("{:.2f}s {} - Finished".format(time.perf_counter() - self._start,
                                             self._name))

def load_transcripts(path_fn, gtf_fn):
    '''
    Load all transcripts and get all paths and genomic coordinates.
    @param path_fn: File containing gene paths and exon coordinates.
    @param gtf_fn: The standard GTF annotation file, contains transcript<->gene
                    graph mapping.
    @return: path, coord. first one is {graph ID -> list of (transcript ID, path)},
                          second one is {graph ID -> dict of {exon ID -> coordinate}}.
    '''
    #  path_file = "/home/hongyuz/data/simu_Full_GEU_new_sequence_header.fa"
    #  Transcripts = ReadGTF("/home/congm1/savanna/savannacong33/NCBI/gencode.v26.annotation.gtf")
    transcripts = ReadGTF(gtf_fn)
    GeneTransMap, TransGeneMap = Map_Gene_Trans(transcripts)
    ret_path = dict()
    ret_coord = dict()
    with open(path_fn) as f:
        for line in f:
            tr, lp = line.split()
            tr = tr[1:]
            # Get ID of corresponding gene graph
            if "_simu" in tr:
                gname = tr[:tr.find("_sim")]
            else:
                gname = TransGeneMap[tr]
            if gname not in ret_path:
                ret_path[gname] = []
                ret_coord[gname] = dict()
            # Parse input
            path = []
            _chr = ""
            for s in lp[1:-1].split(')('):
                tmp = s.split(',')
                exon = int(tmp[0])
                path.append(exon)
                if _chr == "":
                    _chr = tmp[1]
                else:
                    assert _chr == tmp[1]
                cc = (int(tmp[2]), int(tmp[3]))
                if exon in ret_coord[gname]:
                    #  print(gname, exon, ret_coord[gname][exon], cc)
                    assert ret_coord[gname][exon] == cc
                else:
                    ret_coord[gname][exon] = cc
            ret_path[gname].append((tr, path))

    return ret_path, ret_coord

def load_salmon_quant(fn):
    '''
    Load quantification result of Salmon.
    @param fn: Filename for quant.sf.
    @return: dictionary: transcript ID -> copy number.
    '''
    ret = dict()
    with open(fn) as f:
        f.readline()
        for line in f:
            s = line.split()
            ret[s[0]] = float(s[3])
    return ret

def load_matrix_likelihood(fn):
    '''
    Load likelihood estimates from matrices.
    @param fn: Filename.
    @return dictionary: transcript ID -> triplet of (start, end, likelihood).
    '''
    from utils import ReadCorrection_matrix
    with Tick("Load Matrices from file"):
        mat = ReadCorrection_matrix(fn)
    ret = dict()
    with Tick("Parse Matrices") as tick:
        for gid in mat:
            cur = mat[gid]
            ll = []
            for x in range(cur.shape[0]):
                for y in range(cur.shape[1]):
                    if cur[x][y] > 0:
                        assert x > y
                        ll.append((y, x, cur[x][y]))
            ret[gid] = ll
            #  tick.tick(gid, cur.shape)
    return ret

class EfflenCalculator:
    '''
    This is a wrapper over real calculator and it's basically one function.
    The reason this is a class is to do introspections.

    Calculate Effective Length given the following: flow graph specification(
    including all fragments), gene coordinate specification and transcript
    specification(via path_fn and gtf_fn), copy numbers and likelihood matrices.
    @param graph_data_fn: data file generated by trie_conversion.
    @param path_fn, gtf_fn: external data specification.
    @param cns: dictionary of {Transcript ID: Copy number}.
    @param lls: dictionary of {Transcript ID: List of (start, end, likelihood)}.
    @param ll_gen: Generator of lls.items(). To provide a streaming interface.
    @return: dictionary of {Gene graph ID: List of efflen on flow graph edges}.
    '''

    def __init__(self, graph_data_fn, path_fn, gtf_fn, cns, ll_gen):
        with Tick("Initialization"):
            graphs, eq_classes = load_data(graph_data_fn)
            paths, coords = load_transcripts(path_fn, gtf_fn)
            # Bookkeeping (pickling uses)
            self._paths = paths
            self._cns = cns
        with Tick("Efflen Generation"):
            self.flow_efflen, self.gene_effgroup, self.gene_efflen, self.discards = self._work(
            graphs, paths, coords, cns, ll_gen)
        print(self._count_discards())

    def _work(self, s_graphs, s_paths, s_coords, s_cns, gen_lls):
        '''
        Main working routine.
        '''
        self._cn_group = dict()
        #  self._ll_group = dict()
        # BUild the reverse lookup table
        rev_lookup = dict()
        for gid in s_paths:
            for tid, path in s_paths[gid]:
                rev_lookup[tid] = (gid, path)
        d_cns = dict()
        d_lls = dict()
        # Build the transcript->likelihood dict
        for tid, ll in gen_lls:
            gid, path = rev_lookup[tid]
            if gid not in d_lls:
                d_cns[gid] = defaultdict(float)
                d_lls[gid] = defaultdict(dict)
            if tid not in s_cns:
                print("WARNING - Key not found in abundance file:", tid)
                d_cns[gid][tid] = 0
            else:
                d_cns[gid][tid] = s_cns[tid]
            d_lls[gid][tid] = self._unpack_transcript(ll, path, s_coords[gid])
        self._cn_group = d_cns
        flow_efflen = dict()
        gene_efflen = dict()
        gene_effgroup = dict()
        discard_paths = dict()
        for gid in s_graphs:
            graph = s_graphs[gid]
            #  eq_classes = self.eq_classes[gid]
            paths = s_paths[gid]
            coords = s_coords[gid]
            copy_nums = d_cns.get(gid, dict())
            likelihoods = d_lls.get(gid, dict())
            _f, _gg, _g, _d = self._work_graph(graph, paths, coords,
                                                  copy_nums, likelihoods)
            flow_efflen[gid] = _f
            gene_effgroup[gid] = _gg
            gene_efflen[gid] = _g
            discard_paths[gid] = _d
        return flow_efflen, gene_effgroup, gene_efflen, discard_paths

    def _unpack_transcript(self, ll, path, coords):
        '''
        Convert a transcript likelihood to L_p.
        @param ll: the list of (start, end, likelihood).
        @param path: the gene path.
        @param coords: the coordinates.
        @return: a dictionary of {subpaths : likelihood}.
        '''
        coord_map = []  # Position -> Exon (Backref to Path) for this transcript
        ret = defaultdict(float)
        for i in range(len(path)):
            for _ in range(coords[path[i]][1] - coords[path[i]][0]):
                coord_map.append(i)
        for start, end, likelihood in ll:
            gpath = tuple(path[coord_map[start]:coord_map[end - 1] + 1])
            ret[gpath] += likelihood
        return ret


    def _work_graph(self, graph, paths, coords, cns, lls):
        '''
        Figure this mess out for each gene graph.
        @param ...: look into __init__. eq_classes is not explicitly used here.
        @return: efflen for flow graph edge, efflen for gene graph paths, list
                    of discarded paths due to not matching anything.
        '''
        # Regroup L_p for different transcripts
        gpath_groups = dict()
        for tid, ll in lls.items():
            for gpath, val in ll.items():
                if gpath not in gpath_groups:
                    gpath_groups[gpath] = defaultdict(float)
                gpath_groups[gpath][tid] = val
        # Distribute c_p to F_e
        discarded_weights = []
        fpath_weight = defaultdict(float)
        gpath_weight = dict()
        for gpath, vals in gpath_groups.items():
            # zero out check
            total_cp = sum(cns[x] for x in gpath_groups[gpath])
            if total_cp < 1e-10:
                val = sum(vals.values()) / len(vals)
            else:
                val = sum(cns[key] * val for key, val in \
                          gpath_groups[gpath].items()) / total_cp
            gpath_weight[gpath] = val
            elist = graph.match(list(gpath))
            if len(elist) == 0:
                discarded_weights.append((gpath, val))
            else:
                for e in elist:
                    fpath_weight[e] += val
        return fpath_weight, gpath_groups, gpath_weight, discarded_weights

    def _count_discards(self, threshold = 0.1):
        '''
        Print a summary of discarded likelihoods for each gene.
        '''
        s = ""
        all_discarded = 0
        all_used = 0
        for gid in self.flow_efflen:
            used = sum(self.gene_efflen[gid].values())
            discarded = sum(v[1] for v in self.discards[gid])
            all_discarded += discarded
            all_used += used
            if discarded + used < 1e-15:
                continue
            ratio = discarded / (discarded + used)
            if ratio > threshold:
                s += "{}: {:.4f}/{:.4f} {:.2f}% | ".format(gid, discarded, used, 100 * ratio)
        s += "[Total] {:.2f} / {:.2f} {:.2f}%".format(all_discarded, all_used,
                                100 * all_discarded / (all_discarded + all_used))
        return s

def unpack_sparse_matrix(fn):
    '''
    Improved version of sparse matrix unpacking.
    Note: now yields values instead of returning a dictionary
    '''
    with Tick("Unpack Sparse Likelihood") as tick:
        with open(fn, "rb") as f:
            tc = struct.unpack('i', f.read(4))[0]
            for dx in range(tc):
                nlen = struct.unpack('i', f.read(4))[0]
                slen = struct.unpack('i', f.read(4))[0]
                name = f.read(nlen).decode('utf-8')
                data_stream = f.read(16 * slen)
                ll = []
                for i in range(slen):
                    offset = 16 * i
                    _start = struct.unpack('i', data_stream[offset:offset+4])[0]
                    _end = struct.unpack('i', data_stream[offset+4:offset+8])[0]
                    _val = struct.unpack('d', data_stream[offset+8:offset+16])[0]
                    ll.append((_start, _end, _val))
                    assert _start < _end
                if dx % 50 == 0:
                    tick.tick(name, slen, dx, '/', tc)
                    #  if dx == 500:
                        #  break
                yield name, ll


if __name__ == "__main__":
    path_file = "/home/hongyuz/data/simu_Full_GEU_new_sequence_header.fa"
    gtf_file = sys.argv[1]
    graph_data_file = sys.argv[2]
    salmon_result_file = sys.argv[3]
    #  matrix_file = "/mnt/disk33/user/congm1/small_graphsalmon/simu01/salmon/gs_corrections_matrix_100.dat"
    sparse_matrix_file = sys.argv[4]
    outfile = sys.argv[5]

    with Tick("Load Salmon Quant"):
        cns = load_salmon_quant(salmon_result_file)
    #  lls = load_matrix_likelihood(matrix_file)
    ll_gen = unpack_sparse_matrix(sparse_matrix_file)
    ec = EfflenCalculator(graph_data_file, path_file, gtf_file, cns, ll_gen)
    with Tick("Pickling Data"):
        #with open("/home/hongyuz/data/5Kgenes.efflen.dat", "wb") as f1:
        #    pickle.dump(ec.flow_efflen, f1)
        with open(outfile, "wb") as f2:
            pickle.dump(ec, f2)
    #  from IPython import embed
    #  embed()
