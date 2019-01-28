import time
import io
import collections
import pickle
def module_desc():
    '''
    This module transforms input gene graph and hyperedge/equivalance classes
    into new graph based on prefix trie.
    Refer to the design documents for reasoning and detailed specification of
    the model and a bit of discussion.
    '''
    pass

class PrefixTrie:
    '''
    A proper trie that runs slowly but okay - optimizations to come later.
    self.nodes: vID -> List of vertex (representing hyperedge).
    self.node_lookup: Tuple of vertex -> vID.
    self.edges: eID -> (S, T).
    self.edge_lookup: (S, T) -> eID.
    self.out_dest: S -> (target exon -> T).
    '''
    def __init__(self, _edges, hedges, name = ""):
        '''
        Initialization routine that prunes unreachable vertices and edges.
        @param _edges: list of (normal) edges. This dictates possible transition.
        @param hedges: list of hyperedges. This dictates possible states.
        @param name: Name of this graph, for bookkeeping purposes.
        '''
        # Preprocessing
        self._edges = _edges
        self._hedges = hedges
        old_n = max(d[1] for d in _edges) + 1
        hpath = list([x] for x in range(old_n))
        for l in hedges:
            if list(l)[:-1] not in hpath:
                hpath.append(list(l)[:-1])
        self._name = name
        self.nodes = []
        self.node_lookup = dict()
        self.n = 0
        # Build the initial node list
        for p in hpath:
            for i in range(1, len(p) + 1): # Iterate over all prefixes
                kw = tuple(p[:i])
                if kw not in self.node_lookup:
                    self.nodes.append(p[:i])
                    self.node_lookup[kw] = self.n
                    self.n += 1
        # Build the helper out-table
        self._uouts = collections.defaultdict(list)
        for u, v in _edges:
            # check duplicate edges
            assert v not in self._uouts[u]
            self._uouts[u].append(v)
        # Purge the vertices and rebuild the node list
        self._visited = [False] * len(self.nodes)
        self._traverse(self.node_lookup[(0,)])
        self.nodes = list(self.nodes[c] for c in range(self.n) if self._visited[c])
        self.n = len(self.nodes)
        self.node_lookup = dict()
        for i in range(self.n):
            self.node_lookup[tuple(self.nodes[i])] = i
        # Perform topological sort and shuffle the node list
        self.nodes = self._topo_sort()
        self.node_lookup = dict()
        for i in range(self.n):
            self.node_lookup[tuple(self.nodes[i])] = i
        # Index all edges
        self.edges = []
        self.m = 0
        self.out_dest = dict()
        self.edge_lookup = dict()
        self.edge_by_trans = dict()
        for d in range(old_n):
            self.edge_by_trans[d] = []
        for i in range(self.n):
            l = self.nodes[i]
            self.out_dest[i] = dict()
            for mt in self._uouts[l[-1]]:
                target = self._next(i, mt)
                self.edges.append((i, target))
                self.out_dest[i][mt] = target
                self.edge_lookup[(i, target)] = self.m
                self.edge_by_trans[mt].append(self.m)
                self.m += 1

    def _next(self, node_id, trans):
        '''
        Helper function that calculates next node to move to.
        '''
        lt = self.nodes[node_id][:] + [trans]
        while tuple(lt) not in self.node_lookup:  # stop if node exists
            lt = lt[1:]  # move to next prefix - should be fail link
            assert len(lt) >= 1
        return self.node_lookup[tuple(lt)]

    def _traverse(self, node_id):
        '''
        Helper function that traverses the graph in DFS order.
        '''
        if self._visited[node_id]:
            return
        self._visited[node_id] = True
        for mv in self._uouts[self.nodes[node_id][-1]]:
            self._traverse(self._next(node_id, mv))

    def _topo_sort(self):
        '''
        Sort the nodes in topological order. Return the new list of nodes.
        '''
        in_deg = collections.defaultdict(int)
        for i in range(self.n):
            last_exon = self.nodes[i][-1]
            for mv in self._uouts[last_exon]:
                t = self._next(i, mv)
                in_deg[t] += 1
        assert in_deg[0] == 0
        q = collections.deque()  # Let's be real, Queue.queue is heavy
        q.append(0)
        ret = []
        while len(q) > 0:
            c = q.popleft()
            ret.append(self.nodes[c][:])
            last_exon = self.nodes[c][-1]
            for mv in self._uouts[last_exon]:
                t = self._next(c, mv)
                in_deg[t] -= 1
                if in_deg[t] == 0:
                    q.append(t)
        return ret

    def get_size(self):
        '''
        Report the size of this trie graph.
        @return (# of vertices, # of edges).
        '''
        return self.n, self.m

    def match(self, l):
        '''
        Match given hyperedge with edges in the trie.
        '''
        assert len(l) >= 1
        #  assert l is list
        d = len(l)
        ret = []
        for e in self.edge_by_trans[l[-1]]:  # First match the transition
            s = self.edges[e][0]
            if (d == 1) or (self.nodes[s][-(d-1):] == l[:-1]): # Then match the state
                # when d == 1, last zero nodes don't work anymore but we know
                # it's a certain match
                ret.append(e)
        return ret

    def walk(self, t):
        '''
        Walk a transcript through the new graph. Return list of vertexes and
        edges visited.
        @param t: the transcript as list of vertices(in old graph). Must start
                    with 0 but don't need to end with n-1.
        @return: (list of vertices, list of edges).
        '''
        assert t[0] == 0
        cur = self.node_lookup[(t[0],)]
        r0 = [cur]
        r1 = []
        for j in t[1:]:
            nx = self.out_dest[cur][j]
            r0.append(nx)
            r1.append(self.edge_lookup[(cur, nx)])
            cur = nx
        return r0, r1

    def _summary(self):
        '''
        Provide an elaborate summary of this graph. Only for debugging purposes.
        @return: A formatted string that can be printed out.
        '''
        f = io.StringIO()
        print("Name:", self._name, file=f)
        print("Edges:", self._edges, file=f)
        print("Hyperedges:", self._hedges, file=f)
        for i in range(self.n):
            s = "V{} = {} | ".format(i, self.nodes[i])
            for tr, target in self.out_dest[i].items():
                s += "{}->V{}({}, E{}) ".format(tr, target,
                                               self.nodes[target],
                                               self.edge_lookup[(i, target)])
            print(s, file=f)
        ret = f.getvalue()
        f.close()
        return ret

class MergedTrie:
    '''
    This corresponds to a "megagraph" - a Trie that corresponds to union of
    multiple gene graphs.
    Mostly for sanity check purposes, but can be used for other things.
    '''
    def __init__(self, d_edges, d_hedges, d_pts = None, name = ""):
        self._d_pts = d_pts
        self._d_edges = d_edges
        self._d_hedges = d_hedges
        self.vert_map = dict()
        self.traceback = dict()
        full_edges = []
        full_hedges = []
        self.n = 2
        self.traceback[0] = ('X', 0)
        self.traceback[1] = ('X', 1)
        for gname in d_edges:
            edges = d_edges[gname]
            hedges = d_hedges.get(gname, [])
            tn = max(k[1] for k in edges) + 1
            # Establish vertex mapping
            self.vert_map[gname] = dict()
            for i in range(tn):
                self.vert_map[gname][i] = self.n
                self.traceback[self.n] = (gname, i)
                self.n += 1
            cmap = self.vert_map[gname]
            # New edges - virtual source/sink
            full_edges.append([0, cmap[0]])
            full_edges.append([cmap[tn - 1], 1])
            for _s, _t in edges:
                full_edges.append([cmap[_s], cmap[_t]])
            # New hyperedges - don't need presence of 0s
            for path in hedges:
                full_hedges.append(list(cmap[d] for d in path))
        self.g = PrefixTrie(full_edges, full_hedges)

    @classmethod
    def from_PT(cls, d):
        d_edges = dict()
        d_hedges = dict()
        for k in d:
            d_edges[k] = d[k]._edges
            d_hedges[k] = d[k]._hedges
        return cls(d_edges, d_hedges, d)

    def transform_eq_class(self, orig_eq_classes):
        ret = []
        for ll in orig_eq_classes:
            nc = []
            for item in ll:
                if item.gid in self._d_edges:
                    obj = RType(item.gid, item.vpath, item.epath,
                                self.match(item.gid, item.vpath), item.weights)
                    nc.append(obj)
            if len(nc) > 0:
                ret.append(nc)
        return ret

    def PT2PT_edge_map(self, eid):
        '''
        Map an edge in the merged PrefixTree graph to the edge in invididual
        PrefixTree graph.
        @return (gid, eid).
        special case: If edge is from supersource or to supersink, gid is the
        gene graph id and eid will be -1(source) or -2(sink).
        '''
        v0, v1 = self.edges[eid]
        # map to merged gene graph
        p0 = self.nodes[v0]
        p1 = self.nodes[v1]
        if p0[0] == 0:  # supersource edge
            tmp = p1[0]
            return (self.traceback[tmp][0], -1)

        if p1[0] == 1:  # supersink edge
            tmp = p0[-1]
            return (self.traceback[tmp][0], -2)

        # map to invididual gene graph
        tt = p1[-1]
        tb = self.traceback[tt]
        gid = tb[0]
        trans = tb[1]
        gpath = list(self.traceback[x][1] for x in p0)

        assert self._d_pts is not None
        # map to invididual PT
        id_v0 = self._d_pts[gid].node_lookup[tuple(gpath)]
        id_v1 = self._d_pts[gid].out_dest[id_v0][trans]
        eid = self._d_pts[gid].edge_lookup[(id_v0, id_v1)]

        return (gid, eid)

    @property
    def nodes(self):
        return self.g.nodes

    @property
    def edges(self):
        return self.g.edges

    def translate_fw(self, gname, nodes):
        '''
        Translate between vertices in old graph and vertices in merged graph
        (WARNING - IT'S NOT THE PREFIX TREE!)
        @param gname: Gene graph name.
        @param nodes: List of nodes to translate.
        @return: translated nodes.
        '''
        m = self.vert_map[gname]
        return list(m[c] for c in nodes)

    def match(self, gname, l):
        '''
        Match a path on a particular gene graph to edges in new graph.
        Same as PrefixTrie.match except it also requires gene graph ID.
        '''
        return self.g.match(self.translate_fw(gname, l))

    def walk(self, gname, t):
        '''
        Match a transcript on a particular gene graph to path in new graph.
        Same as PrefixTrie.walk except the new virtual source and virtual sink
        is also added.
        '''
        ll = [0]
        for i in t:
            ll.append(self.vert_map[gname][i])
        ll.append(1)
        return self.g.walk(ll)

    class MegaPT_Wrapper:
        '''
        As if you are using a PrefixTrie.
        '''
        def __init__(self, mt, gid):
            self.mt = mt
            self.gid = gid
            self.nodes = self.mt.nodes
            self.edges = self.mt.edges

        def match(self, l):
            return self.mt.match(self.gid, l)

        def walk(self, l):
            return self.mt.walk(self.gid, l)

    def get_wrappers(self):
        '''
        Get a dictionary (graph_id -> pseudo-PrefixTrie) that wraps over this
        MegaPT.
        '''
        ret = dict()
        for gid in self._d_edges:
            ret[gid] = self.MegaPT_Wrapper(self, gid)
        return ret

    def _summary(self):
        '''
        Print a formatted summary.
        '''
        f = io.StringIO()
        for gid in self._d_edges:
            print("- GID = {}".format(gid), file=f)
            print("- Edges =", self._d_edges[gid], file=f)
            print("- Hyperedges =", self._d_hedges.get(gid, []), file=f)
            print("- Vertex Mapping =", self.vert_map[gid], file=f)
        print(self.g._summary(), file=f)
        ret = f.getvalue()
        f.close()
        return ret

def load_graphs(fn):
    '''
    Load gene graphs from filename.
    @param fn: filename to read graph.
    @return: dictionary of (gene graph name -> list of edges).
    '''
    ret = dict()
    gname = ""
    current = []
    with open(fn) as f:
        for line in f:
            s = line.split()
            if s[0] == "graph":
                if gname != "":
                    ret[gname] = current
                    current = []
                gname = s[1]
            elif s[0] == "edge":
                current.append((int(s[2]), int(s[3])))
    ret[gname] = current
    return ret

def load_eq_classes(fn):
    '''
    Load all Equivalence Classes aka Hyperedges.
    @param fn: filename to eq_classes.
    @return: dictionary of (gene graph name -> list of hyperedges in node #).
    '''
    ret = dict()
    with open(fn) as f:
        f.readline()  # skip header
        for line in f:
            s = line.split()
            lg = s[0][1:-1].split(',')
            ll = s[1][2:-2].split('),(')
            assert len(lg) == len(ll)
            for x in range(len(lg)):
                g = lg[x]
                l = ll[x].split(',')
                l = list(int(x) for x in l)
                if g not in ret:
                    ret[g] = []
                if l not in ret[g]:
                    ret[g].append(l)
    return ret

RType = collections.namedtuple("RType", ["gid", "vpath", "epath",
                                            "new_edges", "weights"])
def translate_eq_classes(tries, fn, as_dict = False, allow_discards = False):
    '''
    Load all Equivalence classes again and translate each of them.
    @param tries: The dictionary of gene graph -> PrefixTrie.
    @param fn: filename for eq_classes.
    @param as_dict: If set to True a dict will be used instead with key:value
                    as mentioned below. This is to eliminate dependency.
    @param allow_discards: If an exception should be raise in case gene graph
                            ID does not exist.
    @return a list of all equivalance classes, in form of list of named tuple.
    The named tuple have the following fields:
        gid: Gene Graph ID.
        vpath: Old path by vertex.
        epath: Old path by edge.
        new_edges: List of equivalent edges in new graph.
        weights: Weights that is copied from original file.
    '''
    ret = []
    with open(fn) as f:
        f.readline()  # skip header
        for line in f:
            s = line.split()
            lg = s[0][1:-1].split(',')
            lv = s[1][2:-2].split('),(')
            le = s[2][2:-2].split('),(')
            lw = s[5][1:-1].split(',')
            assert len(lg) == len(lv)
            cur = []
            for x in range(len(lg)):
                g = lg[x]
                v = list(int(x) for x in lv[x].split(','))
                if le[x] == "":
                    e = []
                else:
                    e = list(int(x) for x in le[x].split(','))
                w = float(lw[x])
                if g not in tries:
                    if allow_discards:
                        continue
                    else:
                        raise Exception("Can't match gene graph: {}".format(g))
                np = tries[g].match(v)  # Match in the new graph now
                obj = RType(g, v, e, np, w)
                if as_dict:
                    cur.append(obj._asdict)
                else:
                    cur.append(obj)
            if len(cur) > 0:
                ret.append(cur)
    return ret

def _list_eq_classes(ll, gid, unique = True):
    '''
    Function for debugging - outputs all equivalance classes under a gene graph.
    '''
    f = io.StringIO()
    seen = set()
    for l in ll:
        for rec in l:
            if rec.gid == gid:
                kw = tuple(rec.vpath)
                if (kw not in seen) or (not unique):
                    seen.add(kw)
                    print(rec, file=f)
    ret = f.getvalue()
    f.close()
    return ret

def work(graph_fn, eq_class_fn, out_prefix, dep_free = False, debug = False):
    '''
    Main working routine. Reads all stuff and writes the new graph and new
    equivalance classes as pickled objects.
    @param graph_fn: File name for graph_fragstart.txt.
    @param eq_class_fn: Filename for eq_classes.txt.
    @param out_prefix: Prefix for outputted files.
    @param dep_free: By default the pickled object can only be loaded with this
                    module present(but also will support methods from here).
                    If set to True the pickles will be of a different type.
    @param debug: Defaults to False, set to True to output some diagonistics.
    @return: Nothing. The outputs will be pickled and serialized to file.
            out_prefix + _prefix_tries.dat: Dict(graph ID -> PrefixTrie object)
                        or graph ID -> {"nodes": {list of nodes}, "edges": {list of edges}}
            out_prefix + _translated_classes.dat: List of list of Eq-Classes.
                        See document for translate_eq_classes for details.
    '''
    st = time.perf_counter()
    edges = load_graphs(graph_fn)
    hedges = load_eq_classes(eq_class_fn)
    print("{:.2f} File loaded".format(time.perf_counter() - st))
    tries = dict()
    c = 0
    for gname in edges:
        he = hedges.get(gname, [])
        tries[gname] = PrefixTrie(edges[gname], he, name = gname)
        c += 1
        if (c <= 15) and debug:
            with open(gname + ".txt", "w") as f:
                print(tries[gname]._summary(), file=f)
    print("{:.2f} Graph constructed".format(time.perf_counter() - st))
    translated = translate_eq_classes(tries, eq_class_fn, dep_free)
    print("{:.2f} Equivalance Classes translated".format(time.perf_counter() - st))
    if debug:
        with open("new_eq_samples.txt", "w") as f:
            for d in range(1000):
                print("ID =", d, file=f)
                for item in translated[d]:
                    print(item, file=f)
    with open(out_prefix + "_prefix_tries.dat", "bw") as f:
        if dep_free:
            t = dict()
            for k, v in tries.items():
                t[k] = {"nodes": v.nodes, "edges": v.edges}
            pickle.dump(t, f)
        else:
            pickle.dump(tries, f)
    with open(out_prefix + "_translated_classes.dat", "bw") as f:
        pickle.dump(translated, f)
    print("{:.2f} Data dumped".format(time.perf_counter() - st))

def work_merged(graph_fn, eq_class_fn, out_prefix, debug = False):
    '''
    Main working routine with a merged prefix tree. Parameters are identical to
    work().
    '''
    st = time.perf_counter()
    edges = load_graphs(graph_fn)
    hedges = load_eq_classes(eq_class_fn)
    print("{:.2f} File loaded".format(time.perf_counter() - st))
    if debug:
        new_edges = dict()
        new_hedges = dict()
        i = 0
        for gid in edges:
            i += 1
            if i == 20:
                break
            new_edges[gid] = edges[gid]
            new_hedges[gid] = hedges[gid]
        edges = new_edges
        hedges = new_hedges
    mt = MergedTrie(edges, hedges)
    print("{:.2f} Graph constructed".format(time.perf_counter() - st))
    translated = translate_eq_classes(mt.get_wrappers(), eq_class_fn,
                                      allow_discards=debug)
    print("{:.2f} Equivalance Classes translated".format(time.perf_counter() - st))
    if debug:
        with open("merged_trie.txt", "w") as f:
            print(mt._summary(), file=f)
        with open("new_eq_samples.txt", "w") as f:
            for d in range(min(1000, len(translated))):
                print("ID =", d, file=f)
                for item in translated[d]:
                    print(item, file=f)
    with open(out_prefix + "_prefix_tries.dat", "bw") as f:
        pickle.dump(mt, f)
    with open(out_prefix + "_translated_classes.dat", "bw") as f:
        pickle.dump(translated, f)
    print("{:.2f} Data dumped".format(time.perf_counter() - st))



def load_data(prefix):
    with open(prefix + "_prefix_tries.dat", "rb") as f:
        ret1 = pickle.load(f)
    with open(prefix + "_translated_classes.dat", "rb") as f:
        ret2 = pickle.load(f)
    return ret1, ret2

if __name__ == "__main__":
    #  work("sim2/gs_graph_fragstart.txt", "sim2/eq_classes.txt", "sim2/test_recon",
         #  debug = True)
    work_merged("sim2/gs_graph_fragstart.txt", "sim2/eq_classes.txt", "sim2/test_merged",
                debug = True)
