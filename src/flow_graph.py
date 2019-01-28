import networkx as nx
from math import log
import random
import itertools
import collections

def module_desc():
    '''
    This module handles stuff related to sampling and quantification on flow
    graphs, and also some handy functions to translate between gene graphs and
    flow graphs. I don't know. Those probably should go to the Trie source file
    but that thing is already 600 lines. We will see.
    '''

eps = 1e-7  # placeholder flow value

class FlowGraph:
    '''
    Main class. Initialize with (edges, flows, S, T). Don't care about what does
    the edge of graphs mean, just get the flows done.
    WARNING: We assume the each vertex in G reaches T and is reachable from S,
    even if its flow is zero. Additionally we assume the edges are topologically
    sorted.
    WARNING: If you want to test another flow, you need to build another instance
    of FlowGraph. Changing self.flow will NOT work.
    '''
    def __init__(self, edges, flow, S, T):
        self.edges = edges
        self.flow = flow
        self.S, self.T = S, T
        assert len(self.edges) == len(self.flow)
        self.m = len(self.edges)
        m = len(self.edges)
        self.fw_edges = collections.defaultdict(list)
        self.bw_edges = collections.defaultdict(list)
        for i in range(m):
            u, v = self.edges[i]
            f = self.flow[i]
            self.fw_edges[u].append((i, v, max(f, eps), len(self.fw_edges[u])))
            self.bw_edges[v].append((i, u, max(f, eps), len(self.bw_edges[v])))
        for v in self.fw_edges:
            if v != self.S:
                assert v != self.T
                # At flow balance check don't check those with pity flows
                out_flow = sum(d[2] for d in self.fw_edges[v] if abs(d[2]-eps) > eps / 100)
                in_flow = sum(d[2] for d in self.bw_edges[v] if abs(d[2]-eps) > eps / 100)
                assert abs(out_flow - in_flow) < 10 * eps
        # build the unweighted flow graph
        self._uw_flows = collections.defaultdict(float)
        self._uw_flows[self.S] = 1
        self.uw_fw_edges = collections.defaultdict(list)
        self.uw_bw_edges = collections.defaultdict(list)
        for i in range(m):  # Can do this because edges are topo-sorted
            u, v = self.edges[i]
            f = self._uw_flows[u] / len(self.fw_edges[u])  # Equipartition
            self._uw_flows[v] += f
            self.uw_fw_edges[u].append((i, v, f, len(self.uw_fw_edges[u])))
            self.uw_bw_edges[v].append((i, u, f, len(self.uw_bw_edges[v])))
        # helper variable
        self.full_flow = sum(self.flow[i] for i in range(self.m)
                              if self.edges[i][0] == self.S)

    def _directed_sample(self, v, d, weighted, rng):
        '''
        Sample a path starting from v and ending at T.
        @param v: Starting vertex.
        @param d: Direction of sampling. True for forward, False for backward.
        @param weighted: If path should be sampled by flow value.
        @param rng: random.Random instance to make choices.
        @return list of edges, log-likelihood.
        '''
        ll = []
        lw = 0
        cur = v
        stop = self.T if d else self.S
        while cur != stop:
            if d:
                fd = self.fw_edges[cur] if weighted else self.uw_fw_edges[cur]
            else:
                fd = self.bw_edges[cur] if weighted else self.uw_bw_edges[cur]
            #  print(cur, fd)
            total_weight = sum(x[2] for x in fd)
            sp = rng.random() * total_weight
            for e in fd:
                sp -= e[2]
                if sp <= 0:
                    break
            assert (sp <= 0)
            #  e = rng.choices(fd, list(x[2] for x in fd))[0]

            ll.append(e[0])  # edge ID
            ws = self.fw_edges[cur] if d else self.bw_edges[cur]
            cur = e[1]
            idx = e[3]
            #  print(ws, idx)
            lw += log(ws[idx][2] / sum(x[2] for x in ws))
        return ll, lw


    def sample_path(self, N, weighted, seed = None, vlist = False):
        '''
        Sample S-T path from the graph under Markov model.
        @param N: number of samples.
        @param weighted: If True, path will be sampled by flow values and if
                        False, path will be sampled by randomly picking out-degrees.
        @param seed: Seed to initialize RNG. Leave it to None to initialize with
                        current time.
        @param vlist: Set to True to return a list of vertexes instead.
        @return: elist, log_weights. First one is a list of list indicating the
                    edges in sampled path, second one is a list of *LOG-*likelihood
                    for that path under weighted model(even if weighted = False).
        '''
        rng = random.Random()
        rng.seed(seed)
        ret_list = []
        ret_lw = []
        for i in range(N):
            ll, lw = self._directed_sample(self.S, True, weighted, rng)
            ret_list.append(ll)
            ret_lw.append(lw)
        if vlist:
            ret_list = self.translate_path(ret_list)
        return ret_list, ret_lw

    def sample_centered_path(self, N, x, weighted, seed = None, vlist = False):
        '''
        Sample S-T path from the graph under Markov model, conditioned on an edge
        @param x: index of the edge(consistent with self.edges).
        @param N, weighted, seed, vlist: Same as sample_path.
        @return: elist, log_weights, same as sample_path.
        '''
        rng = random.Random()
        rng.seed(seed)
        ret_list = []
        ret_lw = []
        u, v = self.edges[x]
        # likelihood for the edge u->v
        lw2 = log(max(self.flow[x], eps) / sum(d[2] for d in self.fw_edges[u]))
        for i in range(N):
            ll0, lw0 = self._directed_sample(u, False, weighted, rng)
            ll1, lw1 = self._directed_sample(v, True, weighted, rng)
            ret_list.append(list(reversed(ll0)) + [x] + ll1)
            ret_lw.append(lw0 + lw1 + lw2)
        if vlist:
            ret_list = self.translate_path(ret_list)
        return ret_list, ret_lw

    def diff_flow(self, query):
        '''
        Lower-Bounding OR-Quantification.
        @param query: The list of edges for quantification query.
        '''
        l = set(query)
        g = nx.DiGraph()
        for i in range(self.m):
            if i not in l:
                g.add_edge(self.edges[i][0], self.edges[i][1], capacity = self.flow[i])
                #  print(self.edges[i][0], self.edges[i][1], self.flow[i])
        #  print(g)
        if (self.S in g.nodes) and (self.T in g.nodes):
            f = nx.algorithms.maximum_flow_value(g, self.S, self.T)
        else:
            #  print("WARNING - ", self.S, self.T, "not in graph - stopping")
            f = 0
        return self.full_flow - f

    def split_flow(self, query):
        '''
        Upper-Bounding OR-Quantification.
        @param query: The list of edges for quantification query.
        '''
        l = set(query)
        diverted_flows = collections.defaultdict(float)
        g = nx.DiGraph()
        for i in range(self.m):
            if i in l:
                diverted_flows[self.edges[i][0]] += self.flow[i]
                #  g.add_edge(self.edges[i][0], -1, capacity = self.flow[i])
            else:
                g.add_edge(self.edges[i][0], self.edges[i][1], capacity = self.flow[i])
        for source, fv in diverted_flows.items():
            g.add_edge(source, -1, capacity = fv)
        if (self.S in g.nodes) and (-1 in g.nodes):
            f = nx.algorithms.maximum_flow_value(g, self.S, -1)
        else:
            f = 0
        return f

    def _get_block(self, qlist):
        '''
        Get the blocking subgraph for AND-Quant.
        @param qlist: List of edges for AND-Quantification.
        @return: list of edges in blocking subgraph.
        '''
        s = len(qlist)
        fw_label = collections.defaultdict(lambda : -1)
        bw_label = collections.defaultdict(lambda : s)
        for i in range(s):
            for e in qlist[i]:
                u, v = self.edges[e]
                fw_label[v] = i
                bw_label[u] = i

        # forward loop: determine fw_label = max i: V_i -> x
        for u, v in self.edges:
            fw_label[v] = max(fw_label[v], fw_label[u])

        # backward loop: determine bw_label = min i: x -> U_i
        for u, v in reversed(self.edges):
            bw_label[u] = min(bw_label[u], bw_label[v])

        ret = set()
        # sanity check
        for i in range(s):
            for e in qlist[i]:
                ret.add(e)
                u, v = self.edges[e]
                assert fw_label[u] == i - 1
                assert bw_label[v] == i + 1

        #  print(fw_label)
        #  print(bw_label)
        # check the transition edges: V_i -> u -> v -> U_i+1
        for i in range(self.m):
            u, v = self.edges[i]
            assert fw_label[u] < bw_label[v]
            if bw_label[v] - fw_label[u] == 1:
                ret.add(i)

        return list(ret)

    def block_flow(self, qlist):
        '''
        Upper-Bounding AND-Quantification.
        @param qlist: List of list of edges for each "phase".
        '''
        blocks = set(self._get_block(qlist))
        r = []
        for i in range(self.m):
            if i not in blocks:
                r.append(i)
        return self.full_flow - self.diff_flow(r)

    def split_block_flow(self, qlist):
        '''
        Lower-Bounding AND-Quantification.
        @param qlist: List of list of edges for each "phase".
        '''
        blocks = set(self._get_block(qlist))
        r = []
        for i in range(self.m):
            if i not in blocks:
                r.append(i)
        return self.full_flow - self.split_flow(r)

    def translate_path(self, paths):
        ret = []
        for p in paths:
            ret.append([self.S] + list(self.edges[x][1] for x in p))
            #  c = []
            #  for e in p:
                #  c.append(self.edges[e][0])
                #  c.append(self.edges[e][1])
            #  ret.append(c)
        return ret


if __name__ == "__main__":
    edges = [(1, 2), (1, 3), (1, 4), (2, 3), (2, 5), (3, 4), (4, 5), (4, 6), (5, 6), (6, 7)]
    flows = [1, 2, 1, 1, 0, 3, 2, 2, 2, 4]
    fg = FlowGraph(edges, flows, 1, 7)
    #  print(fg.sample_path(10, weighted = True, vlist = True))
    #  print(fg.sample_path(10, weighted = False, vlist = True))
    #  print(fg.sample_centered_path(10, 4, weighted = False, vlist = True))
    #  print(fg.sample_centered_path(10, 4, weighted = True, vlist = True))
    #  print(fg.sample_centered_path(10, 0, weighted = False, vlist = True))
    #  print(fg.sample_centered_path(10, 0, weighted = True, vlist = True))
    #  print(fg.sample_centered_path(10, 6, weighted = False, vlist = True))
    #  print(fg.sample_centered_path(10, 6, weighted = True, vlist = True))
    print(fg.diff_flow([1, 2, 7]))
    print(fg.split_flow([1, 2, 7]))
    print(fg.block_flow([[1], [7], [9]]))
    print(fg.split_block_flow([[1], [7], [9]]))
