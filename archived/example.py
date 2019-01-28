v0 = Node_t(0, "chr1", -1, -1, True, False, [], [0,1], 1)
v1 = Node_t(1, "chr1", 0, 100, True, True, [0], [2], 1)
v2 = Node_t(2, "chr1", 100, 150, True, True, [1], [3], 1)
v3 = Node_t(3, "chr1", 150, 200, True, True, [2], [4], 1)
v4 = Node_t(4, "chr1", 200, 400, True, True, [3,4], [5,6], 1)
v5 = Node_t(5, "chr1", 400, 600, True, True, [5], [7], 1)
v6 = Node_t(6, "chr1", 600, 800, True, True, [6], [8], 1)
v7 = Node_t(7, "chr1", 2147483647, 2147483647, True, False, [7,8], [], 1)

e0 = Edge_t(0, 0, 1, [])
e1 = Edge_t(1, 0, 2, [])
e2 = Edge_t(2, 1, 3, [])
e3 = Edge_t(3, 2, 4, [])
e4 = Edge_t(4, 3, 4, [])
e5 = Edge_t(5, 4, 5, [])
e6 = Edge_t(6, 4, 6, [])
e7 = Edge_t(7, 5, 7, [])
e8 = Edge_t(8, 6, 7, [])

vNodes = [v0, v1, v2, v3, v4, v5, v6, v7]
vEdges = [e0, e1, e2, e3, e4, e5, e6, e7, e8]
g = GeneGraph_t("test", vNodes, vEdges)

opt = CVXOPTObject_hyper(g, [])
opt.PathGroups = [([1]), ([3]), ([4]), ([5])]
opt.Weights = np.array([1, 2, 20, 50])
H, bound, A, lb = opt.constraints()
cvxopt.solvers.options['show_progress'] = True
sol = cvxopt.solvers.cp(opt.F, G=(-A), h=(-lb), A=H, b=bound)