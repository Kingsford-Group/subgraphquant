#!/bin/python

import sys
import numpy as np

floatpoint_error = 1e-15

def ReadSalmonQuant(filename):
	NumReads = {}
	Expression = {}
	fp = open(filename, 'r')
	linecount = 0
	for line in fp:
		linecount += 1
		if linecount == 1:
			continue
		strs = line.strip().split("\t")
		NumReads[strs[0]] = float(strs[4])
		Expression[strs[0]] = float(strs[4]) / float(strs[2])
	fp.close()
	return NumReads, Expression


def ReadAbundanceGap(filename):
	AbundanceGap = {}
	fp = open(filename, 'r')
	for line in fp:
		if line[0] == '#':
			continue
		strs = line.strip().split("\t")
		lb = float(strs[1])
		ub = float(strs[2])
		if lb < 0:
			lb = 0
		if ub < 0:
			ub = 0
		if lb > ub:
			ub = lb
		AbundanceGap[strs[0]] = (lb, ub)
	fp.close()
	return AbundanceGap


def NormalizeAbundanceGap(AbundanceGap, NumReads):
	'''
	normalization of salmon flow:
		the flow of salmon is calculated by transcript numreads / effective length (sum of flow * efflen over all transcripts = numreads)
	normalization of graph salmon flow:
		sum of (flow * efflen) = numreads
	different samples have different numreads
	post-normalize lb and ub is equivalent to pre-normalize the flow and bound using the normalized flows
	normalize by the numreads so that different samples have the same number of numreads
	'''
	sumreads = np.sum(list(NumReads.values()))
	AbundanceGap = {t:(v[0] / sumreads * 1e6, v[1] / sumreads * 1e6) for t,v in AbundanceGap.items()}
	return AbundanceGap


def NormalizeExpression(Expression, NumReads):
	sumreads = np.sum(list(NumReads.values()))
	Expression = {t:v / sumreads * 1e6 for t,v in Expression.items()}
	return Expression


def AdjustBounds(AbundanceGap, Expression):
	for t,v in Expression.items():
		lb, ub = AbundanceGap[t]
		if lb < 0 or lb > v:
			lb = min(0, v)
		if ub < lb or ub < v:
			ub = max(lb, v)
		AbundanceGap[t] = (lb, ub)
	return AbundanceGap


def ReadCenterMetadata(filename):
	'''
	this is for GEUVADIS dataset
	'''
	CenterSampleMap = {}
	fp = open(filename, 'r')
	for line in fp:
		if line[0] == '#':
			continue
		strs = line.strip().split("\t")
		if strs[1] in CenterSampleMap:
			CenterSampleMap[strs[1]].append(strs[-1])
		else:
			CenterSampleMap[strs[1]] = [strs[-1]]
	fp.close()
	return CenterSampleMap


def ReadTreatmentMetadata(filename):
	'''
	this is for MCF10 dataset
	'''
	TreatmentSampleMap = {}
	fp = open(filename, 'r')
	for line in fp:
		strs = line.strip().split("\t")
		if strs[7] in TreatmentSampleMap:
			TreatmentSampleMap[strs[7]].append(strs[4])
		else:
			TreatmentSampleMap[strs[7]] = [strs[4]]
	fp.close()
	return TreatmentSampleMap


def ProcessUncertaintyMatrix(GroupSampleMap, folder, id_prefix = "", bound_suffix = "prefixgraph/gs_maxflow_bound.txt"):
	transnames = None
	Expression_pergroup = {}
	Lb_pergroup = {}
	Ub_pergroup = {}
	for c,samples in GroupSampleMap.items():
		mat_exp = None
		mat_lb = None
		mat_ub = None
		for s in samples:
			salmonfile = folder + "/" + id_prefix + s + "/quant.sf"
			boundfile = folder + "/" + id_prefix + s + "/" + bound_suffix
			NumReads, Expression = ReadSalmonQuant(salmonfile)
			Expression = NormalizeExpression(Expression, NumReads)
			Bounds = ReadAbundanceGap(boundfile)
			Bounds = NormalizeAbundanceGap(Bounds, NumReads)
			# Bounds = AdjustBounds(Bounds, Expression)
			print("finish reading {}:{}".format(c, s))
			# order the transcript names by dictionary order
			tmptransnames = list(Expression.keys())
			tmptransnames.sort()
			if transnames is None:
				transnames = tmptransnames
			else:
				assert( len(transnames) == len(tmptransnames) and np.all([transnames[i] == tmptransnames[i] for i in range(len(transnames))]) )
			# process matrix
			vec_exp = np.array([Expression[t] for t in transnames]).reshape( (len(transnames), 1) )
			vec_lb = np.array([Bounds[t][0] for t in transnames]).reshape( (len(transnames), 1) )
			vec_ub = np.array([Bounds[t][1] for t in transnames]).reshape( (len(transnames), 1) )
			if mat_exp is None:
				mat_exp = vec_exp
				mat_lb = vec_lb
				mat_ub = vec_ub
			else:
				mat_exp = np.hstack( (mat_exp, vec_exp) )
				mat_lb = np.hstack( (mat_lb, vec_lb) )
				mat_ub = np.hstack( (mat_ub, vec_ub) )
		Expression_pergroup[c] = mat_exp
		Lb_pergroup[c] = mat_lb
		Ub_pergroup[c] = mat_ub
	return transnames, Expression_pergroup, Lb_pergroup, Ub_pergroup


def MinLowerBounds_varylambda(salmon_exps, lbs):
	'''
	Assuming Salmon must be (1-lambda) proportion, the lower bounds becomes (1-lambda) * salmon + lambda * lb
	This function calculates the minimum lower bounds from multiple samples with varying lambda proportions.
	The min lower bounds of the samples is a piece-wise linear functions with respect to lambda.
	input:
		- salmon_exps: Salmon expression of transcripts from multiple samples
		- lbs: the G2 lower bounds of the transcript from multiple samples
	output:
		- region
		- interceptions: the interception of each linear functions in the piece-wise linear curves
		- slopes: the slope of each linear functions in the piece-wise linear curves
	'''
	assert( len(salmon_exps) == len(lbs) )
	assert( np.all([lbs[i] <= salmon_exps[i] for i in range(len(salmon_exps))]) )
	# variables to record which sample has the min lb, and the corresponding interception and slope of the lb line
	indexes = []
	region_start = []
	region_end = []
	interceptions = []
	slopes = []
	# from lambda = 0, find the minimum salmon expression
	s = 0
	region_start.append(s)
	i = np.argmin(salmon_exps)
	indexes.append( i )
	interceptions.append( salmon_exps[i] )
	slopes.append( lbs[i] - salmon_exps[i] )
	# record the possible index of lines that can cross with the current min
	set_possible_indexes = set(list(range(len(salmon_exps)))) - set([i])
	while s < 1:
		# find the region_end: either 1 or the first crossing with other lines
		crossing = []
		for j in set_possible_indexes:
			# parallel case
			if lbs[j] - salmon_exps[j] == lbs[i] - salmon_exps[i]:
				continue
			# find x axis of the crossing
			crossing.append( (j, 1.0 * (salmon_exps[j] - salmon_exps[i]) / (lbs[i] - salmon_exps[i] - lbs[j] + salmon_exps[j]) ) )
		crossing.sort(key = lambda x:x[1])
		# remove the crossing that before s: these lines are not possible to be below the current line
		for (k, x) in crossing:
			if x < s:
				set_possible_indexes.remove(k)
		crossing = [ (k,x) for (k,x) in crossing if x > s ]
		# find the crossing
		if len(crossing) == 0 or crossing[0][1] >= 1:
			region_end.append( 1 )
			s = 1
		else:
			region_end.append(crossing[0][1])
			s = region_end[-1]
			i = crossing[0][0]
			region_start.append( s )
			indexes.append( i )
			interceptions.append( salmon_exps[i] )
			slopes.append( lbs[i] - salmon_exps[i] )
			set_possible_indexes.remove( i )
	assert(len(indexes) == len(region_start))
	assert(len(region_start) == len(region_end))
	assert(len(region_end) == len(interceptions))
	assert(len(interceptions) == len(slopes))
	# assertion about the monotonic decreasing of min lb as lambda increases
	for i in range(1, len(region_start)):
		minlb_start = interceptions[i-1] + slopes[i-1] * region_start[i-1]
		minlb_end = interceptions[i] + slopes[i] * region_start[i]
		assert(minlb_start >= minlb_end - floatpoint_error)
	regions = [(region_start[i], region_end[i]) for i in range(len(region_start))]
	return regions, interceptions, slopes


def MaxUpperBounds_varylambda(salmon_exps, ubs):
	'''
	Similar to MinLowerBounds_varylambda, but calculates the maximum upper bounds
	'''
	assert( len(salmon_exps) == len(ubs) )
	assert( np.all([ubs[i] >= salmon_exps[i] for i in range(len(salmon_exps))]) )
	# variables to record which sample has the max ub, and the corresponding interception and slope of the ub line
	indexes = []
	region_start = []
	region_end = []
	interceptions = []
	slopes = []
	# from lambda = 0, find the maximum salmon expression
	s = 0
	region_start.append(s)
	i = np.argmax(salmon_exps)
	indexes.append( i )
	interceptions.append( salmon_exps[i] )
	slopes.append( ubs[i] - salmon_exps[i] )
	# record the possible index of lines that can cross with the current min
	set_possible_indexes = set(list(range(len(salmon_exps)))) - set([i])
	while s < 1:
		# find the region_end: either 1 or the first crossing with other lines
		crossing = []
		for j in set_possible_indexes:
			# parallel case
			if ubs[j] - salmon_exps[j] == ubs[i] - salmon_exps[i]:
				continue
			# find x axis of the crossing
			crossing.append( (j, 1.0 * (salmon_exps[j] - salmon_exps[i]) / (ubs[i] - salmon_exps[i] - ubs[j] + salmon_exps[j]) ) )
		crossing.sort(key = lambda x:x[1])
		# remove the crossing that before s: these lines are not possible to be above the current line
		for (k, x) in crossing:
			if x < s:
				set_possible_indexes.remove(k)
		crossing = [ (k,x) for (k,x) in crossing if x > s ]
		# find the crossing
		if len(crossing) == 0 or crossing[0][1] > 1:
			region_end.append( 1 )
			s = 1
		else:
			region_end.append(crossing[0][1])
			s = region_end[-1]
			i = crossing[0][0]
			region_start.append( s )
			indexes.append( i )
			interceptions.append( salmon_exps[i] )
			slopes.append( ubs[i] - salmon_exps[i] )
			set_possible_indexes.remove( i )
	assert(len(indexes) == len(region_start))
	assert(len(region_start) == len(region_end))
	assert(len(region_end) == len(interceptions))
	assert(len(interceptions) == len(slopes))
	# assertion about the monotonic increasing of max ub
	for i in range(1, len(region_start)):
		maxub_start = interceptions[i-1] + slopes[i-1] * region_start[i-1]
		maxub_end = interceptions[i] + slopes[i] * region_start[i]
		assert(maxub_start <= maxub_end + floatpoint_error)
	regions = [(region_start[i], region_end[i]) for i in range(len(region_start))]
	return regions, interceptions, slopes


def MeanLowerBounds_varylambda(salmon_exps, lbs):
	assert(len(salmon_exps) == len(lbs))
	interception = np.mean(salmon_exps)
	slope = np.mean([lbs[i] - salmon_exps[i] for i in range(len(lbs))])
	return interception, slope


def MeanUpperBounds_varylambda(salmon_exps, ubs):
	assert(len(salmon_exps) == len(ubs))
	interception = np.mean(salmon_exps)
	slope = np.mean([ubs[i] - salmon_exps[i] for i in range(len(ubs))])
	return interception, slope


def CalculateIValue_minmax(region1_min, interceptions1_min, slopes1_min, region2_min, interceptions2_min, slopes2_min, \
	region1_max, interceptions1_max, slopes1_max, region2_max, interceptions2_max, slopes2_max):
	'''
	This functions calculate the minimum lambda value that make the DE transcripts uncertainty range overlap between the case and control.
	input:
		- interceptions1_min, slopes1_min: the function to record the min lower bounds of the transcript in the case group, with varying lambda
		- interceptions1_max, slopes1_max: the function to record the max upper bounds of the transcript in the case group
		- interceptions2_min, slopes2_min: the function to record the min lower bounds of the transcript in the control group
		- interceptions2_max, slopes2_max: the function to record the max upper bounds of the transcript in the control group
	output:
		- IValue: the minimum lambda to make the case bounds and control bounds overlap
	'''
	union_breaks = [x[0] for x in region1_min] + [x[0] for x in region1_max] + [x[0] for x in region2_min] + [x[0] for x in region2_max]
	union_breaks = list(set(union_breaks))
	union_breaks.sort()
	union_region = [ (union_breaks[j], union_breaks[j+1]) for j in range(len(union_breaks)-1) ] + [ (union_breaks[-1], 1) ]
	# duplicate the interceptions and slopes for the union_breaks
	# min lb of case
	new_int1_min = []
	new_slope1_min = []
	idx = 0
	for j in range(len(union_region)):
		if union_region[j][0] >= region1_min[idx][1]:
			idx += 1
		assert(idx < len(region1_min) and union_region[j][0] >= region1_min[idx][0] and union_region[j][1] <= region1_min[idx][1])
		new_int1_min.append( interceptions1_min[idx] )
		new_slope1_min.append( slopes1_min[idx])
	# min lb of control
	new_int2_min = []
	new_slope2_min = []
	idx = 0
	for j in range(len(union_region)):
		if union_region[j][0] >= region2_min[idx][1]:
			idx += 1
		assert(idx < len(region2_min) and union_region[j][0] >= region2_min[idx][0] and union_region[j][1] <= region2_min[idx][1])
		new_int2_min.append( interceptions2_min[idx] )
		new_slope2_min.append( slopes2_min[idx])
	# max ub of case
	new_int1_max = []
	new_slope1_max = []
	idx = 0
	for j in range(len(union_region)):
		if union_region[j][0] >= region1_max[idx][1]:
			idx += 1
		assert(idx < len(region1_max) and union_region[j][0] >= region1_max[idx][0] and union_region[j][1] <= region1_max[idx][1])
		new_int1_max.append( interceptions1_max[idx] )
		new_slope1_max.append( slopes1_max[idx])
	# max ub of control
	new_int2_max = []
	new_slope2_max = []
	idx = 0
	for j in range(len(union_region)):
		if union_region[j][0] >= region2_max[idx][1]:
			idx += 1
		assert(idx < len(region2_max) and union_region[j][0] >= region2_max[idx][0] and union_region[j][1] <= region2_max[idx][1])
		new_int2_max.append( interceptions2_max[idx] )
		new_slope2_max.append( slopes2_max[idx])
	# find the first region in the union_region that case and control uncertainty range overlaps
	IV = 1
	for i in range(len(union_region)):
		lb1_end = new_int1_min[i] + new_slope1_min[i] * union_region[i][1]
		ub1_end = new_int1_max[i] + new_slope1_max[i] * union_region[i][1]
		lb2_end = new_int2_min[i] + new_slope2_min[i] * union_region[i][1]
		ub2_end = new_int2_max[i] + new_slope2_max[i] * union_region[i][1]
		# if not overlap
		if max(lb1_end, lb2_end) > min(ub1_end, ub2_end):
			continue
		# if overlap, find the first lambda of overlapping
		else:
			lb1_start = new_int1_min[i] + new_slope1_min[i] * union_region[i][0]
			ub1_start = new_int1_max[i] + new_slope1_max[i] * union_region[i][0]
			lb2_start = new_int2_min[i] + new_slope2_min[i] * union_region[i][0]
			ub2_start = new_int2_max[i] + new_slope2_max[i] * union_region[i][0]
			# if the case and control overlap at the start of the region
			if max(lb1_start, lb2_start) <= min(ub1_start, ub2_start):
				IV = union_region[i][0]
				break
			# if case is below control: ub1_start < lb2_start
			elif ub1_start < lb2_start:
				x = 1.0 * (new_int1_max[i] - new_int2_min[i]) / (new_slope2_min[i] - new_slope1_max[i])
				assert(x >= union_region[i][0] and x <= union_region[i][1])
				IV = x
				break
			# if control is below case: ub2_start < lb1_start
			elif ub2_start < lb1_start:
				x = 1.0 * (new_int1_min[i] - new_int2_max[i]) / (new_slope2_max[i] - new_slope1_min[i])
				assert(x >= union_region[i][0] and x <= union_region[i][1])
				IV = x
				break
			else:
				print("Error: the two uncertainty ranges overlap before this region")
				sys.exit()
	return IV


# def CalculateIValue_mean(interception1_min, slope1_min, interception2_min, slope2_min, interception1_max, slope1_max, interception2_max, slope2_max):
# 	# lb and ub at the region start
# 	lb1_start = interception1_min
# 	ub1_start = interception1_max
# 	lb2_start = interception2_min
# 	ub2_start = interception2_max
# 	assert( lb1_start <= ub1_start and lb2_start <= ub2_start )
# 	# lb and ub at the region end
# 	lb1_end = interception1_min + slope1_min
# 	ub1_end = interception1_max + slope1_max
# 	lb2_end = interception2_min + slope2_min
# 	ub2_end = interception2_max + slope2_max
# 	assert( lb1_end <= ub1_end and lb2_end <= ub2_end )
# 	# calculate IV
# 	IV = 1
# 	# overlap in the beginning
# 	if max(lb1_start, lb2_start) <= min(ub1_start, ub2_start):
# 		IV = 0
# 	elif ub1_start < lb2_start:
# 		# not overlap at all through all lambda values
# 		if ub1_end < lb2_end:
# 			IV = 2
# 		# overlap in the middle
# 		else:
# 			x = (interception2_min - interception1_max) / (slope1_max - slope2_min)
# 			assert(x >= 0 and x <= 1)
# 			IV = x
# 	elif ub2_start < lb1_start:
# 		# not overlap at all through all lambda values
# 		if ub2_end < lb1_end:
# 			IV = 2
# 		# overlap in the middle
# 		else:
# 			x = (interception2_max - interception1_min) / (slope1_min - slope2_max)
# 			assert(x >= 0 and x <= 1)
# 			IV = x
# 	return IV


def CalculateIValue_mean(interception1_min, slope1_min, interception2_min, slope2_min, interception1_max, slope1_max, interception2_max, slope2_max, quantile
= 0):
	# lb and ub at the region start
	lb1_start = interception1_min
	ub1_start = interception1_max
	lb2_start = interception2_min
	ub2_start = interception2_max
	# lb and ub at the region end
	lb1_end = interception1_min + slope1_min
	ub1_end = interception1_max + slope1_max
	lb2_end = interception2_min + slope2_min
	ub2_end = interception2_max + slope2_max
	# adjust for quantile for start
	q1_start = (ub1_start - lb1_start) * quantile
	q2_start = (ub2_start - lb2_start) * quantile
	lb1_start += q1_start
	ub1_start -= q1_start
	assert(lb1_start <= ub1_start)
	lb2_start += q2_start
	ub2_start -= q2_start
	assert(lb2_start <= ub2_start)
	# adjust for quantile for end
	q1_end = (ub1_end - lb1_end) * quantile
	q2_end = (ub2_end - lb2_end) * quantile
	lb1_end += q1_end
	ub1_end -= q1_end
	assert(lb1_end <= ub1_end)
	lb2_end += q2_end
	ub2_end -= q2_end
	assert(lb2_end <= ub2_end)
	# adjust for quantile for slope
	q_slope1_min = lb1_end - lb1_start
	q_slope1_max = ub1_end - ub1_start
	q_slope2_min = lb2_end - lb2_start
	q_slope2_max = ub2_end - ub2_start
	# calculate IV
	IV = 1
	# overlap in the beginning
	if max(lb1_start, lb2_start) <= min(ub1_start, ub2_start):
		IV = 0
	elif ub1_start < lb2_start:
		# not overlap at all through all lambda values
		if ub1_end < lb2_end:
			IV = 2
		# overlap in the middle
		else:
			x = (lb2_start - ub1_start) / (q_slope1_max - q_slope2_min)
			assert(x >= 0 and x <= 1)
			IV = x
	elif ub2_start < lb1_start:
		# not overlap at all through all lambda values
		if ub2_end < lb1_end:
			IV = 2
		# overlap in the middle
		else:
			x = (ub2_start - lb1_start) / (q_slope1_min - q_slope2_max)
			assert(x >= 0 and x <= 1)
			IV = x
	return IV

def WriteIV_minmax(outputfile, transnames, Expression_pergroup, Lb_pergroup, Ub_pergroup):
	assert(len(Expression_pergroup) == 2 and len(Lb_pergroup) == 2 and len(Ub_pergroup) == 2)
	groups = list(Expression_pergroup.keys())
	groups.sort()
	# process and write IV for each transcript
	fp = open(outputfile, 'w')
	fp.write("# Name\t{}Mean\t{}Mean\tIV\n".format(groups[0], groups[1]))
	for i in range(len(transnames)):
		t = transnames[i]
		c1_exp = Expression_pergroup[groups[0]][i,:]
		c1_lb = Lb_pergroup[groups[0]][i,:]
		c1_ub = Ub_pergroup[groups[0]][i,:]
		c2_exp = Expression_pergroup[groups[1]][i,:]
		c2_lb = Lb_pergroup[groups[1]][i,:]
		c2_ub = Ub_pergroup[groups[1]][i,:]
		# min lb of c1
		region1_min, interceptions1_min, slopes1_min = MinLowerBounds_varylambda(c1_exp, c1_lb)
		# max ub for c1
		region1_max, interceptions1_max, slopes1_max = MaxUpperBounds_varylambda(c1_exp, c1_ub)
		# min lb for c2
		region2_min, interceptions2_min, slopes2_min = MinLowerBounds_varylambda(c2_exp, c2_lb)
		# max ub for c2
		region2_max, interceptions2_max, slopes2_max = MaxUpperBounds_varylambda(c2_exp, c2_ub)
		# calculate IV
		IV = CalculateIValue_minmax(region1_min, interceptions1_min, slopes1_min, region2_min, interceptions2_min, slopes2_min, region1_max, interceptions1_max, slopes1_max, region2_max, interceptions2_max, slopes2_max)
		fp.write("{}\t{}\t{}\t{}\n".format(t, np.mean(c1_exp), np.mean(c2_exp), IV))
	fp.close()


def WriteIV_mean(outputfile, transnames, Expression_pergroup, Lb_pergroup, Ub_pergroup):
	assert(len(Expression_pergroup) == 2 and len(Lb_pergroup) == 2 and len(Ub_pergroup) == 2)
	groups = list(Expression_pergroup.keys())
	groups.sort()
	# process and write IV for each transcript
	fp = open(outputfile, 'w')
	fp.write("# Name\t{}Mean\t{}Mean\tIV\tIV25\n".format(groups[0], groups[1]))
	for i in range(len(transnames)):
		t = transnames[i]
		c1_exp = Expression_pergroup[groups[0]][i,:]
		c1_lb = Lb_pergroup[groups[0]][i,:]
		c1_ub = Ub_pergroup[groups[0]][i,:]
		c2_exp = Expression_pergroup[groups[1]][i,:]
		c2_lb = Lb_pergroup[groups[1]][i,:]
		c2_ub = Ub_pergroup[groups[1]][i,:]
		# min lb of c1
		interception1_min, slope1_min = MeanLowerBounds_varylambda(c1_exp, c1_lb)
		# max ub for c1
		interception1_max, slope1_max = MeanUpperBounds_varylambda(c1_exp, c1_ub)
		# min lb for c2
		interception2_min, slope2_min = MeanLowerBounds_varylambda(c2_exp, c2_lb)
		# max ub for c2
		interception2_max, slope2_max = MeanUpperBounds_varylambda(c2_exp, c2_ub)
		# calculate IV
		IV = CalculateIValue_mean(interception1_min, slope1_min, interception2_min, slope2_min, interception1_max, slope1_max, interception2_max, slope2_max)
		IV25 = CalculateIValue_mean(interception1_min, slope1_min, interception2_min, slope2_min, interception1_max, slope1_max, interception2_max, slope2_max, quantile = 0.25)
		fp.write("{}\t{}\t{}\t{}\t{}\n".format(t, np.mean(c1_exp), np.mean(c2_exp), IV, IV25))
	fp.close()


def WriteExampleCurve(outputfile, indexes, GroupSampleMap, transnames, Expression_pergroup, Lb_pergroup, Ub_pergroup):
	fp = open(outputfile, 'w')
	fp.write("# Name\tGroup\tSample\treference_proportion\tlb\tub\n")
	for idx in indexes:
		for c,samples in GroupSampleMap.items():
			tname = transnames[idx]
			exp = Expression_pergroup[c][idx,:]
			lb = Lb_pergroup[c][idx,:]
			ub = Ub_pergroup[c][idx,:]
			assert( len(exp) == len(samples) )
			# writing the begining and ending of line of individual samples
			for i in range(len(samples)):
				s = samples[i]
				fp.write("{}\t{}\t{}\t0\t{}\t{}\n".format(tname, c, s, exp[i], exp[i]))
				fp.write("{}\t{}\t{}\t1\t{}\t{}\n".format(tname, c, s, lb[i], ub[i]))
			# writing the line of the mean of the group
			fp.write("{}\t{}\t{}\t0\t{}\t{}\n".format(tname, c, "Mean", np.mean(exp), np.mean(exp)))
			fp.write("{}\t{}\t{}\t1\t{}\t{}\n".format(tname, c, "Mean", np.mean(lb), np.mean(ub)))
	fp.close()


if __name__=="__main__":
	if len(sys.argv) == 1:
		print("python ComputeIValue.py <output folder as in metarun.sh>")
	else:
		script_folder = sys.argv[0]
		folder = sys.argv[1]

		# get the meta file directory from script_folder
		script_folder_split = script_folder.split("/")
		if len(script_folder_split) > 2:
			script_folder = "/".join(script_folder_split[:-2])
		elif len(script_folder_split) == 2:
			script_folder = "./"
		elif len(script_folder_split) == 1:
			script_folder = "../"
		metafile = script_folder + "data/Metadata_mcf10.txt"
		print(metafile)
		print(folder + "/MCF10/IValue_mean_ext.txt")
		print(folder + "/MCF10/IValue_curve_example.txt")

		# read range of optima output
		TreatmentSampleMap = ReadTreatmentMetadata(metafile)
		transnames, Expression_pertreatment, Lb_pertreatment, Ub_pertreatment = ProcessUncertaintyMatrix(TreatmentSampleMap, folder + "/MCF10")
		WriteIV_mean(folder + "/MCF10/IValue_mean_ext.txt", transnames, Expression_pertreatment, Lb_pertreatment, Ub_pertreatment)

		# draw example expression curve with changing proportion of reference transcript abundances
		# unreliable, unreliable, reliable but contradictary, reliable and consistent, reliable and consistent
		tnames = ["ENST00000010404.6", "ENST00000396832.5", "ENST00000377482.9", "ENST00000451137.6", "ENST00000337393.9"]
		indexes = [transnames.index(tname) for tname in tnames]
		WriteExampleCurve(folder + "/MCF10/IValue_curve_example.txt", indexes, TreatmentSampleMap, transnames, Expression_pertreatment, Lb_pertreatment, Ub_pertreatment)
