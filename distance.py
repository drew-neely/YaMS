
def print2darr(arr) :
	arr = [', '.join(['None    ' if i is None else f'{i:.2e}' for i in a]) for a in arr]
	arr = ',\n '.join(arr)
	print(f'[{arr}]')

def subseqs(seq, n=2) :
	sub_counts = {}
	for i in range(len(seq) - n + 1) :
		seq_str = ''.join([str(na) for na in seq[i:i+n]])
		if seq_str in sub_counts :
			sub_counts[seq_str] += 1
		else :
			sub_counts[seq_str] = 1
	return sub_counts

def pairwise_distance(seq1, seq2) :
	sub1 = subseqs(seq1)
	sub2 = subseqs(seq2)
	closeness = 0
	for subseq in set(sub1.keys()).intersection(set(sub2.keys())) :
		closeness += min(sub1[subseq], sub2[subseq])
	return 1 / closeness

# get distance matrix
def get_dmat(seqs) :
	mat = []
	for _ in seqs :
		mat.append([None] * len(seqs))
	for i, seq1 in enumerate(seqs) :
		for q, seq2 in enumerate(seqs) :
			if q == i :
				mat[i][q] = 0
			elif q > i :
				dist = pairwise_distance(seq1, seq2)
				mat[i][q] = dist
				mat[q][i] = dist
	return mat

# get q matrix
def get_qmat(dmat) :
	n = len(dmat)
	qmat = []
	for _ in range(n) :
		qmat.append([None] * n)
	for i in range(n) :
		for j in range(n) :
			qmat[i][j] = (n - 2) * dmat[i][j]
			for k in range(n) :
				qmat[i][j] -= dmat[i][k]
				qmat[i][j] -= dmat[j][k]
	return qmat

# imlements neigbor joining as described at :
# 	wikipedia.org/wiki/Neighbor_joining
def join_neighbors(seqs) :
	dmat = get_dmat(seqs)
	print2darr(dmat)
	while len(seqs) > 3 :
		qmat = get_qmat(dmat)
		min_val, min_pair = None, None
		for i in range(len(qmat)) :
			for j in range(len(qmat)) :
				if i != j and (min_val is None or qmat[i][j] < min_val) :
					min_val, min_pair = qmat[i][j], (i, j)
		assert min_pair is not None
		min_pair = tuple(sorted(min_pair, reverse=True))
		seqs.append(tuple([seqs.pop(a) for a in min_pair]))
		for row in dmat :
			row.append(None)
		dmat.append([None] * (len(seqs) + 2))
		for k in range(len(dmat) - 1) :
			if k not in min_pair :
				dist = 0.5 * (dmat[min_pair[0]][k] + dmat[min_pair[1]][k] - dmat[min_pair[0]][min_pair[1]])
				dmat[len(dmat)-1][k] = dist
				dmat[k][len(dmat)-1] = dist
		dmat[len(dmat)-1][len(dmat)-1] = 0
		for i in min_pair :
			dmat.pop(i)
		for row in dmat :
			for i in min_pair :
				row.pop(i)
	return seqs
	
