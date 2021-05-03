from copy import copy
from sequence import Sequence
from scoring import get_score
from multiprocessing import Pool, cpu_count
from needleman_wunsch import needleman_wunsch_ss
from copy import deepcopy
from timer import Timer

class Node :
		def __init__(self, id, seq=None, left=None, right=None) :
			self.seq = seq
			self.id = id
			self.left = left
			self.right = right
			self.up = None
		
		def __repr__(self) :
			return f'<{self.id} : {None if self.left is None else self.left.id}, {None if self.right is None else self.right.id}, {None if self.up is None else self.up.id}>'

class Tree :
	def __init__(self, tuple_tree) :
		self.nodes = []
		self.leaf_pairs = []
		self.root = self.build_tree(tuple_tree)

	def build_tree(self, tuple_tree) :
		if isinstance(tuple_tree, Sequence) :
			node = Node(len(self.nodes), seq=tuple_tree)
		elif isinstance(tuple_tree, tuple) and len(tuple_tree) == 2:
			left = self.build_tree(tuple_tree[0])
			right = self.build_tree(tuple_tree[1])
			node = Node(len(self.nodes), left=left, right=right)
			left.up = node
			right.up = node
			if isinstance(tuple_tree[0], Sequence) and isinstance(tuple_tree[1], Sequence) :
				self.leaf_pairs.append((left, right))
		else :
			raise ValueError(f"tuple_tree contains invalid entry :\n\t{tuple_tree}")
		self.nodes.append(node)
		return node

def print2darr(arr) :
	arr = [', '.join(['None    ' if i is None else f'{i:.2e}' for i in a]) for a in arr]
	arr = ',\n '.join(arr)
	print(f'[{arr}]')

def subseqs(seq, n=None) :
	if n is None :
		n = 2
	sub_counts = {}
	for i in range(len(seq) - n + 1) :
		seq_str = ''.join([str(na) for na in seq[i:i+n]])
		if seq_str in sub_counts :
			sub_counts[seq_str] += 1
		else :
			sub_counts[seq_str] = 1
	return sub_counts

def pairwise_distance_kmer(seq1, seq2, k) :
	if k is None :
		k = 2
	closeness = 0
	for _ in range(2, k+1) :
		sub1 = subseqs(seq1, k)
		sub2 = subseqs(seq2, k)
		for subseq in set(sub1.keys()).intersection(set(sub2.keys())) :
			closeness += min(sub1[subseq], sub2[subseq])
	return 1 / closeness

# get distance matrix with kmer distances
def get_dmat_kmers(seqs, k) :
	mat = []
	for _ in seqs :
		mat.append([None] * len(seqs))
	for i, seq1 in enumerate(seqs) :
		for q, seq2 in enumerate(seqs) :
			if q == i :
				mat[i][q] = 0
			elif q > i :
				dist = pairwise_distance_kmer(seq1, seq2, k)
				mat[i][q] = dist
				mat[q][i] = dist
	return mat

def pairwise_distance_nw(seq1, seq2, i, q) :
	seq1 = deepcopy(seq1)
	seq2 = deepcopy(seq2)
	needleman_wunsch_ss(seq1, seq2)
	return (get_score(seq1, seq2), i, q)

def get_dmat_nw(seqs, threads) :
	mat = []
	for _ in seqs :
		mat.append([None] * len(seqs))
	for i in range(len(seqs)) :
		mat[i][i] = 0
	if threads == None :
		threads = cpu_count()
	pairs = []
	for i in range(len(seqs)) :
		for q in range(i+1, len(seqs)) :
			pairs.append((seqs[i], seqs[q], i, q))
	with Pool(processes=threads) as pool :
		dists = pool.starmap(pairwise_distance_nw, pairs)
	for dist, i, q in dists :
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
def join_neighbors(seqs, dmat) :
	assert dmat is not None, "join neighbors must be given a distance matrix"
	assert len(seqs) > 1, "Can't do neighbor joining on only one element"
	seqs = copy(seqs) # don't mes up the original array
	while len(seqs) > 2 :
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
	seqs = (seqs[0], seqs[1]) # no need to do last iteration with only 2 things in the list
	return Tree(seqs)
	
