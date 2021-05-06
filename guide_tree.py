from copy import copy
from os import remove
from sequence import NucleicAcid, Sequence
from scoring import get_score
from multiprocessing import Pool, cpu_count
from needleman_wunsch import needleman_wunsch_ss, print2darr
from copy import deepcopy
from timer import Timer
from itertools import product

"""
		Guide tree classes
"""

class Node :
		def __init__(self, id, seq=None, left=None, right=None) :
			# container variables
			self.seq = seq
			self.id = id
			self.left = left
			self.right = right
			self.up = None
			
			# memoization for methods variables
			self.sub_leafs_seqs = [self.seq] if self.is_leaf() else None
		
		def __repr__(self) :
			return f'<{self.id} : {None if self.left is None else self.left.id}, {None if self.right is None else self.right.id}, {None if self.up is None else self.up.id}>'

		def __str__(self) :
			if self.is_leaf() :
				return f"<LNode {repr(self.seq)}>"
			else :
				return f"<INode {'P' if self.seq is not None else ''}({str(self.left)[7:-1]}, {str(self.right)[7:-1]})>"


		def is_leaf(self) :
			if self.left is None and self.right is None and isinstance(self.seq, Sequence) :
				return True
			assert isinstance(self.left, Node) and isinstance(self.left, Node)
			return False
		
		# determines if subtrees begining at nodes have the same structure
		def __eq__(self, other) :
			assert isinstance(other, Node)
			if self.is_leaf() != other.is_leaf() :
				return False
			elif self.is_leaf() :
				if self.seq.name == other.seq.name :
					return True
				else :
					return False
			elif (self.left == other.left and self.right == other.right) \
					or (self.left == other.right and self.right == other.left) :
				return True
			else :
				return False

		# goes down subtree structure until a node with a sequences is found 
		# 	- subtrees to that node have external links, and so cannot have external gaps
		def remove_external_gaps(self) :
			if self.seq is not None : # remove and stop if node has seq

				self.seq.remove_external_gaps()
			elif not self.is_leaf() :
				self.left.remove_external_gaps()
				self.right.remove_external_gaps()
			else :
				assert False

		def start_pairs(self) :
			if self.seq is not None :
				return []
			elif not self.is_leaf() :
				if self.left.seq is not None and self.right.seq is not None :
					return [(self.left, self.right)]
				else :
					return self.left.start_pairs() + self.right.start_pairs()
			else :
				assert False

		def alignments_to_completion(self) :
			if self.seq is None :
				assert not self.is_leaf()
				return 1 + self.left.alignments_to_completion() + self.right.alignments_to_completion()
			else :
				return 0



class Tree :
	def __init__(self, tuple_tree) :
		self._str = str(tuple_tree)
		self.nodes = []
		self.start_pairs = []
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
				self.start_pairs.append((left, right))
		else :
			raise ValueError(f"tuple_tree contains invalid entry :\n\t{tuple_tree}")
		self.nodes.append(node)
		return node

	# returns the number of internal nodes whose subtrees changed
	def extract_reusable(self, old) :
		assert isinstance(old, Tree)
		self_internal_nodes  = [n for n in self.nodes  if not n.is_leaf()]
		old_internal_nodes = [n for n in old.nodes if not n.is_leaf()]
		# note: equality_mat[i][q] = self_internal_nodes[i] == old_internal_nodes[q]
		equality_mat = [[None for _ in range(len(old_internal_nodes))] for _ in range(len(self_internal_nodes))]
		pairs = []
		for (si, sn), (oi, on) in product(enumerate(self_internal_nodes), enumerate(old_internal_nodes)) :
			if Tree._enter_equality(self_internal_nodes, old_internal_nodes, equality_mat, si, oi) :
				pairs.append((sn, on))

		for self_node, old_node in pairs :
			self_node.seq = old_node.seq

		self.start_pairs = self.root.start_pairs()
		return len(self_internal_nodes) - len(pairs)

	# must be called before tree is iterated
	def remove_external_gaps(self) :
		self.root.remove_external_gaps()

	def alignments_to_completion(self) :
		return self.root.alignments_to_completion()

	def __str__(self) :
		return f"<TREE {str(self.root)[7:-1]}>"

	@staticmethod
	def _enter_equality(self_nodes, other_nodes, eq_mat, i, q) :
		if eq_mat[i][q] is not None :
			return eq_mat[i][q]
		self_node = self_nodes[i]
		other_node = other_nodes[q]
		res = None
		# compare children that are leafs
		self_child_leaf = sorted([n.seq for n in [self_node.left, self_node.right] if n.is_leaf()], key=lambda x: x.name)
		other_child_leaf = sorted([n.seq for n in [other_node.left, other_node.right] if n.is_leaf()], key=lambda x: x.name)
		if self_child_leaf != other_child_leaf :
			res = False
		elif len(self_child_leaf) == 2 :
			assert len(other_child_leaf) == 2
			res = True
		else :
			# compare children that are internals
			self_child_internals = [self_nodes.index(n) for n in [self_node.left, self_node.right] if not n.is_leaf()]
			other_child_internals = [other_nodes.index(n) for n in [other_node.left, other_node.right] if not n.is_leaf()]
			n = len(self_child_internals)
			assert n == len(other_child_internals) and n != 0
			if n == 1 :
				res = Tree._enter_equality(self_nodes, other_nodes, eq_mat, self_child_internals[0], other_child_internals[0])
			elif n == 2 :
				res = (Tree._enter_equality(self_nodes, other_nodes, eq_mat, self_child_internals[0], other_child_internals[0]) and \
						Tree._enter_equality(self_nodes, other_nodes, eq_mat, self_child_internals[1], other_child_internals[1])) or \
						(Tree._enter_equality(self_nodes, other_nodes, eq_mat, self_child_internals[0], other_child_internals[1]) and \
						Tree._enter_equality(self_nodes, other_nodes, eq_mat, self_child_internals[1], other_child_internals[0]))
			else : 
				assert False
		# handle result
		assert res is not None
		assert eq_mat[i][q] is None
		if res : # if a pair is found, neither element can be part of another pair
			for x in range(len(self_nodes)) :
				assert eq_mat[x][q] is not True
				eq_mat[x][q] = False
			for x in range(len(other_nodes)) :
				assert eq_mat[i][x] is not True
				eq_mat[i][x] = False
			eq_mat[i][q] = True
		else :
			eq_mat[i][q] = False
		return res




def print2darr(arr) :
	arr = [' '.join([(f'{i:.2e}' if isinstance(i, float) else str(i)).ljust(7) for i in a]) for a in arr]
	arr = ',\n '.join(arr)
	print(f'[{arr}]')

"""
		kmer distance functions
"""

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

"""
		needleman-wunsch distance functions
"""

def pairwise_distance_nw(seq1, seq2, i, q) :
	seq1 = deepcopy(seq1)
	seq2 = deepcopy(seq2)
	needleman_wunsch_ss(seq1, seq2)
	return (1/get_score(seq1, seq2), i, q)

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

"""
		score distance functions
			- Not to be used for initial guide tree construction, but may be used for iterative improvements
"""

def get_dmat_scores(seqs, scoring_config) :
	mat = []
	for _ in seqs :
		mat.append([None] * len(seqs))
	for i in range(len(seqs)) :
		mat[i][i] = 0
	for i in range(len(seqs)) :
		for q in range(i+1, len(seqs)) :
			dist = 1 / get_score(seqs[i], seqs[q])
			mat[i][q] = dist
			mat[q][i] = dist
	return mat


"""
		pmat distance
			- Not to be used for initial guide tree construction, but may be used for iterative improvements
"""

def pdist(seq1, seq2) :
	diffs = 0
	l = 0
	for i in range(max(len(seq1), len(seq2))) :
		na1 = seq1[i] if i < len(seq1) else NucleicAcid.GAP
		na2 = seq2[i] if i < len(seq2) else NucleicAcid.GAP
		if na1 != na2 :
			diffs += 1
			l += 1
		elif na1 != NucleicAcid.GAP or na2 != NucleicAcid.GAP :
			l += 1
	return diffs / l


def get_dmat_pdist(seqs) :
	mat = []
	for _ in seqs :
		mat.append([None] * len(seqs))
	for i in range(len(seqs)) :
		mat[i][i] = 0
	for i in range(len(seqs)) :
		for q in range(i+1, len(seqs)) :
			score = pdist(seqs[i], seqs[q])
			mat[i][q] = score
			mat[q][i] = score
	return mat



"""
		Neighbor joingin functions
"""

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
	seqs = copy(seqs) # don't mess up the original array
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
	
