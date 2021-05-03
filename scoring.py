from sequence import NucleicAcid, Profile, Sequence
from itertools import combinations

class Scoring_Config :
	# matrix is 2d array such that matrix[i][j] is the score of NucleicAcid[i] alligned with NucleicAcid[j]
	def __init__(self, matrix, open_gap, cont_gap) :
		self.matrix = matrix
		self.open_gap = open_gap
		self.cont_gap = cont_gap

	def __getitem__(self, index) :
		if not isinstance(index, tuple) and len(index) == 2 :
			raise ValueError("Index must be tuple of len 2")
		if not (isinstance(index[0], NucleicAcid) and isinstance(index[1], NucleicAcid)) :
			raise ValueError("Index must be of form (NucleicAcid, NucleicAcid)")
		if not (index[0].is_acid() and index[1].is_acid()) :
			raise ValueError("Index Nucleic acids must be acid types")
		return self.matrix[index[0].value][index[1].value]

blast_config = Scoring_Config([
	[ 5,  -4, -4, -4],
	[ -4,  5, -4, -4],
	[ -4, -4,  5, -4],
	[ -4, -4, -4, 5 ],
], 9, 1)

def get_score(_seqs, scoring_config=blast_config) :
	seqs = []
	for seq in _seqs :
		if isinstance(seq, Profile) :
			seqs += seq.seqs
		elif isinstance(seq, Sequence) :
			seqs.append(seq)
		else :
			assert False, f'input "{seq}"'
	score = 0
	for i in range(max([len(s) for s in seqs])) :
		total_gap_penalty = 0
		nas = []
		for seq in seqs :
			if i < len(seq) :
				na = seq[i]
			else :
				na = NucleicAcid.GAP
			if na == NucleicAcid.GAP :
				if i != 0 and (i-1 >= len(seq) or seq[i-1] == NucleicAcid.GAP):
					total_gap_penalty += scoring_config.cont_gap
				else :
					total_gap_penalty += scoring_config.open_gap
			else :
				nas.append(na)
		score -= total_gap_penalty
		nas = combinations(nas, 2) # generate all combintions of length 2
		for na_pair in nas :
			score += scoring_config[na_pair]
	return score
