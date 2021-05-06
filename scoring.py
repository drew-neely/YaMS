from sequence import NucleicAcid, Profile, Sequence
from itertools import combinations


def get_match_score(counts, scoring_config) :
	score = 0
	for i in range(4) :
		count = counts[i]
		if count > 0 :
			for q in range(i+1, 4) :
				count2 = counts[q]
				score += count * count2 * scoring_config.matrix[i][q]
			if count > 1 :
				score += (count * (count - 1) // 2) * scoring_config.matrix[i][i]
	return score

def get_gap_score(og, cg, scoring_config) :
	return - (og * scoring_config.open_gap + cg * scoring_config.cont_gap)

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

def flatten_args(*_seqs) :
	seqs = []
	for seq in _seqs :
		if isinstance(seq, Profile) :
			seqs += seq.seqs
		elif isinstance(seq, Sequence) :
			seqs.append(seq)
		elif isinstance(seq, list) or isinstance(seq, tuple):
			for e in seq :
				seqs += flatten_args(e)
		else :
			assert False, f'input type: {type(seq)} -\n"{seq}"'
	return seqs

def get_score(*_seqs, scoring_config=blast_config) :
	seqs = flatten_args(_seqs)
	prof = Profile(seqs)
	score = 0
	last_valid_col = -1
	for i in range(max([len(s) for s in seqs])) :
		nas, og, cg = prof[i]
		if og + cg == len(prof) : # only gaps at i - nothing to do
			assert nas == [0,0,0,0]
		else : # at least one not gap at i
			if i - 1 != last_valid_col : # had to skip columns
				og, cg = 0, 0
				for seq in seqs :
					if seq[i] == NucleicAcid.GAP and last_valid_col >= 0 and seq[last_valid_col] == NucleicAcid.GAP :
						cg += 1
					elif seq[i] == NucleicAcid.GAP :
						og += 1
			
			score += get_gap_score(og, cg, scoring_config) + get_match_score(nas, scoring_config)
			last_valid_col = i
	return score
