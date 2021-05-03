from sequence import NucleicAcid
from scoring import Scoring_Config, blast_config
from math import factorial
from timer import Timer

def print2darr(arr) :
	arr = [' '.join([(f'{i:.2e}' if isinstance(i, float) else str(i)).ljust(7) for i in a]) for a in arr]
	arr = ',\n '.join(arr)
	print(f'[{arr}]')

# def nCr(n,r):
#     return factorial(n) // factorial(r) // factorial(n-r)

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


DIAG_ARROW = 0
LEFT_ARROW = 1
UP_ARROW = 2

def needleman_wunsch_ss(seq1, seq2, scoring_config=blast_config) :
	Timer.start("needleman_wunsch_ss")
	n, m = len(seq1)+1, len(seq2)+1 # matrices are m x n (m rows, n columns)
	smat = [None] * (n * m) # (columns are elements of seq1, row are elements of seq2)
	gmat = [None] * (n * m)
	hmat = [None] * (n * m)
	dmat = [None] * (n * m) # direction of selection

	"""
			Initialize first row and column
	"""
	smat[0] = 0
	for row in range(1, m) :
		score = - scoring_config.open_gap - (row - 1) * scoring_config.cont_gap
		smat[row * n] = score
		gmat[row * n] = score
		hmat[row * n] = score
		dmat[row * n] = UP_ARROW
	for col in range(1, n) :
		score = - scoring_config.open_gap - (col - 1) * scoring_config.cont_gap
		smat[col] = score
		gmat[col] = score
		hmat[col] = score
		dmat[col] = LEFT_ARROW
	
	"""
		Fill in tables and arrows
	"""
	for row in range(1, m) :
		for col in range(1, n) :
			# Fill in matricies at (row, col) - It is guranteed tht (row-1,col),(row-1,col-1), and (ro1, col-1) exist			

			diag_idx = (row-1)*n + (col-1)
			left_idx = row * n + (col-1)
			up_idx = (row - 1) * n + col
			row_col_idx = row * n + col

			f = smat[diag_idx] + scoring_config[(seq1[col-1], seq2[row-1])]
			gmat[row_col_idx] = max(smat[up_idx] - scoring_config.open_gap,
									gmat[up_idx] - scoring_config.cont_gap)
			hmat[row_col_idx] = max(smat[left_idx] - scoring_config.open_gap,
									hmat[left_idx] - scoring_config.cont_gap)
			score = max(f, gmat[row_col_idx], hmat[row_col_idx])
			smat[row_col_idx] = score
			if score == f:
				dmat[row_col_idx] = DIAG_ARROW
			elif score == gmat[row_col_idx] :
				dmat[row_col_idx] = UP_ARROW
			elif score == hmat[row_col_idx] :
				dmat[row_col_idx] = LEFT_ARROW
			else :
				assert(False)
	"""
		Parse arrows to add gaps
	"""
	row, col = m - 1, n - 1
	while row > 0 or col > 0 :
		idx = row * n + col
		if dmat[idx] == DIAG_ARROW :
			row -= 1
			col -= 1
		elif dmat[idx] == UP_ARROW :
			seq1.insert_gap(col)
			row -= 1
		elif dmat[idx] == LEFT_ARROW :
			seq2.insert_gap(row)
			col -= 1
		else :
			assert(False)
	"""
		Return score of alignment
	"""
	Timer.end("needleman_wunsch_ss")
	return smat[-1]

# from bottom right corner :
# 	Move left represents adding in

def needleman_wunsch_pp(prof1, prof2, scoring_config=blast_config) :
	Timer.start("needleman_wunsch_pp")
	n, m = prof1.seq_len+1, prof2.seq_len+1 # matrices are m x n (m rows, n columns)
	smat = [None] * (n * m) # (columns are elements of seq1, row are elements of seq2)
	gmat = [None] * (n * m)
	hmat = [None] * (n * m)
	dmat = [None] * (n * m) # direction of selection

	"""
			Initialize first row and column
	"""
	smat[0] = 0
	(nas1_0, og1_0, _) = prof1[0]
	(nas2_0, og2_0, _) = prof2[0]
	# first row
	score1r = (- scoring_config.open_gap) * (len(prof1) - og1_0) + (- scoring_config.cont_gap) * og1_0 \
				+ get_match_score(nas2_0, scoring_config) + (- scoring_config.open_gap) * og2_0
	smat[n] = score1r
	gmat[n] = score1r
	hmat[n] = score1r
	dmat[n] = UP_ARROW
	
	# first col
	score1c = (- scoring_config.open_gap) * (len(prof2) - og2_0) + (- scoring_config.cont_gap) * og2_0 \
				+ get_match_score(nas1_0, scoring_config) + (- scoring_config.open_gap) * og1_0
	smat[1] = score1c
	gmat[1] = score1c
	hmat[1] = score1c
	dmat[1] = LEFT_ARROW

	for row in range(2, m) : # inserting gap in prof1
		(nas2, og2, cg2) = prof2[row-1]
		score = smat[(row-1) * n] + (- scoring_config.cont_gap) * len(prof1) \
				+ get_match_score(nas2, scoring_config) + (- scoring_config.open_gap) * og2 \
				+ (- scoring_config.cont_gap) * cg2
		smat[row * n] = score
		gmat[row * n] = score
		hmat[row * n] = score
		dmat[row * n] = UP_ARROW
	for col in range(2, n) : # inserting a gap in prof2
		(nas1, og1, cg1) = prof1[col-1]
		score = smat[col-1] + (- scoring_config.cont_gap) * len(prof2) \
				+ get_match_score(nas1, scoring_config) + (- scoring_config.open_gap) * og1 \
				+ (- scoring_config.cont_gap) * cg1
		smat[col] = score
		gmat[col] = score
		hmat[col] = score
		dmat[col] = LEFT_ARROW

	"""
			Initialize reference data
	"""
	# store index calculations
	index_info_1 = []
	index_info_2 = []
	for col in range(1, n) :
		index_info_1.append(prof1[col-1])
	for row in range(1, m) :
		index_info_2.append(prof2[row-1])

	# count types of gaps for insertion at each index
	iogs1, icgs1, iogs2, icgs2 = [], [], [], [] # Insert Open Gap, Insert Cont Gap 
	i_gap_scores_1, i_gap_scores_2 = [], []
	for col in range(1, n) :
		icgs1.append(0)
		iogs1.append(0)
		for seq in prof1.seqs :
			if seq[col-1] == NucleicAcid.GAP or (col+1 < len(seq) and seq[col+1] == NucleicAcid.GAP) :
				icgs1[col-1] += 1
			else :
				iogs1[col-1] += 1
		i_gap_scores_1.append(get_gap_score(iogs1[col-1], icgs1[col-1], scoring_config))
	for row in range(1, m) :
		icgs2.append(0)
		iogs2.append(0)
		for seq in prof2.seqs :
			if seq[row-1] == NucleicAcid.GAP or (row+1 < len(seq) and seq[row+1] == NucleicAcid.GAP) :
				icgs2[row-1] += 1
			else :
				iogs2[row-1] += 1
		i_gap_scores_2.append(get_gap_score(iogs2[row-1], icgs2[row-1], scoring_config))

	# count match_scores
	match_scores_1 = []
	match_scores_2 = []
	gap_scores_1 = []
	gap_scores_2 = []
	for col in range(1, n) :
		(nas1, og1, cg1) = index_info_1[col-1]
		match_scores_1.append(get_match_score(nas1, scoring_config))
		gap_scores_1.append(get_gap_score(og1, cg1, scoring_config))
	for row in range(1, m) :
		(nas2, og2, cg2) = index_info_2[row-1]
		match_scores_2.append(get_match_score(nas2, scoring_config))
		gap_scores_2.append(get_gap_score(og2, cg2, scoring_config))
	
	# get match comparison scores for all pairs of indecies
	match_scores_12 = [None] * (m * n)
	nas = [0] * 4
	for row in range(1, m) :
		for col in range(1, n) :
			nas1 = index_info_1[col-1][0]
			nas2 = index_info_2[row-1][0]
			for i in range(4) :
				nas[i] = nas1[i] + nas2[i]
			match_scores_12[row * n + col] = get_match_score(nas, scoring_config)

	"""
		Fill in tables and arrows
	"""
	for row in range(1, m) :
		for col in range(1, n) :
			diag_idx = (row-1)*n + (col-1)
			left_idx = row * n + (col-1)
			up_idx = (row - 1) * n + col
			row_col_idx = row * n + col

			# Fill in matricies at (row, col) - It is guranteed tht (row-1,col),(row-1,col-1), and (ro1, col-1) exist
			(nas1, og1, cg1), (nas2, og2, cg2) = index_info_1[col-1], index_info_2[row-1]

			# calculate reusable variables
			match_score_12 = match_scores_12[row_col_idx]
			match_score_1  = match_scores_1[col-1]
			match_score_2  = match_scores_2[row-1]

			gap_score_12 = get_gap_score(og1 + og2, cg1 + cg2, scoring_config)
			gap_score_1  = gap_scores_1[col-1]
			gap_score_2  = gap_scores_2[row-1]

			i_gap_score_1  = i_gap_scores_1[col-1]
			i_gap_score_2  = i_gap_scores_2[row-1]

			# match prof1 to prof2
			f = smat[diag_idx] + match_score_12 + gap_score_12 # DIAG_ARROW

			# gap in prof1
			gs = smat[up_idx] + match_score_2 + gap_score_2 + i_gap_score_1
			gg = gmat[up_idx] + match_score_2 + gap_score_2 - (scoring_config.cont_gap if row-1 != 0 else scoring_config.open_gap) * len(prof1)
			gmat[row_col_idx] = max(gs, gg) # UP_ARROW

			# gap in prof2
			hs = smat[left_idx] + match_score_1 + gap_score_1 + i_gap_score_2
			hh = hmat[left_idx] + match_score_1 + gap_score_1 - (scoring_config.cont_gap if col-1 != 0 else scoring_config.open_gap) * len(prof2)
			hmat[row_col_idx] = max(hs, hh) # LEFT_ARROW

			# print(f'\tf: {f}, gs {gs}, gg {gg}, hs {hs}, hh {hh}')
			score = max(f, gmat[row_col_idx], hmat[row_col_idx])
			smat[row_col_idx] = score
			if score == f:
				dmat[row_col_idx] = DIAG_ARROW
				# print("\tDIAG")
			elif score == gmat[row_col_idx] :
				dmat[row_col_idx] = UP_ARROW
				# print("\tUP")
			elif score == hmat[row_col_idx] :
				dmat[row_col_idx] = LEFT_ARROW
				# print("\tLEFT")
			else :
				assert(False)
	# print2darr(smat)
	# print()
	# print2darr(dmat)
	# print()

	"""
		Parse arrows to add gaps
	"""
	row, col = m - 1, n - 1
	# print("PATH")
	while row > 0 or col > 0 :
		# print((row - 1, col - 1, smat[row][col], None if dmat[row][col] is None else ['DIAG', 'LEFT', 'UP'][dmat[row][col]]))
		idx = row * n + col
		if dmat[idx] == DIAG_ARROW :
			row -= 1
			col -= 1
		elif dmat[idx] == UP_ARROW :
			prof1.insert_gap(col)
			row -= 1
		elif dmat[idx] == LEFT_ARROW :
			prof2.insert_gap(row)
			col -= 1
		else :
			assert(False)
	"""
		Return score of alignment
	"""
	Timer.end("needleman_wunsch_pp")
	return smat[-1]