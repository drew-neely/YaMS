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
	Timer.start("get_match_score")
	score = 0
	counts = list(counts.items())
	for i in range(len(counts)) :
		na, count = counts[i]
		if count > 1 :
			score += (count * (count - 1) // 2) * scoring_config[(na, na)]
		for q in range(i+1, len(counts)) :
			na2, count2 = counts[q]
			score += count * count2 * scoring_config[(na, na2)]
	Timer.end("get_match_score")
	return score

def get_gap_score(og, cg, scoring_config) :
	return - (og * scoring_config.open_gap + cg * scoring_config.cont_gap)


DIAG_ARROW = 0
LEFT_ARROW = 1
UP_ARROW = 2

def needleman_wunsch_ss(seq1, seq2, scoring_config=blast_config) :
	Timer.start("needleman_wunsch_ss")
	n, m = len(seq1)+1, len(seq2)+1 # matrices are m x n (m rows, n columns)
	smat = [] # (columns are elements of seq1, row are elements of seq2)
	gmat = []
	hmat = []
	dmat = [] # direction of selection
	for _ in range(m) :
		smat.append([None] * n)
		gmat.append([None] * n)
		hmat.append([None] * n)
		dmat.append([None] * n)

	"""
			Initialize first row and column
	"""
	smat[0][0] = 0
	for row in range(1, m) :
		score = - scoring_config.open_gap - (row - 1) * scoring_config.cont_gap
		smat[row][0] = score
		gmat[row][0] = score
		hmat[row][0] = score
		dmat[row][0] = UP_ARROW
	for col in range(1, n) :
		score = - scoring_config.open_gap - (col - 1) * scoring_config.cont_gap
		smat[0][col] = score
		gmat[0][col] = score
		hmat[0][col] = score
		dmat[0][col] = LEFT_ARROW
	
	"""
		Fill in tables and arrows
	"""
	starting_diags = list(zip(range(1,m), [1]*(m-1))) + list(zip([m-1]*(n-2), range(2,n)))
	for row, col in starting_diags :
		while row > 0 and col < n :
			# Fill in matricies at (row, col) - It is guranteed tht (row-1,col),(row-1,col-1), and (ro1, col-1) exist			
			f = smat[row-1][col-1] + scoring_config[(seq1[col-1], seq2[row-1])]
			gmat[row][col] = max(smat[row-1][col] - scoring_config.open_gap,
									gmat[row-1][col] - scoring_config.cont_gap)
			hmat[row][col] = max(smat[row][col-1] - scoring_config.open_gap,
									hmat[row][col-1] - scoring_config.cont_gap)
			score = max(f, gmat[row][col], hmat[row][col])
			smat[row][col] = score
			if score == f:
				dmat[row][col] = DIAG_ARROW
			elif score == gmat[row][col] :
				dmat[row][col] = UP_ARROW
			elif score == hmat[row][col] :
				dmat[row][col] = LEFT_ARROW
			else :
				assert(False)
			row, col = row - 1, col + 1
	"""
		Parse arrows to add gaps
	"""
	row, col = m - 1, n - 1
	while row > 0 or col > 0 :
		if dmat[row][col] == DIAG_ARROW :
			row -= 1
			col -= 1
		elif dmat[row][col] == UP_ARROW :
			seq1.insert_gap(col)
			row -= 1
		elif dmat[row][col] == LEFT_ARROW :
			seq2.insert_gap(row)
			col -= 1
		else :
			assert(False)
	"""
		Return score of alignment
	"""
	Timer.end("needleman_wunsch_ss")
	return smat[m-1][n-1]

# from bottom right corner :
# 	Move left represents adding in

def needleman_wunsch_pp(prof1, prof2, scoring_config=blast_config) :
	Timer.start("needleman_wunsch_pp")
	n, m = prof1.seq_len+1, prof2.seq_len+1 # matrices are m x n (m rows, n columns)
	smat = [] # (columns are elements of seq1, row are elements of seq2)
	gmat = []
	hmat = []
	dmat = [] # direction of selection
	for _ in range(m) :
		smat.append([None] * n)
		gmat.append([None] * n)
		hmat.append([None] * n)
		dmat.append([None] * n)

	"""
			Initialize first row and column
	"""
	smat[0][0] = 0
	(nas1_0, og1_0, _) = prof1[0]
	(nas2_0, og2_0, _) = prof2[0]
	# first row
	score1r = (- scoring_config.open_gap) * (len(prof1) - og1_0) + (- scoring_config.cont_gap) * og1_0 \
				+ get_match_score(nas2_0, scoring_config) + (- scoring_config.open_gap) * og2_0
	smat[1][0] = score1r
	gmat[1][0] = score1r
	hmat[1][0] = score1r
	dmat[1][0] = UP_ARROW
	
	# first col
	score1c = (- scoring_config.open_gap) * (len(prof2) - og2_0) + (- scoring_config.cont_gap) * og2_0 \
				+ get_match_score(nas1_0, scoring_config) + (- scoring_config.open_gap) * og1_0
	smat[0][1] = score1c
	gmat[0][1] = score1c
	hmat[0][1] = score1c
	dmat[0][1] = LEFT_ARROW

	for row in range(2, m) : # inserting gap in prof1
		(nas2, og2, cg2) = prof2[row-1]
		score = smat[row-1][0] + (- scoring_config.cont_gap) * len(prof1) \
				+ get_match_score(nas2, scoring_config) + (- scoring_config.open_gap) * og2 \
				+ (- scoring_config.cont_gap) * cg2
		smat[row][0] = score
		gmat[row][0] = score
		hmat[row][0] = score
		dmat[row][0] = UP_ARROW
	for col in range(2, n) : # inserting a gap in prof2
		(nas1, og1, cg1) = prof1[col-1]
		score = smat[0][col-1] + (- scoring_config.cont_gap) * len(prof2) \
				+ get_match_score(nas1, scoring_config) + (- scoring_config.open_gap) * og1 \
				+ (- scoring_config.cont_gap) * cg1
		smat[0][col] = score
		gmat[0][col] = score
		hmat[0][col] = score
		dmat[0][col] = LEFT_ARROW

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
		(nas1, og1, cg1) = prof1[col-1]
		match_scores_1.append(get_match_score(nas1, scoring_config))
		gap_scores_1.append(get_gap_score(og1, cg1, scoring_config))
	for row in range(1, m) :
		(nas2, og2, cg2) = prof2[row-1]
		match_scores_2.append(get_match_score(nas2, scoring_config))
		gap_scores_2.append(get_gap_score(og2, cg2, scoring_config))
	
	
	"""
		Fill in tables and arrows
	"""
	starting_diags = list(zip(range(1,m), [1]*(m-1))) + list(zip([m-1]*(n-2), range(2,n)))
	for row, col in starting_diags :
		while row > 0 and col < n :
			# Fill in matricies at (row, col) - It is guranteed tht (row-1,col),(row-1,col-1), and (ro1, col-1) exist
			# if row == col  and row % 100 == 0:
			# 	print(row)
			# print(f'({row-1}, {col-1})')
			(nas1, og1, cg1), (nas2, og2, cg2) = index_info_1[col-1], index_info_2[row-1]
			# print(f"\tacids_1: {nas1}, og1: {og1}, cg1: {cg1}")
			# print(f"\tacids_2: {nas2}, og2: {og2}, cg1: {cg2}")
			nas = nas1 + nas2
			iog1, icg1 = iogs1[col-1], icgs1[col-1]
			iog2, icg2 = iogs2[row-1], icgs2[row-1]

			# iog1, icg1, iog2, icg2 = 0, 0, 0, 0 # Insert Open Gap, Insert Cont Gap 
			# for seq in prof1.seqs :
			# 	if seq[col-1] == NucleicAcid.GAP or (col+1 < len(seq) and seq[col+1] == NucleicAcid.GAP) :
			# 		icg1 += 1
			# 	else :
			# 		iog1 += 1
			# for seq in prof2.seqs :
			# 	if seq[row-1] == NucleicAcid.GAP or (row+1 < len(seq) and seq[row+1] == NucleicAcid.GAP) :
			# 		icg2 += 1
			# 	else :
			# 		iog2 += 1

			# print(f'\topen1: {iog1}, cont1: {icg1}')
			# print(f'\topen2: {iog2}, cont2: {icg2}')
			# print(f'\tmatch_score_all: {get_match_score(nas1 + nas2, scoring_config)}, gap_score_all: {get_gap_score(og1+og2, cg1+cg2, scoring_config)}')
			# print(f'\tmatch_score_1: {get_match_score(nas1, scoring_config)}, gap_score_1: {get_gap_score(og1, cg1, scoring_config)}')
			# print(f'\tmatch_score_2: {get_match_score(nas2, scoring_config)}, gap_score_2: {get_gap_score(og2, cg2, scoring_config)}')
			# print(f'\tsmat[row-1][col-1]: {smat[row-1][col-1]}')
			# print(f'\tsmat[row-1][col]: {smat[row-1][col]}, gmat[row-1][col]: {gmat[row-1][col]}')
			# print(f'\tsmat[row][col-1]: {smat[row-1][col-1]}, hmat[row][col-1]: {hmat[row][col-1]}')


			# calculate reusable variables
			match_score_12 = get_match_score(nas, scoring_config)
			match_score_1  = match_scores_1[col-1]
			match_score_2  = match_scores_2[row-1]

			gap_score_12 = get_gap_score(og1 + og2, cg1 + cg2, scoring_config)
			gap_score_1  = gap_scores_1[col-1]
			gap_score_2  = gap_scores_2[row-1]

			i_gap_score_1  = i_gap_scores_1[col-1]
			i_gap_score_2  = i_gap_scores_2[row-1]

			# match prof1 to prof2
			f = smat[row-1][col-1] + match_score_12  + gap_score_12 # DIAG_ARROW

			# gap in prof1
			gs = smat[row-1][col] + match_score_2 + gap_score_2 + i_gap_score_1
			gg = gmat[row-1][col] + match_score_2 + gap_score_2 - (scoring_config.cont_gap if row-1 != 0 else scoring_config.open_gap) * len(prof1)
			gmat[row][col] = max(gs, gg) # UP_ARROW

			# gap in prof2
			hs = smat[row][col-1] + match_score_1 + gap_score_1 + i_gap_score_2
			hh = hmat[row][col-1] + match_score_1 + gap_score_1 - (scoring_config.cont_gap if col-1 != 0 else scoring_config.open_gap) * len(prof2)
			hmat[row][col] = max(hs, hh) # LEFT_ARROW

			# print(f'\tf: {f}, gs {gs}, gg {gg}, hs {hs}, hh {hh}')
			score = max(f, gmat[row][col], hmat[row][col])
			smat[row][col] = score
			if score == f:
				dmat[row][col] = DIAG_ARROW
				# print("\tDIAG")
			elif score == gmat[row][col] :
				dmat[row][col] = UP_ARROW
				# print("\tUP")
			elif score == hmat[row][col] :
				dmat[row][col] = LEFT_ARROW
				# print("\tLEFT")
			else :
				assert(False)
			row, col = row - 1, col + 1
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
		if dmat[row][col] == DIAG_ARROW :
			row -= 1
			col -= 1
		elif dmat[row][col] == UP_ARROW :
			prof1.insert_gap(col)
			row -= 1
		elif dmat[row][col] == LEFT_ARROW :
			prof2.insert_gap(row)
			col -= 1
		else :
			assert(False)
	"""
		Return score of alignment
	"""
	Timer.end("needleman_wunsch_pp")
	return smat[m-1][n-1]