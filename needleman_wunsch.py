from scoring import blast_config

def print2darr(arr) :
	arr = [' '.join([(f'{i:.2e}' if isinstance(i, float) else str(i)).ljust(7) for i in a]) for a in arr]
	arr = ',\n '.join(arr)
	print(f'[{arr}]')

DIAG_ARROW = 0
LEFT_ARROW = 1
UP_ARROW = 2

def needleman_wunsch(seq1, seq2, scoring_config=blast_config) :
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
	return smat[m-1][n-1]