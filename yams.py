from argparse import ArgumentParser
from sequence import Profile, read_fasta, Sequence
from scoring import get_score
from guide_tree import join_neighbors, Tree, get_dmat_kmers, get_dmat_nw
from needleman_wunsch import needleman_wunsch_ss, needleman_wunsch_pp
from timer import Timer
from align import align_tree
import signal

signal.signal(signal.SIGINT, lambda x, y: exit(0))

parser = ArgumentParser(description='Beep boop align sequence')
parser.add_argument('infile', help='input file in .fasta format')
parser.add_argument('outfile', help='output file name')
parser.add_argument('-no-guide-tree', '-ngt', help='will use guide tree based on order rather than distance if specified', action="store_true")
parser.add_argument('-kmer-distance', '-kdist', help='k value to use in kmer distance calculation', type=int, default=None)
parser.add_argument('-score-distance', '-sdist', help='use needleman wunsch alignment to calculate distance', action='store_true')

parser.add_argument('-threads', '-t', help='Number of threads to use - default is the number of cpus avaliable', type=int, default=None)
args = parser.parse_args()

# print('Args=', args)
if args.kmer_distance != None and args.score_distance :
	print("Cannot do both kmer distance and needleman wunsch distance calculations")
	exit()



if __name__ == "__main__" :
	"""
		Read sequences from file
	"""
	sequences = read_fasta(args.infile)


	"""
		Print initial info
	"""
	starting_score = get_score(sequences)
	print(f'Starting score: {starting_score:,}')
	print('Begining Progressive MSA')
	print('------------\n')
	Timer.start("alignment")
	

	"""
		Calculate distance matrix (dmat)
	"""
	dmat = None
	if not args.no_guide_tree : # don't waste time on dmat if we aren't gonna use it
		if args.score_distance :
			dmat = get_dmat_nw(sequences, threads=args.threads)
		else :
			dmat = get_dmat_kmers(sequences, args.kmer_distance)


	"""
		build guide tree
	"""
	tree = None
	if args.no_guide_tree :
		tree = sequences
		while len(tree) > 1 :
			for i in range(len(tree) // 2) :
				tree[i] = (tree[i], tree.pop(i+1))
		tree = Tree(tree[0])
	else :
		tree = join_neighbors(sequences, dmat)


	"""
		Do the alignment
	"""
	sequences = align_tree(tree, threads=args.threads)

	
	"""
		Print results
	"""
	Timer.end("alignment")
	print('------------------------------------------------')
	print(f'score: {starting_score:,} -> {get_score(sequences):,}')
	print(f'time_elapsed: {Timer.log("alignment", log=False):5.2f}')

	"""
		Output result to file
	"""
	output = '\n\n'.join([str(seq) for seq in sequences])
	with open(args.outfile, 'w') as file:
		file.write(output)