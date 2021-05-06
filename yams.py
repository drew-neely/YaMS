
from argparse import ArgumentParser
from sequence import Profile, read_fasta, Sequence
from scoring import get_score
from guide_tree import join_neighbors, Tree, get_dmat_kmers, get_dmat_nw, get_dmat_scores, get_dmat_pdist
from needleman_wunsch import needleman_wunsch_ss, needleman_wunsch_pp
from timer import Timer
from align import align_tree

import signal
from multiprocessing import cpu_count

signal.signal(signal.SIGINT, lambda x, y: exit(0))

parser = ArgumentParser(description='Beep boop align sequence')
parser.add_argument('infile', help='input file in .fasta format')
parser.add_argument('outfile', help='output file name', nargs='?')
parser.add_argument('-no-guide-tree', '-ngt', help='will use guide tree based on order rather than distance if specified', action="store_true")
parser.add_argument('-kmer-distance', '-kdist', metavar='k', help='k value to use in kmer distance calculation', type=int)
parser.add_argument('-score-distance', '-sdist', help='use needleman wunsch alignment to calculate distance', action='store_true')
parser.add_argument('-iterate-to-convergence', "-converge", help="iterate until convergence is reached", action='store_true')
parser.add_argument('-iterative-depth', '-depth', metavar='d', help='number of iterative steps to perform', type=int)
parser.add_argument('-no-iterative', '-noit', help='don\'t perform the iterative step at all', action='store_true')

parser.add_argument('-threads', '-t', metavar='t', help='Number of threads to use - default is the number of cpus avaliable', type=int, default=None)
args = parser.parse_args()

# print('Args=', args)
if args.kmer_distance != None and args.score_distance :
	print("Error: Cannot do both kmer distance and needleman wunsch distance calculations")
	exit()

if (args.iterate_to_convergence and args.iterative_depth is not None) or \
		(args.iterate_to_convergence and args.no_iterative) or \
		(args.iterative_depth is not None and args.no_iterative) :
	print("Error: Only one iterative depth option may be specified")
	exit()

if args.outfile is None :
		print("Warning: no output file specified - result will be thrown away")

if args.iterate_to_convergence :
	args.iterative_depth = -1
elif args.no_iterative :
	args.iterative_depth = 0
elif args.iterative_depth is None:
	print("Warning: No iterative dpeth specified - defaulting to iterate until convergence")
	args.iterative_depth = -1
elif args.iterative_depth < -1 :
	print(f"Error: Invalid iterative depth: {args.iterative_depth}")
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
		Iteratively improve
	"""
	changed_nodes = None
	for d in range(args.iterative_depth if args.iterative_depth >= 0 else 10000) :
		print("------------")
		dmat = get_dmat_pdist(sequences)
		new_tree = join_neighbors(sequences, dmat)
		new_changed_nodes = new_tree.extract_reusable(tree)
		print("changed_nodes", changed_nodes, new_changed_nodes)
		if new_changed_nodes == 0 :
			print("Iteration converged via guide tree convergence")
			break
		elif changed_nodes is not None and new_changed_nodes >= changed_nodes: 
			print("Iteration converged via lack of progression in guide tree")
			break
		else :
			print(f"beginning iterative step {d}")
			changed_nodes = new_changed_nodes
			new_tree.remove_external_gaps()
			sequences = align_tree(new_tree, threads=args.threads)
			tree = new_tree
	
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
	if args.outfile is not None :
		output = '\n\n'.join([str(seq) for seq in sequences])
		with open(args.outfile, 'w') as file:
			file.write(output)