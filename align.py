import argparse
from distance import pairwise_distance
from sequence import Profile, read_fasta, Sequence
from scoring import get_score
from distance import pairwise_distance, join_neighbors
from needleman_wunsch import needleman_wunsch_ss, needleman_wunsch_pp
from timer import Timer

parser = argparse.ArgumentParser(description='Beep boop align sequence')
parser.add_argument('infile', help='input file in .fasta format')
parser.add_argument('outfile', help='output file name')
parser.add_argument('-genetic', '-g', action='store_true', 
	help='align with genetic algorithm')
parser.add_argument('-nw_pp_test', action='store_true',
	help='run profile profile alignment test')
args = parser.parse_args()
# print('Args=', args)

sequences = read_fasta(args.infile)


def nw_pp_test(s1, s2) :
	assert s1 < s2 # pop in backwards index order for simplicity
	s2 = sequences.pop(s2)
	s1 = sequences.pop(s1)
	alignment_function = None
	if isinstance(s1, Sequence) and isinstance(s2, Sequence) :
		# s1 = Profile([s1])
		# s2 = Profile([s2])
		alignment_function = needleman_wunsch_ss
	elif isinstance(s1, Profile) and isinstance(s2, Sequence) :
		s2 = Profile([s2])
		alignment_function = needleman_wunsch_pp
	elif isinstance(s1, Sequence) and isinstance(s2, Profile) :
		s1 = Profile([s1])
		alignment_function = needleman_wunsch_pp
	elif isinstance(s1, Profile) and isinstance(s2, Profile) :
		alignment_function = needleman_wunsch_pp
	else :
		assert(False)

	print(f"Performing alignment between {repr(s1)} and {repr(s2)}")
	print(f"Starting score: {get_score([s1, s2])}")
	# print(s1)
	# print(s2)
	nw_score = alignment_function(s1, s2)
	score = get_score([s1, s2])
	print(f"nw_score = {nw_score}, score = {score}")

	Timer.log(alignment_function.__name__)
	if alignment_function == needleman_wunsch_pp :
		Timer.log("get_match_score")

	# print(s1)
	# print(s2)

	if nw_score != score :
		print("--------\n\tFAIL\nScores don't match\n--------")

	if isinstance(s1, Sequence) and isinstance(s2, Sequence) :
		sequences.append(Profile([s1, s2]))
	elif isinstance(s1, Profile) and isinstance(s2, Profile) :
		sequences.append(Profile(s1.seqs + s2.seqs))
	else :
		assert(False)


if __name__ == "__main__" :

	if args.nw_pp_test :
		print(f'starting score: {get_score(sequences)}')
		while len(sequences) > 1 :
			i = len(sequences) // 2
			while i > 0 :
				nw_pp_test(0, 1)
				print()
				i -= 1
		print('------------------------------------------------')
		print("PASSED")
		print(f'ending score: {get_score(sequences)}')
	
	output = '\n\n'.join([str(seq) for seq in sequences])
	with open(args.outfile, 'w') as file:
		file.write(output)
		
				

	# print("needlman score", needleman_wunsch_ss(sequences[0], sequences[1]))
	# print("score", get_score(sequences[0:2]))
	# print(sequences[0])
	# print(sequences[1])

	# print("Before alignment")
	# print("score", get_score(sequences))
	# print(sequences[0])
	# print(sequences[1])
	# print(sequences[2])
	# print()

	# prof1 = Profile(sequences[:2])
	# prof2 = Profile(sequences[2:])

	# print("After alignment")
	# print("needlman score", needleman_wunsch_pp(prof1, prof2))
	# print("score", get_score(sequences))
	# print(sequences[0])
	# print(sequences[1])
	# print(sequences[2])

	# neighbor_tree = join_neighbors(sequences)
	# for n in neighbor_tree :
	# 	print(n)





