import argparse
from distance import pairwise_distance
from sequence import read_fasta, Sequence
from scoring import get_score
from distance import pairwise_distance, join_neighbors
from needleman_wunsch import needleman_wunsch

parser = argparse.ArgumentParser(description='Beep boop align sequence')
parser.add_argument('infile', help='input file in .fasta format')
parser.add_argument('outfile', help='output file name')
parser.add_argument('-genetic', '-g', action='store_true', 
	help='align with genetic algorithm')
args = parser.parse_args()
# print('Args=', args)

sequences = read_fasta(args.infile)

if __name__ == "__main__" :
	print("Before alignment")
	print("score", get_score(sequences[0:2]))
	print(sequences[0])
	print(sequences[1])
	print()
	print("After alignment")
	print("needlman score", needleman_wunsch(sequences[0], sequences[1]))
	print("score", get_score(sequences[0:2]))
	print(sequences[0])
	print(sequences[1])

	# neighbor_tree = join_neighbors(sequences)
	# for n in neighbor_tree :
	# 	print(n)





