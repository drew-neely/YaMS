import argparse
from sequence import read_fasta, Sequence
from scoring import get_score

parser = argparse.ArgumentParser(description='Beep boop align sequence')
parser.add_argument('infile', help='input file in .fasta format')
parser.add_argument('outfile', help='output file name')
parser.add_argument('-genetic', '-g', action='store_true', 
	help='align with genetic algorithm')
args = parser.parse_args()
# print('Args=', args)

sequences = read_fasta(args.infile)

if __name__ == "__main__" :
	print(get_score(sequences))



