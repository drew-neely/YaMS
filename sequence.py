from enum import Enum, unique
from copy import deepcopy
from typing import Iterable


code = ['a', 'c', 'g', 't', '-', 'n',]

@unique
class NucleicAcid(Enum):
	ADENINE = 0
	CYTOSINE = 1
	GUANINE = 2
	THYMINE = 3
	GAP = 4
	ANY = 5

	def __str__(self) :
		return code[self.value]



class Sequence :
	def __init__(self, name, seq) :
		self.name = name
		self.seq = bytearray([code.index(c) for c in seq.lower()])

	def __str__(self) :
		seq_str = ''.join([code[a] for a in self.seq])
		return f'>{self.name}\n{seq_str}\n\n'

	def __getitem__(self, index) :
		if isinstance(index, slice) :
			return [NucleicAcid(a) for a in self.seq[index]]
		else :
			return NucleicAcid(self.seq[index])

	def __setitem__(self, index, acid) :
		if isinstance(index, slice) :
			if isinstance(acid, Iterable) :
				self.seq[index] = [a.value for a in acid]
			else :
				print(index, acid.value)
				self.seq[index] = [acid.value]
		else :
			self.seq[index] = acid.value

	def copy(self) :
		return deepcopy(self)
		

# filename => name of fasta file containing one or more sequences
# returns => list of sequence objects with each object corresponding to one sequence from the file
def read_fasta(filename) :
	with open(filename, 'r') as file :
		contents = file.read()
	contents = contents.strip()
	entries = [entry for entry in contents.split('>') if len(entry) > 0]
	sequences = []
	for entry in entries :
		name = entry[:entry.find('\n')]
		seq = ''.join(entry[entry.find('\n')+1:].split())
		sequences.append(Sequence(name, seq))
	return sequences