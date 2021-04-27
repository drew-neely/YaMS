from enum import Enum, unique
from copy import deepcopy
from typing import Iterable
from collections import Counter


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

	def __repr__(self) :
		return self.name

	def is_acid(self) :
		return self.value < 4

class Profile :
	def __init__(self, seqs) :
		self.seqs = seqs
		self.seq_len = len(self.seqs[0])
		for seq in seqs :
			if len(seq) != self.seq_len :
				raise ValueError("All sequences in profile must be same length")
	
	def __getitem__(self, index) :
		if isinstance(index, slice) :
			return [self[i] for i in range(index.start, index.stop, index.step)]
		else :
			entries = [0] * 4
			open_gaps = 0
			cont_gaps = 0
			for seq in self.seqs :
				na = seq[index] if index < len(seq) else NucleicAcid.GAP
				if na == NucleicAcid.GAP and index == 0 :
					open_gaps += 1
					continue
				prev_na = seq[index-1] if index != 0 and index - 1< len(seq) else NucleicAcid.GAP
				if na == NucleicAcid.GAP :
					if prev_na == NucleicAcid.GAP :
						cont_gaps += 1
					else :
						open_gaps += 1
				else :
					entries[na.value] += 1
			return (entries, open_gaps, cont_gaps)

	def insert_gap(self, index) :
		for seq in self.seqs :
			seq.insert_gap(index)
	
	def __len__(self) :
		return len(self.seqs)

	def __repr__(self) :
		return f"<PROF {' '.join([seq.name for seq in self.seqs])}>"
	
	def __str__(self) :
		return '\n'.join([str(seq) for seq in self.seqs])
				

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
				self.seq[index] = [acid.value]
		else :
			self.seq[index] = acid.value

	def insert_gap(self, index) :
		if index <= len(self.seq) :
			self.seq.insert(index, NucleicAcid.GAP.value)

	def __len__(self) :
		return len(self.seq)

	def __repr__(self) :
		return f'<SEQ {self.name}>'

	def __str__(self) :
		return f">{self.name}\n{''.join([str(a) for a in self])}"

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