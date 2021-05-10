YaMS - **Y**et **A**nother **M**ultiple **S**equence aligner
============================================================

Usage
-----

YaMS is a generic DNA sequence aligner that takes multiple sequences in
a single FASTA as input and outputs an alignment of these sequences.

    python yams.py <path to input fasta> <path to output fasta>

YaMS should be run with PyPy3.7. YaMS will work being run by CPython,
but will run about 10 times slower.

PyPy3.7 can be obtained from the PyPy download site:
<a href="https://www.pypy.org/download.html" class="uri">https://www.pypy.org/download.html</a>

positional arguments:

| Argument |         Description         |
|:--------:|:---------------------------:|
|  infile  | input file in .fasta format |
|  outfile |       output file name      |

optional arguments:

<table>
<colgroup>
<col style="width: 32%" />
<col style="width: 67%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: center;">Flag(s)</th>
<th style="text-align: center;">Description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: center;">-h, –help</td>
<td style="text-align: center;">show help message and exit</td>
</tr>
<tr class="even">
<td style="text-align: center;">-no-guide-tree, -ngt</td>
<td style="text-align: center;">will use guide tree based on order rather than distance if specified</td>
</tr>
<tr class="odd">
<td style="text-align: center;">-kmer-distance k, -kdist k</td>
<td style="text-align: center;">k value to use in kmer distance calculation</td>
</tr>
<tr class="even">
<td style="text-align: center;">-score-distance, -sdist</td>
<td style="text-align: center;">use needleman wunsch alignment to calculate distance</td>
</tr>
<tr class="odd">
<td style="text-align: center;">-iterate-to-convergence, -converge</td>
<td style="text-align: center;">iterate until convergence is reached</td>
</tr>
<tr class="even">
<td style="text-align: center;">-iterative-depth d, -depth d</td>
<td style="text-align: center;">number of iterative steps to perform</td>
</tr>
<tr class="odd">
<td style="text-align: center;">-no-iterative, -noit</td>
<td style="text-align: center;">don’t perform the iterative step at all</td>
</tr>
<tr class="even">
<td style="text-align: center;">-threads t, -t t</td>
<td style="text-align: center;">Number of threads to use - default is the number of cpus avaliable</td>
</tr>
</tbody>
</table>

Introduction
============

Sequence alignment is the first step of many analysis that involve DNA,
RNA, or Amino Acids. The goal of sequence alignment is to identify the
locations in which various mutation events occurred since the most
recent shared ancestor. These mutations are insertions, deletions, and
point mutations. Correctly identifying insertions and deletions is done
by adding gaps into the sequences in attempts to shift the rest of the
elements to align better. There are two parts to making a high quality
MSA:

-   Choosing a good scoring function by which to evaluate a possible
    alignment
-   Optimizing the solution within the search space defined by the
    scoring function

An MSA of a set of sequences contains information about the way those
sequences evolved in the subject organisms. The scoring function must be
chosen such that the quality of an MSA represents the likelihood that
the mutations actually occurred the way the MSA describes. As such, the
scoring function must be heavily based in the biology and chemistry of
DNA. Very advanced scoring functions have been created that account for
the structural properties of DNA in fine detail. **<sup>\[1\]</sup> I
chose not to focus in depth on this aspect of MSA. YaMS uses a
sums-of-pairs substitution cost and an affine gap penalty. While a
sums-of-pairs model is not entirely accurate in modeling the biology of
DNA mutations, it is much simpler than other more accurate models,
allowing for increased performance and faster runtimes.

The bulk of the time and effort put into YaMS was spent on creating an
alignment that attempts to maximize the score given by the scoring
function. YaMS follows the MSA model used in MUSCLE which uses several
steps to construct and refine an alignment. The first step is a
progressive alignment in which sequences are incorporated into the MSA
in order of similarity. The second is an iterative progressive alignment
that attempts to optimize the order in which sequences are included.
MUSCLE has a third refinement step aiming to resolve some of the
inaccuracies of progressive alignment, but this is not included in YaMS.
**<sup>\[2\]</sup>

Implementation
--------------

The input to a YaMS is a single FASTA file that contains two or more DNA
sequences. The sequences are identified by their names in the FASTA
file, so the names must be unique. Once the FASTA file is read into
memory and parsed, the alignment process begins. Both stage output a
valid alignment, but the second attempts to refine and improve the
alignment output by the first.

![Flow chart describing the stages and intermediate steps of a YaMS
alignment](Resources/Flow_chart.png)

### Pairwise Alignment

The fundamental process that drives a progressive alignment is a
pairwise alignment (an alignment of two sequences). The pairwise
alignment algorithm used in YaMS is an implementation of Needleman
Wunsch modified in 3 ways. YaMS contains two implementations of modified
Needleman-Wunsch algorithms. The first is Needleman-Wunsch with affine
gap penalty which is used for sequence-sequence alignments. The second
is a profile-profile Needleman-Wunsch with affine gap penalty used for
profile-profile and profile-sequence alignments. The later is also
capable of doing sequence-sequence alignments, but the former has better
performance so is used when possible.

#### Needleman-Wunsch

The Needleman-Wunsch algorithm uses a scoring function that assigns a
value to each possible pair of Nucleic Acids and a constant gap penalty.

![The basic Needleman-Wunsch algorithm: g = gap penalty, s(x,y) = score
of aligning nucleic acids x and y.
**<sup>\[4\]</sup>](Resources/nw_level1.png)

Given two sequences x and y, where *x*<sub>*i*</sub>, *y*<sub>*i*</sub>
are the *i*<sup>*t**h*</sup> nucleic acids in each sequence, the
algorithm works by computing a score matrix, V, and a direction matrix,
D. *V*<sub>*i*, *j*</sub> is the best possible score of the alignment
between the first i elements of x and the first j elements of y. Thus,
the *i*<sup>*t**h*</sup> index in the first row and column represent
aligning the *i*<sup>*t**h*</sup> element of one sequence to the first
element of the other by inserting gaps at the beginning of one of the
sequences. For this reason the first row and column are initialized to
*i* × *g*. For every other index, *V*<sub>*i*, *j*</sub> is a function
of the index to the left, the index above and the index to the left
diagonal. These three indices contain the best score for the three
possible alignments for the previous nucleic acids in the sequence. The
score of aligning *x*<sub>*i*</sub> and *y*<sub>*j*</sub> may be the
score of aligning:

1.  *x*<sub>*i*</sub> to a gap after *y*<sub>*j*</sub>
2.  *y*<sub>*j*</sub> to a gap after *x*<sub>*i*</sub>
3.  *x*<sub>*i*</sub> to *y*<sub>*j*</sub>

The maximum of these options is chosen, and the choice is stored in the
*D*<sub>*i*, *j*</sub> as an arrow referencing the index from which the
maximum value was calculated. If there is more than one maximum of these
three values, then there is more than one optimal alignment. When this
occurs, YaMS just records one of the options in *D*<sub>*i*, *j*</sub>
and ignores the other. The way this is chosen is undefined.

*V*<sub>*s**i**z**e*(*x*), *s**i**z**e*(*y*)</sub> contains the score of
aligning the entire sequence. Using the arrows stored in the direction
matrix, the optimal way to insert gaps into the sequence is encoded in
the path from *D*<sub>*s**i**z**e*(*x*), *s**i**z**e*(*y*)</sub> to
*D*<sub>0, 0</sub>. This path is a list of arrows. The path is iterated
over from end to beginning. For a LEFT or RIGHT arrow, a gap is inserted
in the corresponding location in the correct sequence. For a DIAG arrow,
the corresponding locations in the sequences are matched or mismatched.
The result is the alignment between the two sequences that yields the
global optimal score.

#### Needleman-Wunsch with Affine Gap Penalty

Since insertion/deletion events can involve multiple sites, more
accurate sequences can be obtained by distinguishing between an opening
gap and a gap extension. YaMS by default gives the first gap in a
sequence of 1 or more gaps a penalty of -9. Each other gap in the
sequence receives a score of -1. Thus, given the opening gap penalty, d,
and the gap extension penalty e, the total penalty for a sequence of n
gaps is *d* + (*n* − 1) × *e*. The standard Needleman-Wunsch algorithm
does not allow for affine gap penalty in the scoring function, so YaMS
uses an adapted version of the algorithm.**<sup>\[1, 6\]</sup>

The penalty for inserting a gap at *x*<sub>*i*</sub> or
*y*<sub>*j*</sub> (ie. letting *D*<sub>*i*, *j*</sub> equal LEFT or
RIGHT) is dependent on whether the optimal alignment of the sequences up
to that point ends with a gap or not. This makes for 5 possible
alignment choices rather than the 3 in the standard Needleman-Wunsch
algorithm.

1.  *x*<sub>*i*</sub> to a gap after *y*<sub>*j*</sub> given
    *x*<sub>*i* − 1</sub> aligns to *y*<sub>*j*</sub>
2.  *x*<sub>*i*</sub> to a gap after *y*<sub>*j*</sub> given
    *x*<sub>*i* − 1</sub> aligns to a gap
3.  *y*<sub>*j*</sub> to a gap after *x*<sub>*i*</sub> given
    *y*<sub>*j* − 1</sub> aligns to *x*<sub>*i*</sub>
4.  *y*<sub>*j*</sub> to a gap after *x*<sub>*i*</sub> given
    *y*<sub>*j* − 1</sub> aligns to a gap
5.  *x*<sub>*i*</sub> to *y*<sub>*j*</sub>

Choices 1 and 3 involve inserting a gap after a match or mismatch, so
the opening gap penalty should be used. Choices 2 and 4 involve
inserting a gap after a gap, so the gap continuation penalty should be
used. Finally, choice 5 should use the score
*s*(*x*<sub>*i*</sub>, *y*<sub>*j*</sub>).**<sup>\[6\]</sup>

![Needleman-Wunsch algorithm with affine gaps.
**<sup>\[6\]</sup>](Resources/nw_level2.png)

This is implemented through the use of 3 additional matrices F, G and H.
*F*<sub>*i*, *j*</sub> is the score of choice 5. *G*<sub>*i*, *j*</sub>
is the maximum of choices 1 and 2. *H*<sub>*i*, *j*</sub> is the maximum
of choice 3 and 4. As such, matrix F contains the score yielded when a
match is allowed, matrix G contains the scores when a gap is inserted in
y, and matrix H contains the scores when a gap is inserted in x. The
best of the five scores is stored in matrix V, and the choice made is
stored in matrix D represented by an arrow. The optimal alignment is
constructed from matrix D in the same way.

#### Profile-Profile Needleman-Wunsch with Affine Gap Penalty

To perform an MSA, it must be possible to perform an alignment between
more than two sequences. The way YaMS does this is by doing profile
alignments. A profile is a set of sequences that have already been
aligned to each other. It is itself an MSA. The input to the
profile-profile aligner is two profiles *x* and *y*, where
*x*<sub>*i*</sub> and *y*<sub>*i*</sub> are vectors containing the
counts of each amino acid at the *i*<sup>*t**h*</sup> position of the
profile’s sequences. These counts do not include gaps that may exist in
some sequences at the *i*<sup>*t**h*</sup> position.
*x*<sub>*i*</sub><sup>*o*</sup> and *y*<sub>*i*</sub><sup>*o*</sup> are
the number of opening gaps at the *i*<sup>*t**h*</sup> position of the
profiles while *x*<sub>*i*</sub><sup>*e*</sup> and
*y*<sub>*i*</sub><sup>*e*</sup> are the number of gap extension at the
*i*<sup>*t**h*</sup> position of the profiles. Additionally
*x*<sub>*i*</sub><sup>*i**o*</sup> and
*y*<sub>*i*</sub><sup>*i**o*</sup> are defined as the number of opening
gaps that would exist in the profiles at the *i*<sup>*t**h*</sup>
position if a gap where to be inserted at the *i*<sup>*t**h*</sup>
position in every sequence in the profiles. Likewise,
*x*<sub>*i*</sub><sup>*i**e*</sup> and
*y*<sub>*i*</sub><sup>*i**e*</sup> are the number of gap extensions that
would exist at the *i*<sup>*t**h*</sup> position in each profile under
the same condition. Finally, *x*<sup>*s*</sup> and *y*<sup>*s*</sup> are
the number of sequences included in each profile.

When a pofile-profile alignment inserts a gap in one of the profiles at
the *i*<sup>*t**h*</sup> position, it inserts a gap in all of that
profile’s sequences at the *i*<sup>*t**h*</sup> position. Unlike the
sequence-sequence alignment algorithms described above, profile-profile
alignment is not globally optimal. However, this restriction is applied
in order to improve the computability of alignments. Finding an optimal
alignment between n sequences of length m requires
*O*(*m*<sup>*n*</sup>) time and space.**<sup>\[3\]</sup> For all
reasonably sized inputs, this is not acceptable. A single alignment
between two profiles containing any number of sequences of length m
requires only *O*(*m*<sup>2</sup>). Thus, profile-profile alignments are
a necessary tradeoff between accuracy and performance.

![Profile-Profile Needleman-Wunsch algorithm with affine
gaps](Resources/nw_level3.png)

The two differences between a profile-profile and sequence-sequence
alignment in YaMS are the way the match scores are calculated and the
way gap types are identified. Since each profile may contain multiple
nucleic acids at each index, the score of aligning two indices is
calculated by taking the sum of scores of all pairs of nucleic acids
that are aligned at the same index. The function
*m**m*<sub>*s*</sub>*c**o**r**e* takes the counts of each amino acid and
computes this score combinatorially. Much like the sequence-sequence
version of Needleman-Wunsch, the profile-profile version builds 4
matrices that contain the maximum score of aligning every pair of
nucleic acids under 5 different conditions. However, since the profile
inputs are MSAs themselves, there may already be gaps in some sequences.
Thus, inserting a gap at a given position in a profile may be inserting
an open gap in some the sequences and inserting a gap extension in
others. Additionally, aligning one index of a profile to an index in
another may still result in gaps at those indices. This algorithm takes
this into account by counting and weighing existing gaps as well as
determining the types of gaps that would be inserted at each place in
each sequence.

I could not find in depth descriptions of the profile-profile alignment
algorithms used by existing sequence aligners, so it is possible that
extending Needleman-Wunsch in this way is unique to YaMS.

#### Pairwise Alignment Optimizations

The vast majority of the time spent performing a YaMS MSA is spent
performing pairwise sequence alignments, so most of the optimization
work I performed was in this area. A single pairwise alignment consists
of an initialization stage where all of the matrices are created and
initialized, followed by a loop iterating over the indices of the
matrices which consumes most of the time. The simplest optimization that
gave a speedup of about 2 was to calculate and store all reusable values
(ie. amino acid counts at each position) before the loop, rather than
recalculating them for every position. The most complex and also most
effective optimizations implemented in YaMS are cache locality
optimizations.

The sequence-sequence and profile-profile alignments algorithms create
either 2 or 5 matrices. The size of each of these matrices is the
product of the length of the input sequences. Since input sequences may
frequently be several thousand acids long, these matrices can become
very large, and require hundreds of megabytes of memory. This is far
larger than the cache available to most average desktop computers.
Originally, the matrices were implemented as lists of lists. In python,
this means that the rows of each matrix were scattered throughout
memory. Additionally, the matrices were originally iterated diagonally
from left to right and bottom to top (filling in the matrices from the
top left to the bottom right). I change the matrices to be implemented
as a single list with manual row-major indexing and changing the
iteration order to be row by row (left to right and top to bottom). This
vastly improved the cache locality by ensuring that the rows of the
matrices were close together and clustering sequentially accessed memory
locations. These changes alone yielded a speedup of 20.

Additionally, I chose to use PyPy instead of CPython to execute the
program. CPython is the standard python interpreter, and is much more
widespread and commonly used. However, without having to modify the code
at all PyPy yields a speedup of about 10. It does this by performing JiT
compilation on the most time intensive parts of a program. This reduces
the overhead introduced by regular interpreters like CPython, and gives
the benefits of micro-optimazations that are possible when compiling.

Between all of these optimizations, I was able to reduce the time for a
single profile-profile optimization for sequences of length 2400 from
about 25 minutes to 4 seconds. This is a significant and necessary
improvement required if YaMS is to have any practical application.

### Progressive Alignment

The first stage of the alignment process is a progressive alignment. A
pairwise distance matrix is computed using kmers distance approximation.
Then a guide tree is built from this matrix using nearest neighbor
clustering.**<sup>\[5\]</sup> The guide tree is a binary tree used to
determine the order in which sequences will be aligned. Each leaf node
of the guide tree represents a single sequence while each internal node
represents a profile. The profile represented by each internal node is
an MSA of all the sequences at the leafs of the node’s respective
subtree. The final MSA is represented by the profile at the root of the
tree. The iteration of the guide tree is parallelized as follows.

-   The initial pairs of sequences to align are computed by iterating
    the guide tree and recording all sibling leaf nodes.
-   The main thread sends each initial pair on a pipe to to a pool of
    worker threads. It then waits for responses.
-   A worker picks up a pair of sequences to align, performs the
    pairwise alignment, combines the pair to create a new profile, then
    returns the new profile to the main thread.
-   The main thread receives the result of the alignment and stores it
    in the corresponding internal node in the tree.
-   If the internal node’s sibling has already been aligned, the main
    thread dispatches the sibling pair to the worker threads.
-   The alignment terminates once the tree’s root has been aligned.

The guide tree level parallelism is limited. For inputs containing a
smaller number of sequences, some threads are likely to be starved early
in the computation. As the part of the guide tree being worked on gets
narrower there become fewer sister pairs. Additional limitations to the
speedup that guide tree level parallelism can provide lies in the fact
that sequence alignments are very memory heavy operations, meaning that
on many computers, the memory access speed will become the bottleneck
for the speed of the program rather than the number of cores executing
alignments. From simple benchmarks on a single machine I estimate that
running on 4 cores gives a speedup of about 2 only while none of the
threads are starved.

### Iterative Progressive Alignment

The second stage of the alignment is an iterative progressive alignment
aiming to resolve some of the errors introduced during the first
progressive alignment. Kmers distance approximation is not a very good
measure of distance due its simplicity, but it is one of the best
available to use on unaligned sequences. However, with an approximate
alignment available, it is possible to more accurately asses the
pairwise distance between sequences. YaMS may use either pdistance or
the inverse of the score depending on which flags are specified.
Pdistance is the proportion of sites that differ. The inverse of the
objective function score may also be used. I suspect from minimal
testing that pdistance performs slightly better, so it is the default
option.

A new distance matrix is created using the pairwise distance function of
choice, and a new guide tree is constructed using nearest neighbor
clustering on the distance matrix.**<sup>\[5\]</sup> This new guide tree
is then used for a new progressive alignment. However, since the new
guide tree frequently has a similar structure to the old one, all nodes
in the new guide tree with subtrees that have identical structure to
that of a node in the old guide tree, need not be recomputed. Comparing
the old guide tree to the new and copying over the alignments that have
already been performed greatly reduces the number of alignments needing
to be run.

This process is repeated until either a specified number of iterations
have been completed or the number of subtrees that changed from one
iteration to the last stops decreasing.**<sup>\[2\]</sup>

Discussion
----------

YaMS is a fully functional multiple sequence aligner capable of
producing alignments that look pretty okay. From visual examination of
the results, a human can pick out regions that should have some simple
and apparently obvious modifications. I believe this to be a result of
the flaws of pairwise alignment. Every MSA created via pairwise
alignment contains at least one pair of sequences that were originally
aligned to each other without regard for any of the other sequences in
the MSA. Since there is no mechanism by which gaps may be relocated or
removed in progressive alignment, once a gap is inserted in these early
alignments, it will be a part a part of the result. MUSCLE uses a
refinement stage that contains a mechanism through which gaps may be
removed.**<sup>\[2\]</sup> YaMS could possibly benefit extensively from
the addition of a refinement stage.

As far as performance is concerned, YaMS runs at a similar speed to
MUSCLE. This is a good indicator of the quality of the YaMS
implementation. MUSCLE is written in C++ which is much faster than
python when it comes to scientific algorithms that are very CPU and
memory intensive. I am, however, unsure if converting YaMS to C++ would
yield significant speedup. PyPy may already be providing most of the
potential speedup from using compiled code rather than python. The
biggest area of improvement I imagine a C++ implementation giving is in
compressing the size of the arrays used in the pairwise alignment
algorithms and thus improving cache locality. Python represents
everything, including integers, as objects, making storing an integer
require significantly more space than would be required in a language
such as C++. Reducing the size of the arrays increases the proportion of
the array that can be stored in a cache at any given time. However, I am
unfamiliar with what PyPy may already be doing in this area.

References
----------

\[1\] Chatzou, M., Magus, C., & Chang, J. (2015). Multiple sequence
alignment modeling: Methods and applications. Oxford Academic. doi:
<a href="https://doi.org/10.1093/bib/bbv099" class="uri">https://doi.org/10.1093/bib/bbv099</a>

\[2\] Edgar, R. C. (2004). MUSCLE: A multiple sequence alignment method
with reduced time and space complexity. BMC Bioinformatics, 5(1), 113.
doi: 10.1186/1471-2105-5-113

\[3\] Generalized Dynamic Programming for Multiple Sequence Alignment.
(n.d.). Retrieved May 10, 2021, from
<a href="http://readiab.org/book/0.1.2/2/3" class="uri">http://readiab.org/book/0.1.2/2/3</a>

\[4\] Multiple sequence alignment. (2021, May 07). Retrieved May 10,
2021, from
<a href="https://en.wikipedia.org/wiki/Multiple_sequence_alignment" class="uri">https://en.wikipedia.org/wiki/Multiple_sequence_alignment</a>

\[5\] Neighbor joining. (2021, April 23). Retrieved May 10, 2021, from
<a href="https://en.wikipedia.org/wiki/Neighbor_joining" class="uri">https://en.wikipedia.org/wiki/Neighbor_joining</a>

\[6\] Sequence Alignment Lecture. (n.d.). Retrieved May 10, 2021, from
<a href="https://web.stanford.edu/class/cs262/presentations/lecture2.pdf" class="uri">https://web.stanford.edu/class/cs262/presentations/lecture2.pdf</a>
