from multiprocessing import Queue, Process, cpu_count
from scoring import get_score, blast_config
from sequence import Sequence, Profile
from guide_tree import Tree, Node
from needleman_wunsch import needleman_wunsch_ss, needleman_wunsch_pp
from timer import Timer
import sys

"""
	Worker  ----- to_main -----> Main
      ^                           /
       \------- to_worker -------

	messages sent on to_worker channel
		(seq, seq, id) - the seqs (can be Sequence or Profile) of will be aligned
	
	messages sent on to_main channel
		(Profile, id, perf_dict) - the aligned profile and the id are returned

"""


def worker(in_channel, out_channel, config=blast_config) :
	try :
		while True :
			param = in_channel.get()
			if param == None : # end process on None input
				return
			assert isinstance(param, tuple), "worker input must be tuple"
			assert len(param) == 3, "worker input must be tuple of len 2"
			s1, s2, up_id = param
			assert isinstance(s1, Sequence) or isinstance(s1, Profile), f"worker input s1 must be tuple of 2 seqs: {repr(s1)}"
			assert isinstance(s2, Sequence) or isinstance(s2, Profile), f"worker input s2 must be tuple of 2 seqs: {repr(s1)}"
			
			start_score = get_score(s1, s2, scoring_config=config)
			names = (repr(s1), repr(s2))

			alignment_function = None

			if isinstance(s1, Sequence) and isinstance(s2, Sequence) :
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
				assert(False), "invalid pair of sequences"

			nw_score = alignment_function(s1, s2, scoring_config=config)
			score = get_score(s1, s2, scoring_config=config)
			time_elapsed = Timer.log(alignment_function.__name__, log=False)
			prof = None

			if isinstance(s1, Sequence) and isinstance(s2, Sequence) :
				prof = Profile([s1, s2])
			elif isinstance(s1, Profile) and isinstance(s2, Profile) :
				prof = Profile(s1.seqs + s2.seqs)
			else :
				assert(False)

			perf_dict = {'start_score': start_score, 'names': names, 'nw_score': nw_score, 'score': score, 'time_elapsed': time_elapsed}
			output = (prof, up_id, perf_dict)
			out_channel.put(output)
	except KeyboardInterrupt :
		return


def print_progress_bar (iteration, total):
	"""
	Call in a loop to create terminal progress bar
	@params:
		iteration   - Required  : current iteration (Int)
		total       - Required  : total iterations (Int)
	"""
	# if iteration != 0 :
	# 	sys.stdout.write("\r")
	fraction = f"{iteration} / {total} completed"
	filled_length = int(100 * iteration // total)
	bar = '█' * filled_length + '-' * (100 - filled_length)
	print(f'alignments |{bar}| {fraction}', end = '\r')
	# Print New Line on Complete
	if iteration == total: 
		print()

def align_tree(tree, threads=None) :
	assert isinstance(tree, Tree)
	if threads == None :
		threads = cpu_count()
	out_channel = Queue()
	in_channel = Queue()
	workers = []
	for _ in range(threads) :
		w = Process(target=worker, args=(out_channel, in_channel))
		w.start()
		workers.append(w)

	for node1, node2 in tree.start_pairs :
		assert node1.up is node2.up and node1.up is not None
		out_channel.put((node1.seq, node2.seq, node1.up.id))
	
	# progress bar
	num_todo = tree.alignments_to_completion()
	completed = 0
	print_progress_bar(completed, num_todo)

	cpu_time = 0
	while tree.root.seq is None :
		prof, parent_id, perf_dict = in_channel.get()
		# print(f'finished alignment between :\n\t{perf_dict["names"][0]}\n\t\tand\n\t{perf_dict["names"][1]}')
		# print(f'\t\tscore: {perf_dict["start_score"]:,} -> {perf_dict["score"]:,} (nw score: {perf_dict["nw_score"]:,})')
		# print(f'\t\ttime elapsed: {perf_dict["time_elapsed"]:5.2f}s')
		completed += 1
		print_progress_bar(completed, num_todo)

		cpu_time += perf_dict["time_elapsed"]
		parent = tree.nodes[parent_id]
		parent.seq = prof
		if parent.up is not None :
			grandparent = parent.up
			uncle = grandparent.left if grandparent.right.id == parent_id else grandparent.right
			if uncle.seq is not None :
				out_channel.put((uncle.seq, parent.seq, grandparent.id))

	for _ in workers : # tell all workers to die
		out_channel.put(None)
	for w in workers : # wait for workers to die
		w.join()
	print("------------")
	print("align_tree finished")
	print(f"\tcpu_time: {cpu_time:5.2f}")
	print(f"\tfinal score: {get_score(tree.root.seq.seqs):,}")
	return tree.root.seq.seqs
	

