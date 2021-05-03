from time import time
from collections import OrderedDict

start_times = []

named_times = OrderedDict()

class Timer :

	@staticmethod
	def start(id=None) :
		if id is not None :
			if not id in named_times :
				named_times[id] = (0, None)
			assert named_times[id][1] is None, f"Starting timer {id} twice"
			named_times[id] = (named_times[id][0], time())
		else :
			global start_times
			start_times.append(time())
	
	@staticmethod
	def end(id=None) :
		t = time()
		if id is not None and id in named_times :
			total, start = named_times[id]
			assert start is not None, f"Stopping unstarted timer {id}"
			named_times[id] = (total + t - start, None)
		else :
			global start_times
			assert len(start_times) != 0
			t = t - start_times.pop()
			if id is not None :
				print(f'Time elapsed: {t:5.2}s - {id}')
			else :
				print(f'Time elapsed: {t:5.2}s')

	@staticmethod
	def log(id, log=True, reset=True) :
		assert id in named_times, f"timer for {id} not found"
		total, start = named_times[id]
		assert start is None, f"Logging unstopped timer {id}"
		if(log) :
			print(f'Time elapsed for timer {id}: {total:5.2}s')
		if reset :
			del named_times[id]
		return total

	@staticmethod
	def log_all(log=True, reset=True) :
		res = {}
		for id in list(named_times.keys()) :
			total = Timer.log(id, log=log, reset=reset)
			res[id] = total
		return res
