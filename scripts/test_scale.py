#!/usr/bin/python

import sys
import re
from popen2 import popen2

# example
# python test_scale.py "./bin/test_solver --task mult_dense -d --iters 100 --dim %d -t %d" "100,5000"

cmd = sys.argv[1]
dims = sys.argv[2]

dims = map(int, dims.split(","))
dims = range(*dims)

r = re.compile(r".* ([0-9]+\.[0-9]+).*")

tries = 3

def run_single(c):
	print >> sys.stderr, "run '%s'" % c
	cout, cin = popen2(c)
	cin.close()
	a = None
	for line in cout:
		m = r.match(line)
		#print "line === ", line
		if m:
			a = float(m.groups()[0])
	
	return a

def run_benchmark(t):
	ret = []
	for dim in dims:
		c = cmd % (dim, t)
		r = []
		for tr in xrange(tries):
			r.append(run_single(c))
		ret.append(min(r))

	return ret

def build_scale(r1, r):
	ret = []
	for t1, t2 in zip(r1, r):
		print "%f/%f" % (t1, t2)
		ret.append(t1/t2)
	return ret

def output_scale(fname, scale):
	f = open(fname, "w")
	for dim, s in zip(dims, scale):
		print >> f, "%d %f" % (dim, s)

res1 = run_benchmark(1)
res2 = run_benchmark(2)
res3 = run_benchmark(3)
res4 = run_benchmark(4)

scale2 = build_scale(res1, res2)
scale3 = build_scale(res1, res3)
scale4 = build_scale(res1, res4)

output_scale("scale2.txt", scale2)
output_scale("scale3.txt", scale3)
output_scale("scale4.txt", scale4)

output_scale("time1.txt", res1)
output_scale("time2.txt", res2)
output_scale("time3.txt", res3)
output_scale("time4.txt", res4)

