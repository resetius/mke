#!/usr/bin/python

import sys
import re
from popen2 import popen2

cmd = sys.argv[1]

dims = sys.argv[2]
dims = map(int, dims.split(","))
dims = range(*dims)

r = re.compile(".* ([0-9]+\.[0-9]+).*")

def run_benchmark(t):
	ret = []
	for dim in dims:
		c = cmd % (dim, t)
		print >> sys.stderr, "run '%s'" % c
		cout, cin = popen2(c)
		cin.close()
		for line in cout:
			m = r.match(line)
			if m:
				a = float(m.groups()[0])
				#print "->>>>%s: %f" % (line, a)
				ret.append(a)

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
res4 = run_benchmark(4)

scale2 = build_scale(res1, res2)
scale4 = build_scale(res1, res4)

output_scale("scale2.txt", scale2)
output_scale("scale4.txt", scale4)

