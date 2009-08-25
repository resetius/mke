#!/usr/bin/lua

function extract(old, str, pattern)
	for w in string.gmatch(str, pattern) do
		return w
	end
	return old
end

function read_errors(f)
	local e1 = nil
	local e2 = nil
	local e3 = nil

	for line in f:lines() do
		e1 = extract(e1, line, "bv: norm ([0-9]*\.[0-9]*[eE]?[-+][0-9]*)")
	end

	print ("e1", e1)

	return e1
end

function run_test(path, tp)
	local exe1 = path .. "/test_barvortex"
	local exe2 = path .. "/sphere --type " .. tp
	print("run in " .. path)
	os.execute(exe2 .. " --coord local --iter 4 > ss4.txt")
	os.execute(exe2 .. " --coord local --iter 5 > ss5.txt")
	local e11 = read_errors(io.popen(exe1 .. " -f ss4.txt --task test --time 0.01 -v 0", "r"))
	local e12 = read_errors(io.popen(exe1 .. " -f ss5.txt --task test --time 0.01 -v 0", "r"))
	local k = e11 / e12
	print ("e11 / e12 =", k)
	if (k < 2) then
		print ("e11 / e12 < 2")
		os.exit(-1)
	end
	if (tonumber(e11) > 1e-4 or tonumber(e12) > 1e-5) then
		print ("e11 > 1e-4 or e12 > 1e-5")
		os.exit(-1)
	end
end

run_test(arg[1], arg[2])
os.exit(0)

