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

	print (f)
	for line in f:lines() do
		print (line)
		e1 = extract(e1, line, "L1: invert  err=([0-9]*.[0-9]*[eE]?[-+][0-9]*)")
		e2 = extract(e2, line, "L2: laplace err=([0-9]*.[0-9]*[eE]?[-+][0-9]*)", "%1")
		e3 = extract(e3, line, "L3: laplace err=([0-9]*.[0-9]*[eE]?[-+][0-9]*)", "%1")
	end

	print ("e1", e1)
	print ("e2", e2)
	print ("e3", e3)

	if (math.abs(e3) > 1e-7) then
		print (e3 .. "> 1e-7")
		os.exit(-1)
	end

	return e1, e2, e3
end

function run_test()
	local path = arg[1]
	local exe1 = path .. "/mke_test_laplace test_laplace" -- "/test_laplace"
	local exe2 = path .. "/mke_mesh_rectangle"
	print("run in " .. path)
	os.execute(exe2 .. " 0 0 1 1 3 > r3.txt")
	os.execute(exe2 .. " 0 0 1 1 4 > r4.txt")
	os.execute(exe2 .. " 0 0 1 1 5 > r5.txt")
	local e11, e21, e31 = read_errors(io.popen(exe1 .. " -f r3.txt -d", "r"))
	local e12, e22, e32 = read_errors(io.popen(exe1 .. " -f r4.txt -d", "r"))
	local k = e11 / e12
	print ("e11 / e12 =", k)
	if (k < 3.5) then
		print ("e11 / e12 < 3.5")
		os.exit(-1)
	end
	local e13, e23, e33 = read_errors(io.popen(exe1 .. " -f r5.txt -d", "r"))
	local k = e12 / e13
	print ("e12 / e13 =", k)
	if (k < 3.5) then
		print ("e12 / e13 < 3.5")
		os.exit(-1)
	end
end

run_test()
os.exit(0)

