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

	print (f)
	for line in f:lines() do
		e1 = extract(e1, line, "answer nev: U = ([0-9]*.[0-9]*[eE]?[-+][0-9]*)")
		e2 = extract(e2, line, "answer nev: V = ([0-9]*.[0-9]*[eE]?[-+][0-9]*)")
	end

	print ("e1", e1)
	print ("e2", e2)

	return e1, e2
end

function run_test()
	local path = arg[1]
	local exe1 = path .. "/mke_test_system_laplace"
	local exe2 = path .. "/mke_mesh_rectangle"
	print("run in " .. path)
	os.execute(exe2 .. " 0 0 1 1 3 > r3.txt")
	os.execute(exe2 .. " 0 0 1 1 4 > r4.txt")
	os.execute(exe2 .. " 0 0 1 1 5 > r5.txt")
	local e11, e21 = read_errors(io.popen(exe1 .. " r3.txt", "r"))
	local e12, e22 = read_errors(io.popen(exe1 .. " r4.txt", "r"))
	local k1 = e11 / e12
	local k2 = e21 / e22
	print ("e11 / e12 =", k1)
	print ("e21 / e22 =", k2)
	if (k1 < 3 or k2 < 3) then
		print ("e11 / e12 < 3 or e21 / e22 < 3")
		os.exit(-1)
	end
	local e13, e23 = read_errors(io.popen(exe1 .. " r5.txt", "r"))
	local k1 = e12 / e13
	local k2 = e22 / e23
	print ("e12 / e13 =", k1)
	print ("e22 / e23 =", k2)
	if (k1 < 3 or k2 < 3) then
		print ("e12 / e13 < 3 or e22 / e23 < 3")
		os.exit(-1)
	end
end

run_test()
os.exit(0)

