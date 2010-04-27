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
		e1 = extract(e1, line, "jacobian err=([0-9]*\.[0-9]*[eE]?[-+][0-9]*)")
	end

	print ("e1", e1)

	return e1
end

function run_test(path, tp)
	local exe1 = path .. "/mke_test_barvortex"
	local exe2 = path .. "/mke_mesh_sphere --type "  .. tp
	print("run in " .. path)
	os.execute(exe2 .. " --coord local --iter 4 > ss4.txt")
	local e11 = read_errors(io.popen(exe1 .. " -f ss4.txt --task jacobian", "r"))
	if (tonumber(e11) > 1) then
		print ("e11 > 1")
		os.exit(-1)
	end
end

run_test(arg[1], arg[2])
os.exit(0)

