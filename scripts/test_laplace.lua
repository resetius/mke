#!/usr/bin/lua

function read_errors(f)
	local e1 = 1.0
	local e2 = 1.0
	local e3 = 1.0

	while true do
		local line = f.read("*line")
		if not line then break end
	end

	if (math.abs(e3) > 1e-7) then
		print (e3 .. "> 1e-7")
		os.exit(-1)
	end

	return e1, e2, e3
end

function run_test()
	local path = arg[1]
	local exe1 = path .. "/test_laplace"
	local exe2 = path .. "/rectangle"
	print("run in " .. path)
	os.execute(exe2 .. " 0 0 1 1 3 > r3.txt")
	os.execute(exe2 .. " 0 0 1 1 4 > r4.txt")
	os.execute(exe2 .. " 0 0 1 1 5 > r5.txt")
	local e11, e21, e31 = read_errors(io.popen(exe1 .. " r3.txt", r))
	local e12, e22, e32 = read_errors(io.popen(exe1 .. " r4.txt", r))
	if (e11 / e12 < 3.5) then
		print ("e11 / e12 < 3.5")
		os.exit(-1)
	end
	local e13, e23, e33 = read_errors(io.popen(exe1 .. " r5.txt", r))
	if (e12 / e13 < 3.5) then
		print ("e12 / e13 < 3.5")
		os.exit(-1)
	end
end

run_test()
os.exit(0)
