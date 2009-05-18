#!/usr/bin/lua
path = arg[1]
exe1 = path .. "/test_laplace"
exe2 = path .. "/rectangle"
print(path)
os.execute(exe2 .. " 0 0 1 1 3 > r3.txt")
os.execute(exe2 .. " 0 0 1 1 4 > r4.txt")
os.execute(exe2 .. " 0 0 1 1 5 > r5.txt")
f = io.popen(exe1 .. " r3.txt", r)
f = io.popen(exe1 .. " r4.txt", r)
f = io.popen(exe1 .. " r5.txt", r)

