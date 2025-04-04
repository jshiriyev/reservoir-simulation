import timeit

SETUP_CODE = """

import numpy as np

from vectorpy import starsplit
from vectorpy import starsplit_numpy
from vectorpy import starsplit_npvect

from vectorcy import starsplit as starsplitcy

N = 1000000

np.random.seed(2021)
part1 = np.random.randint(1,5,N)
part2 = np.round_(np.random.rand(N)*10000,1)

fstring = ("{}*{}")

def star_merge_(repeat,value):

	if repeat==1:
		return f"{value}"
	else:
		return f"{repeat}*{value}"

vmerge = np.vectorize(star_merge_)

result_array = vmerge(part1,part2)

result_list = result_array.tolist()

result_string = " ".join(result_list)

"""

RUN_CODE_1 = """

pylist = starsplit(result_list)

print(pylist[:10])

print(len(pylist))

"""

RUN_CODE_2 = """

nparray = starsplit_numpy(result_array)

print(nparray[:10])

print(nparray.size)

"""

# RUN_CODE_3 = """

# ints,flts = starsplit_npvect(result_array)

# nparray = np.repeat(flts,ints)

# print(nparray[:10])

# """

RUN_CODE_4 = """

pylist = starsplitcy(result_list)

print(pylist[:10])

print(len(pylist))

"""

py1 = timeit.timeit(RUN_CODE_1,setup=SETUP_CODE,number=1)
py2 = timeit.timeit(RUN_CODE_2,setup=SETUP_CODE,number=1)
# py3 = timeit.timeit(RUN_CODE_3,setup=SETUP_CODE,number=1)
py4 = timeit.timeit(RUN_CODE_4,setup=SETUP_CODE,number=1)

print(f"Python List runs in {py1}")
print(f"Numpy Array runs in {py2}")
# print(f"Vectorized Numpy Array runs in {py3}")
print(f"Cythonized Python List runs in {py4}")