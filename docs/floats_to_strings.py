import timeit

SETUP_CODE = """

import numpy as np

from vectorpy import float2str

N = 1000000

np.random.seed(2021)
# part1 = np.random.randint(1,5,N)
part2 = np.round_(np.random.rand(N)*10000,3)

result_list = part2.tolist()

print(result_list[:10])

"""

RUN_CODE_1 = """

pylist = float2str(result_list,"{:.0f}")

print(pylist[:10])

"""

py1 = timeit.timeit(RUN_CODE_1,setup=SETUP_CODE,number=1)

print(f"Python List runs in {py1}")