from io import StringIO 

import re

import matplotlib.pyplot as plt

import numpy as np

import dirsetup

from cypy.vectorpy import starsplit

Nx,Ny,Nz = (1,1,2)

Nsurf_nodes = (Nx+1)*(Ny+1)

coord = np.empty((Nsurf_nodes*6,))
zcorn = np.empty((2*Nx*2*Ny*2*Nz,))

comment = None
endsection = "/"

headers = ["COORD","ZCORN"]
running = [coord,zcorn]

test_lines = "\
\n\
\n\
COORD\n\
0  0  0 0  0  10000\n\
10 0  0 10 0  10000\n\
0  10 0 0  10 10000\n\
10 10 0 10 10 10000\n\
/\
\n\
\n\
ZCORN\n\
4*2003 8*2015 4*2032\n\
/\
"

file = StringIO(test_lines)

text = file.read()

sections = text.split(endsection)

for section in sections:

	if section=="":
		continue

	section = section.strip("\n")

	section = re.sub('\n+',' ',section)
	section = re.sub(' +',',',section)

	array = section.split(",")
	
	section_keyword = array.pop(0)

	index = headers.index(section_keyword)

	running[index][:] = np.array(starsplit(array))

coord = running[0].reshape((Nsurf_nodes,6))

zcorn = running[1] #((2*Nx*2*Ny*2*Nz,))

coord1 = coord*1
coord2 = coord*1

coord1[:,2] = zcorn[ :Nsurf_nodes]
coord1[:,5] = zcorn[  Nsurf_nodes:2*Nsurf_nodes]

coord2[:,2] = zcorn[2*Nsurf_nodes:3*Nsurf_nodes]
coord2[:,5] = zcorn[3*Nsurf_nodes:]
            
fig = plt.figure(figsize=(8,8))

ax = fig.add_subplot(projection='3d')

ax.scatter(*coord1[:,:3].T)
ax.scatter(*coord1[:,3:].T)
ax.scatter(*coord2[:,:3].T)
ax.scatter(*coord2[:,3:].T)

ax.set_xlabel("x-axis")
ax.set_ylabel("y-axis")
ax.set_zlabel("z-axis")

plt.show()
