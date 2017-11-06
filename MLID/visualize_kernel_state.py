import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Read in data
with open("kernel_state.data") as f:
    content = f.readlines()
content = [[float(n) for n in line.split()] for line in content]

with open("kernel_u.data") as f:
    content2 = f.readlines()[1:]
content2 = [[float(n) for n in line.split()] for line in content2]

fig = plt.figure()
ax = Axes3D(fig)
fig2 = plt.figure()
ax2 = Axes3D(fig2)

threshold = 0.1

for i in range(0, len(content)):
    x = float(content[i][0])
    y = float(content[i][1])
    z = float(content[i][2])

    a = float(content2[i][0])
    b = float(content2[i][1])
    t = float(content2[i][2])

    if float(content[i][-1]) == 0.0:
        ax.scatter(x, y, z, c='r')
        ax2.scatter(a, b, t, c='r')
    elif float(content[i][-1]) < threshold:
        ax.scatter(x, y, z, c='b')
        ax2.scatter(a, b, t, c='b')

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('theta')

ax2.set_xlabel('a')
ax2.set_ylabel('b')
ax2.set_zlabel('t')

#ax2.set_xlim(-4.8, -4.7)
#ax2.set_ylim(-2.0, 1.0)
#ax2.set_zlim(-0.3, 0.1)

#ax.set_aspect('equal')

plt.show()
print ""