import numpy as np
import matplotlib.pyplot as plt
from FD import transfer_function


a1 = 4.9407
b1 = 4.9824
t1 = 5.0023
a2 = 5.0014
b2 = 4.9977
t2 = 5.0023
a3 = 5.1511
b3 = 5.0130
t3 = 5.0023

t = 0.0
delta_t = 0.1
s0 = np.array([0.0, 0.0, 0.0, 0.0, 0.0])
s = s0
u = np.zeros(3)
traj = np.zeros((1, 2))
while t <= t1:
    u[0] = a1
    u[1] = b1
    u[2] = t + delta_t
    s = transfer_function(s, u)
    position = np.zeros((1, 2))
    position[0, 0] = s[0]
    position[0, 1] = s[1]
    np.append(traj, position, axis=1)

while t <= t2:
    u[0] = a2
    u[1] = b2
    u[2] = t + delta_t
    s = transfer_function(s, u)
    position = np.zeros((1, 2))
    position[0, 0] = s[0]
    position[0, 1] = s[1]
    np.append(traj, position, axis=1)

while t <= t3:
    u[0] = a3
    u[1] = b3
    u[2] = t + delta_t
    s = transfer_function(s, u)
    position = np.zeros((1, 2))
    position[0, 0] = s[0]
    position[0, 1] = s[1]
    np.append(traj, position, axis=1)

plt.figure()
plt.plot(traj[:, 0], traj[:, 1])

plt.show()
print ""
