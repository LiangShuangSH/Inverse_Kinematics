import numpy as np
import matplotlib.pyplot as plt
from FD import transfer_function


def visualize(s0, a1, b1, t1, a2, b2, t2, a3, b3, t3):

    t = 0.0
    delta_t = 0.0001

    u = np.zeros(3)
    traj1 = np.zeros((1, 2))
    traj2 = np.zeros((1, 2))
    traj3 = np.zeros((1, 2))

    while t <= t1:
        u[0] = a1
        u[1] = b1
        u[2] = t
        t += delta_t
        s = transfer_function(s0, u)
        position = np.zeros((1, 2))
        position[0, 0] = s[0]
        position[0, 1] = s[1]
        traj1 = np.append(traj1, position, axis=0)

    s1 = s

    while t <= t1 + t2:
        u[0] = a2
        u[1] = b2
        u[2] = t - t1
        t += delta_t
        s = transfer_function(s1, u)
        position = np.zeros((1, 2))
        position[0, 0] = s[0]
        position[0, 1] = s[1]
        traj2 = np.append(traj2, position, axis=0)

    s2 = s

    while t <= t1 + t2 + t3:
        u[0] = a3
        u[1] = b3
        u[2] = t - t1 - t2
        t += delta_t
        s = transfer_function(s2, u)
        position = np.zeros((1, 2))
        position[0, 0] = s[0]
        position[0, 1] = s[1]
        traj3 = np.append(traj3, position, axis=0)

    return traj1, traj2, traj3

s0 = np.array([0.0, 0.0, 0.0, 0.0, 0.0])

trajo1, trajo2, trajo3 = visualize(s0, 1.0376, 0.00508147, 1.98614, 1.93246, 0.10882, 0.85757, 2.18772, -1.66413, 0.233161)
traj1, traj2, traj3 = visualize(s0, 1.02978,-0.0143787,2.07238,1.90052,0.0882537,0.864398,2.17105,-1.6757,0.24804)
traji1, traji2, traji3 = visualize(s0, 1.18854,-0.00709698,1.97316,1.89831,-0.0451131,0.820281,2.18084,-1.65262,0.189247)

plt.figure()

plt.plot(trajo1[:, 0], trajo1[:, 1], 'm.')
plt.plot(trajo2[:, 0], trajo2[:, 1], 'm.')
plt.plot(trajo3[:, 0], trajo3[:, 1], 'm.')

plt.plot(traj1[:, 0], traj1[:, 1], 'r.')
plt.plot(traj2[:, 0], traj2[:, 1], 'b.')
plt.plot(traj3[:, 0], traj3[:, 1], 'y.')

plt.plot(traji1[:, 0], traji1[:, 1], 'k.')
plt.plot(traji2[:, 0], traji2[:, 1], 'k.')
plt.plot(traji3[:, 0], traji3[:, 1], 'k.')

plt.axis('equal')
plt.show()
print ""
