import math
import random
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.special as special
from datetime import datetime


def normalize_angle(angle):
    new_angle = angle
    while new_angle <= -math.pi:
        new_angle += 2.0 * math.pi
    while new_angle > math.pi:
        new_angle -= 2.0 * math.pi
    return new_angle


def draw_trajectory(S):
    data_num = S.shape[0]
    x = np.zeros(data_num)
    y = np.zeros(data_num)
    for i in range(0, data_num):
        x[i] = S[i,0]
        y[i] = S[i,1]
    plt.figure()
    plt.plot(x,y,'r.')
    plt.axis('equal')
    plt.axis([-1.0,1.0,-0.2,1.5])
    plt.draw()


def transfer_function(s, u):
    a = u[0]
    b = u[1]
    t = u[2]

    v = s[3]
    w = s[4]
    theta = s[2]

    s_n = np.array(s, copy=True)

    # The easy ones
    dz = 0.5*b*t*t + w*t
    dv = a*t
    dw = b*t

    # b = 0 and c = 0 case -> just linear acceleration
    if abs(b) < 0.001 and abs(w) < 0.001:
        l = 0.5*a*t*t + v*t
        dx = math.cos(theta) * l
        dy = math.sin(theta) * l
    # b = 0 case -> formula singularity covered
    elif abs(b) < 0.001:
        dx = a * (math.cos(w * t + theta) - math.cos(theta)) / (theta * theta) + ((a * t + v) * math.sin(w * t + theta) - v * math.sin(theta)) / w
        dy = a * (math.sin(w * t + theta) - math.sin(theta)) / (theta * theta) - ((a * t + v) * math.cos(w * t + theta) - v * math.cos(theta)) / w
    else:
        flipped = False
        if b < 0:
            b = -b
            w = -w
            flipped = True
        sb = math.sqrt(b)
        pb15 = math.pow(b, 1.5)
        gamma = math.cos(0.5 * w * w / b - theta)
        sigma = math.sin(0.5 * w * w / b - theta)
        SPI = math.pow(math.pi, 0.5)

        s1, c1 = special.fresnel((w + b*t) / (sb*SPI))
        s0, c0 = special.fresnel(w / (sb*SPI))
        S = s1 - s0
        C = c1 - c0

        dx = SPI * (b * v - a * w) * (sigma * S + gamma * C) / pb15 + (a / b) * (math.sin(0.5 * b * t * t + w * t + theta) - math.sin(theta))
        dy = SPI * (b * v - a * w) * (gamma * S - sigma * C) / pb15 - (a / b) * (math.cos(0.5 * b * t * t + w * t + theta) - math.cos(theta))

        if flipped:
            c2d = math.cos(2 * theta)
            s2d = math.sin(2 * theta)
            dxt = c2d * dx + s2d * dy
            dy = s2d * dx - c2d * dy
            dx = dxt

    s_n[0] += dx
    s_n[1] += dy
    s_n[2] += dz
    s_n[3] += dv
    s_n[4] += dw

    return s_n


def Euler_Integrate(s, u):
    T = 10**(-6)
    a = u[0]
    b = u[1]
    theta = s[2]
    v = s[3]
    w = s[4]

    s_n = np.array(s, copy=True)
    s_n[0] += v * math.cos(theta) * T
    s_n[1] += v * math.sin(theta) * T
    s_n[2] += w * T
    s_n[3] += a * T
    s_n[4] += b * T

    return s_n


def Jacobian(s, u):
    a = u[0]
    b = u[1]
    t = u[2]

    x = s[0]
    y = s[1]
    theta = s[2]
    v = s[3]
    w = s[4]

    # The easy ones
    dxdt = (a*t + v) * math.cos(0.5*b*t*t + w*t)
    dydt = (a*t + v) * math.sin(0.5*b*t*t + w*t)
    dzda = 0.0
    dzdb = 0.5*t*t
    dzdt = b*t + w
    dvda = t
    dvdb = 0.0
    dvdt = a
    dwda = 0.0
    dwdb = t
    dwdt = b

    # b = 0 and w = 0 -> linear acceleration
    if abs(b) < 0.00001 and abs(w) < 0.00001:
        dxda = 0.5 * t * t
        dyda = 0.0
        # Singularity at b = 0
        dxdb = 0.0
        dydb = 1.0
    # b = 0
    elif abs(b) < 0.00001:
        dxda = (w*t*math.sin(w*t) + math.cos(w*t) - 1) / (w*w)
        dyda = (math.sin(w*t) - w*t*math.cos(w*t)) / (w*w)
        # Singularity at b = 0
        dxdb = -math.sin(w*t)
        dydb = math.cos(w*t)
    else:
        sgn = 1.0
        if b < 0.0:
            b = -b
            w = -w
            sgn = -1.0

        sb = math.sqrt(b)
        pb15 = pow(b, 1.5)
        pb25 = pow(b, 2.5)
        pb35 = pow(b, 3.5)
        gamma = math.cos(0.5*w*w/b)
        sigma = math.sin(0.5*w*w/b)
        kappa = math.cos(0.5*t*t*b + w*t)
        zeta = math.sin(0.5*t*t*b + w*t)
        SPI = math.pow(math.pi, 0.5)
        s1,c1 = special.fresnel((w + b*t)/(sb*SPI))
        s0,c0 = special.fresnel(w/(sb*SPI))
        c = c1 - c0
        s = s1 - s0
        dsdb = math.sin((b*t+w)*(b*t+w)/(2.0*b)) * ()


# random generate num data points from state 0
def random_generation_0(num, amin, amax, bmin, bmax, tmax, unum, visualize):
    s0 = np.array([0.0, 0.0, 0.0, 0.0, 0.0])
    records = np.zeros((num, 5))
    records_u = np.zeros((num,3))
    # generating data points
    #random.seed(datetime.now())
    for i in range(0, num):
        s = np.array(s0, copy=True)
        angle_sum = 0.0
        w = 0.0
        for u_idx in range(0, unum):
            a = random.uniform(amin, amax)
            b = random.uniform(bmin, bmax)
            t = random.uniform(0.0, tmax)
            u = np.array([a, b, t])
            s = transfer_function(s, u)
            angle_sum += abs(w*t + 0.5 * b * t*t)
            w = b * t
        if angle_sum >= 2.0 * math.pi:
            i -= 1
            continue
        else:
            records[i, :] = s
            records_u[i, :] = u

    # visualization
    if visualize:
        figure = plt.figure()
        ax = Axes3D(figure)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('theta')
        ax.set_xlim(-5.0, 5.0)
        ax.set_ylim(-5.0, 5.0)
        ax.scatter(records[:, 0], records[:, 1], records[:, 2])
        plt.draw()

    return records, records_u


# Main
def main():
    s0 = np.array([0.0, 0.0, 0.0, -30.1769, -0.431406])
    u1 = np.array([0.40257, -1.25378, 0.0146489])
    u2 = np.array([2.17582, 2.4135, 0.123905])
    u3 = np.array([0.162694, -0.942372, 0.825526])
    s1 = transfer_function(s0, u1)
    s2 = transfer_function(s1, u2)
    s3 = transfer_function(s2, u3)
    print ""

if __name__ == "__main__":
    main()

