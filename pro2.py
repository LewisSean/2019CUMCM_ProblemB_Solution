#!/usr/bin/env python
# _*_ coding:utf-8 _*_
import numpy as np
import matplotlib.pyplot as plt

J = 0.5 * 3.6 * 0.04 + 1/12 * 3.6 * (0.22 ** 2)
beta = 0
v = np.array([0., 0., 0.])
w = 0
o = np.array([0, 0, -0.11])
n = np.array([0, 0, 1])
F = [80, 80, 80, 80, 90, 80, 80, 80]
G = 3.6 * 9.8
f_init = 68.15
seta0 = 0
deltaT = 0.0005
times = int(0.1 / deltaT)
r = 0.2
num = 8
length = 1.7
depth = 0.11

def circle_curving(n, r, o):
    u = np.array([n[2], 0, -n[0]])
    v = np.cross(u, n)
    u_std = u / np.sqrt(u[0] ** 2 + u[1] ** 2 + u[2] ** 2)
    v_std = v / np.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)
    theta = np.linspace(0, np.pi * 2, 721)
    x = o[0] + r * (u_std[0] * np.cos(theta) + v_std[0] * np.sin(theta))
    y = o[1] + r * (u_std[1] * np.cos(theta) + v_std[1] * np.sin(theta))
    z = o[2] + r * (u_std[2] * np.cos(theta) + v_std[2] * np.sin(theta))
    return x, y, z
    '''
    mpl.rcParams['legend.fontsize'] = 10
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot(x, y, z, label='circle curve')
    ax.legend()
    ax.set_zlabel('Z')  # 坐标轴
    ax.set_ylabel('Y')
    ax.set_xlabel('X')
    plt.show()
    '''


def circle_point(x, y, z, num):
    if num == 8:
        p_x = x[0:719:90]
        p_y = y[0:719:90]
        p_z = z[0:719:90]
        res = np.array([p_x, p_y, p_z]).T
        # print(res)
        # plt.scatter(res[:, 0], res[:, 1])
        # plt.scatter(res[1, 0], res[1, 1])
        # plt.show()
        return res
        '''
        ax = plt.subplot(211, projection='3d')
        ax.scatter(p_x, p_y, p_z, label='circle curve')
        ax.legend()
        ax.set_zlabel('Z')  # 坐标轴
        ax.set_ylabel('Y')
        ax.set_xlabel('X')

        ax = plt.subplot(312, projection='3d')
        ax.plot(x, y, z, label='circle curving')
        ax.legend()
        ax.set_zlabel('Z')  # 坐标轴
        ax.set_ylabel('Y')
        ax.set_xlabel('X')

        ax = plt.subplot(212)
        ax.scatter(p_x, p_y, label='2D from Z axis')
        ax.legend()
        ax.set_ylabel('Y')
        ax.set_xlabel('X')
        ax.grid()
        plt.show()
        '''


def get_angle(index, num, length, depth, seta0):  # length: 1.7/2.0   seta0: 初始偏移角
    if num == 8:
        vector = np.array([np.cos(-seta0 / 180 * np.pi), np.sin(-seta0 / 180 * np.pi)])
        part = 360 / num
        for ii in np.arange(1, num):
            vector = np.vstack((vector, np.array( [np.cos((-seta0 + part * ii) / 180 * np.pi), np.sin((-seta0 + part * ii) / 180 * np.pi)] )))
        lenZ = np.tan(np.arcsin(depth / length))
        nf = np.array([vector[index, 0], vector[index, 1], lenZ])
        return nf / np.linalg.norm(nf)     # 返回方向向量



'''
def pro1_v1(M, m, v2_, e):  # e弹性模量：0.7~0.9
    v1 = ((1 - e) * M + 2 * m) * v2_ / ((1 + e) * M)
    return v1
'''

def force_point(nf, f, o, p):  # arg:n 标准化后的f的方向向量    返回在x,y,z轴的受力和力矩   o:圆心 p:圆上作用点
    cosZ = nf[2]
    cosY = nf[1]
    cosX = nf[0]
    # 分力
    fx = cosX * f
    fy = cosY * f
    fz = cosZ * f
    # 分力矩
    op = np.array([p[0] - o[0], p[1] - o[1], p[2] - o[2]])
    fZ, fY, fX = np.array([f * nf[0], f * nf[1], 0]), np.array([f * nf[0], 0, f * nf[2]]), np.array([0, f * nf[1], f * nf[2]])
    opZ, opY, opX = np.array([op[0], op[1], 0]), np.array([op[0], 0, op[2]]), np.array([0, op[1], op[2]])
    tZ, tY, tX = np.cross(opZ, fZ), np.cross(opY, fY), np.cross(opX, fX)

    return [fx, fy, fz], tX+tY+tZ


def init_state(num):  # deltaT: 1/1000
    if num == 10:
        pass
    if num == 8:
        x, y, z = circle_curving(n, r, o)
        pList = circle_point(x, y, z, num)
        nList = []
        fList = []
        tList = []
        F = [80, 80, 80, 80, 90, 80, 80, 80]
        for i in np.arange(0, num):
            nList.append(get_angle(i, 8, length, depth, seta0))
            FT = force_point(nList[i], F[i], o, pList[i])
            fList.append(FT[0])
            tList.append(FT[1])

        T = tList
        F = np.array(fList)
        N = pList
        Ttotal = np.sum(T, 0)
        Ttotal[0] = Ttotal[2] = 0
        Ftotal = np.sum(F, 0)
        Ftotal[2] -= G
        Ftotal[1] = 0
        return Ttotal, Ftotal, N




def iter_state(Ttotal, Ftotal, v, w, deltaT, o, beta):
    o += v * deltaT
    beta += deltaT * w
    w += (Ttotal * (deltaT / J))[1]
    a = Ftotal / 3.6
    v += a * deltaT
    n = np.array([np.sin(beta), 0, np.cos(beta)])   # 圆心法向量
    return n, v, w, o, beta



Ttotal, Ftotal, N = init_state(num)

for i in np.linspace(0, 0.1, times + 1):
    n, v, w, o, beta = iter_state(Ttotal, Ftotal, v, w, deltaT, o, beta)
    x, y, z = circle_curving(n, r, o)
    pList = circle_point(x, y, z, num)
    nList = []
    fList = []
    tList = []
    F = [80, 80, 80, 80, 90, 80, 80, 80]
    for j in np.arange(0, num):
        nList.append(get_angle(j, num, length, depth, seta0))  # 求各个受力点的法向量
        FT = force_point(nList[j], F[j], o, pList[j])
        fList.append(FT[0])
        tList.append(FT[1])

    T = tList
    F = np.array(fList)
    N = pList
    Ttotal = np.sum(T, 0)
    Ttotal[0] = Ttotal[2] = 0
    Ftotal = np.sum(F, 0)
    Ftotal[2] -= G
    Ftotal[1] = 0
    print(beta)

    if(i == 0.06) :
        plt.scatter(N[:, 0], N[:, 1])
        plt.scatter(N[4, 0], N[4, 1])
        plt.grid()
        plt.show()

sita = np.arccos(n[2])
print(sita * 180 / np.pi)















