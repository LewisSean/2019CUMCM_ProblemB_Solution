#!/usr/bin/env python
# _*_ coding:utf-8 _*_
# 文件：pro3.py
import numpy as np
import matplotlib.pyplot as plt

def point_move(N, rx, ry, o, beta):

    for i in np.arange(0, np.size(N, 0)):
        N[i] = [o[0] + rx[i] * np.cos(beta), o[1] + ry[i], o[2] - rx[i] * np.sin(beta)]
    return N


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


def circle_point(x, y, z, num, index):
    if num == 8:
        p_x = x[index:719:90]
        p_y = y[index:719:90]
        p_z = z[index:719:90]
        res = np.array([p_x, p_y, p_z]).T
        return res


def get_angle(index, num, length, depth, seta0):  # length: 1.7/2.0   seta0: 初始偏移角
    if num == 8:
        vector = np.array([np.cos(-seta0 / 180 * np.pi), np.sin(-seta0 / 180 * np.pi)])
        part = - 360 / num
        for ii in np.arange(1, num):
            vector = np.vstack((vector, np.array( [np.cos((-seta0 + part * ii) / 180 * np.pi), np.sin((-seta0 + part * ii) / 180 * np.pi)] )))
        lenZ = np.tan(np.arcsin(depth / length))
        nf = np.array([vector[index, 0], vector[index, 1], lenZ])

        return nf / np.linalg.norm(nf)

def force_point(nf, f, o, p):
    cosZ = nf[2]
    cosY = nf[1]
    cosX = nf[0]
    # 分力
    fx = cosX * f
    fy = cosY * f
    fz = cosZ * f
    # 分力矩
    op = np.array([p[0] - o[0], p[1] - o[1], p[2] - o[2]])
    ttt = np.cross(op, np.array([fx, fy, fz]))
    return [fx, fy, fz], ttt


def init_state(depth):
    o = np.array([0, 0, -depth])
    n = np.array([0, 0, 1])
    FF = [80, 80, 80, 80, 90, 80, 80, 80]
    G = 3.6 * 9.8
    seta0 = 0
    addr = 0
    r = 0.2
    num = 8
    length = 1.7
    if num == 10:
        pass
    if num == 8:
        x, y, z = circle_curving(n, r, o)
        N0 = circle_point(x, y, z, num, addr)
        nList = []
        fList = []
        T = []
        for i in np.arange(0, num):
            nList.append(get_angle(i, 8, length, depth, seta0))
            FT = force_point(nList[i], FF[i], o, N0[i])
            fList.append(FT[0])
            T.append(FT[1])

        F = np.array(fList)
        Ttotal = np.sum(T, 0)
        Ttotal[0] = Ttotal[2] = 0
        Ftotal = np.sum(F, 0)
        Ftotal[2] -= G
        Ftotal[1] = 0
        return Ttotal, Ftotal, N0

def iter_state(Ttotal, Ftotal, v, w, deltaT, o, beta):
    J = 0.5 * 3.6 * 0.04 + 1 / 12 * 3.6 * (0.22 ** 2)
    o += v * deltaT
    beta += deltaT * w
    w += (Ttotal * (deltaT / J))[1]
    a = Ftotal / 3.6
    v += a * deltaT
    return v, w, o, beta

def fun(depth):
    beta = 0
    v = np.array([0., 0., 0.])
    w = 0
    o = np.array([0, 0, -depth])
    FF = [80, 80, 80, 80, 90, 80, 80, 80]
    G = 3.6 * 9.8
    seta0 = 0
    deltaT = 0.0005
    endtime = 0.25
    times = int(endtime / deltaT)
    num = 8
    length = 1.7
    Ttotal, Ftotal, N0 = init_state(depth)
    N = N0
    rx = np.array(N0[:, 0])
    ry = np.array(N0[:, 1])
    set_beta = []
    for i in np.linspace(0, endtime, times + 1):
        v, w, o, beta = iter_state(Ttotal, Ftotal, v, w, deltaT, o, beta)
        N = point_move(N, rx, ry, o, beta)
        nList = []
        fList = []
        T = []
        for j in np.arange(0, num):
            tmp = get_angle(j, num, length, depth, seta0)
            nList.append(tmp)
            FT = force_point(nList[j], FF[j], o, N[j])
            fList.append(FT[0])
            T.append(FT[1])

        F = np.array(fList)
        T = np.array(T)
        Ttotal = np.sum(T, 0)
        Ttotal[0] = Ttotal[2] = 0
        Ftotal = np.sum(F, 0)
        Ftotal[2] -= G
        Ftotal[1] = 0
        set_beta.append(beta)
    set_beta = np.array(set_beta)
    return np.max(set_beta * 180 / np.pi), np.min(set_beta * 180 / np.pi), Ftotal[2]


data = []
for i in np.linspace(0.081, 0.11, 31):
    max_beta, min_beta, FZ = fun(i)
    data.append([i, max_beta - min_beta, FZ])
data = np.array(data)
for i in np.arange(0, np.size(data, 0)):
    print(data[i, 0], data[i, 1], data[i, 2])

