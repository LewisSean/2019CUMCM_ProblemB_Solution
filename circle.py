#!/usr/bin/env python
# _*_ coding:utf-8 _*_

import numpy as np
import matplotlib.pyplot as plt
def circle_curving(n, r, o):
    u = np.array([n[2], 0, -n[0]])
    v = np.cross(u, n)
    u_std = u / np.sqrt(u[0]**2 + u[1]**2 + u[2]**2)
    v_std = v / np.sqrt(v[0]**2 + v[1]**2 + v[2]**2)
    theta = np.linspace(0, np.pi * 2, 721)
    print(theta[1] / np.pi * 180)
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

def circle_point(x ,y ,z ,num):
    if num == 8:
        p_x = x[45:719:90]
        p_y = y[45:719:90]
        p_z = z[45:719:90]
        res = np.array([p_x, p_y, p_z]).T

        print(res)
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


def get_angle(p,index,num, length):    # length: 1.7/2.0
    if num == 8:

        vector = np.array([[np.cos(-22.5 / 180 * np.pi), np.sin(-22.5 / 180 * np.pi)],
                [np.cos(-67.5 / 180 * np.pi), np.sin(-67.5 / 180 * np.pi)],
                [np.cos(-112.5 / 180 * np.pi), np.sin(-112.5 / 180 * np.pi)],
                [np.cos(-157.5 / 180 * np.pi), np.sin(-157.5 / 180 * np.pi)]
        ])
        '''
        vector = np.array([[res[0, 0], res[0, 1]],
                           [res[1, 0], res[1, 1]],
                           [res[2, 0], res[2, 1]],
                           [res[3, 0], res[3, 1]]
        ])
        for i in np.arange(0, np.size(vector, 0)):
            vector[i] = vector[i] / np.linalg.norm(vector[i])
        '''
        if index == 0 or index == 1 or index == 2 or index == 3:
            vec_n = np.array([vector[index, 0], vector[index, 1]])
        else: vec_n = np.array([vector[7 - index, 0], -vector[7 - index, 1]])


        for x in np.linspace(p[0] + length * vec_n[0], p[0] + length * 0.5 * vec_n[0], 2000):
            vector = np.array([x-p[0], vec_n[1]/vec_n[0]*(x - p[0]), - p[2]])   # 绳索向量
            if np.fabs(np.linalg.norm(vector) - length) < 0.0002:
                sign = np.array([[1, -1, 1], [-1, -1, 1], [-1, 1, 1], [1, 1, 1]])
                vector = np.fabs(vector)
                vector = vector * sign[int(index / 2)] / np.linalg.norm(vector)
                return [x, vec_n[1]/vec_n[0] * (x - p[0]) + p[1], 0], vector  # 绳索末端位置 + 方向向量



def pro1_v1(M, m, v2_, e):    #e弹性模量：0.7~0.9
    v1 = ((1 - e) * M + 2 * m) * v2_ / ((1 + e) * M)
    return v1


def force_point(n, f, o, p):   # arg:n 标准化后的f的方向向量    返回在x,y,z轴的受力和力矩   o:圆心 p:圆上作用点
    cosZ = n[2]
    cosY = n[1]
    cosX = n[0]
    # 分力
    fx = cosX * f
    fy = cosY * f
    fz = cosZ * f
    # 分力矩
    op = np.array([p[0] - o[0], p[1] - o[1], p[2] - o[2]])
    fZ, fY, fX = np.array([f * n[0], f * n[1], 0]), np.array([f * n[0], 0, f * n[2]]), np.array([0, f * n[1], f * n[2]])
    opZ, opY, opX = np.array([op[0], op[1], 0]), np.array([op[0], 0, op[2]]), np.array([0, op[1], op[2]])
    tZ, tY, tX = np.cross(opZ, fZ), np.cross(opY, fY), np.cross(opX, fX)

    return [fx, fy, fz], [tX, tY, tZ]

def micro_state(deltaT, F, T, o, v, L, num ,r ,R):
    J = 0.25 * 3.6 * (r**2 + R**2)
    if num == 10:
        pass
    if num == 8:
        [Fx, Fy, Fz] = np.sum(F, axis=0)
        [Tx, Ty, Tz] = np.sum(F, axis=0)
        print([Fx, Fy, Fz])
        print([Tx, Ty, Tz])
        Fz = Fz - 3.6 * 9.8



'''
x,y,z = circle_curving(np.array([1,0,1]), 20, np.array([1,1,-1]))
aa = circle_point(x, y, z, 8)
listX = []
listY = []
listZ = []
for i in np.arange(0, np.size(aa,0)):
    xx, yy = get_angle([1, 1, -1], i, 8, 170)
    listX.append(xx)
    listY.append(yy)
plt.scatter(aa[:,0], aa[:,1])
plt.scatter(np.array(listX).T,np.array(listY).T)
plt.show()
'''

#plt.scatter(circle_point(x,y,z,8)[:,0],circle_point(x,y,z,8)[:,1])
#plt.scatter(circle_point(x,y,z,8)[2,0],circle_point(x,y,z,8)[2,1])
#plt.show()
x, y, z = circle_curving(np.array([0, 0, 1]), 0.2, np.array([0, 0, -0.11]))
pList = circle_point(x, y, z, 8)
nList = []
fList = []
tList = np.array([0, 0, 0])
for i in np.arange(0, np.size(pList, 0), 1):
    nList.append(get_angle(pList[i], i, 8, 1.7, pList)[1])
    aa = force_point(nList[i], 80, [0, 0, -0.11], pList[i])
    fList.append(aa[0])
    tList = np.vstack((tList, aa[1]))
F = np.array(fList)
T = np.array(tList[1: np.size(tList, axis=0), :])

N = nList[0]
for i in np.arange(1, np.size(nList, 0)):
    N = np.vstack((N, nList[i]))

#print(N)
#print(F)
#print(T)
print(np.sum(F, 0))
print(np.sum(T, 0))


plt.scatter(T[0:23, 0], T[0:23, 1])
#plt.scatter(F[:,0],F[:,1])
#plt.scatter(N[:, 0], N[:, 1])
plt.show()
