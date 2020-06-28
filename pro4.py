#!/usr/bin/env python
# _*_ coding:utf-8 _*_
# 文件: pro4.py
import numpy as np
import geatpy as ea
import time

def circle_curving(n, r, o):  #绘制圆轨迹
    u = np.array([n[2], 0, -n[0]])
    v = np.cross(u, n)
    u_std = u / np.sqrt(u[0] ** 2 + u[1] ** 2 + u[2] ** 2)
    v_std = v / np.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)
    theta = np.linspace(0, np.pi * 2, 721)
    x = o[0] + r * (u_std[0] * np.cos(theta) + v_std[0] * np.sin(theta))
    y = o[1] + r * (u_std[1] * np.cos(theta) + v_std[1] * np.sin(theta))
    z = o[2] + r * (u_std[2] * np.cos(theta) + v_std[2] * np.sin(theta))
    return x, y, z


def circle_point(x, y, z, num, index):   #绘制num个受力点
    if num == 8:
        p_x = x[index:719:90]
        p_y = y[index:719:90]
        p_z = z[index:719:90]
        res = np.array([p_x, p_y, p_z]).T
        return res
    elif num == 10:
        p_x = x[index:719:72]
        p_y = y[index:719:72]
        p_z = z[index:719:72]
        res = np.array([p_x, p_y, p_z]).T
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
    vector = np.array([np.cos(-seta0 / 180 * np.pi), np.sin(-seta0 / 180 * np.pi)])
    part = - 360 / num
    for ii in np.arange(1, num):
        vector = np.vstack((vector, np.array( [np.cos((-seta0 + part * ii) / 180 * np.pi), np.sin((-seta0 + part * ii) / 180 * np.pi)] )))
    lenZ = np.tan(np.arcsin(depth / length))
    nf = np.array([vector[index, 0], vector[index, 1], lenZ])
    return nf / np.linalg.norm(nf)


def force_point(nf, f, o, p):   #  arg:n 标准化后的f的方向向量    返回在x,y,z轴的受力和力矩   o:圆心 p:圆上作用点
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


def get_N0_nf_F_T(FF):
    o = np.array([0, 0, -0.15])
    n = np.array([0., 0., 1.])
    G = 3.6 * 9.8
    seta0 = 24
    addr = 48
    r = 0.2
    num = 10
    length = 2.0
    depth = 0.15

    x, y, z = circle_curving(n, r, o)
    N0 = circle_point(x, y, z, num, addr)
    nflist = []
    F = []
    T = []
    for i in np.arange(0, num):
        nflist.append(get_angle(i, num, length, depth, seta0))
        FT = force_point(nflist[i], FF[i], o, N0[i])
        F.append(FT[0])
        T.append(FT[1])

    nflist = np.array(nflist)
    Ttotal = np.sum(T, 0)
    Ftotal = np.sum(F, 0)
    Ftotal[2] -= G
    
    return N0, nflist, Ftotal, Ttotal


def iter_state(Ttotal, Ftotal, v, w, deltaT, o, beta):
    J = 0.5 * 3.6 * 0.04 + 1 / 12 * 3.6 * (0.22 ** 2)

    o += v * deltaT
    beta += deltaT * w
    w += (Ttotal * (deltaT / J))[1]
    a = Ftotal / 3.6
    v += a * deltaT
    return v, w, o, beta

def point_move(N, rx, ry, o, beta):
    for i in np.arange(0, np.size(N, 0)):
        N[i] = [o[0] + rx[i] * np.cos(beta), o[1] + ry[i], o[2] - rx[i] * np.sin(beta)]
    return N

def func(FF):
    beta = 0
    v = np.array([0., 0., 0.])
    w = 0
    o = np.array([0, 0, -0.15])
    G = 3.6 * 9.8
    deltaT = 0.0002
    endtime = 0.24
    times = int(endtime / deltaT)
    num = 10

    N0, nflist, Ftotal0, Ttotal0 = get_N0_nf_F_T(FF)
    N = N0
    Ftotal = Ftotal0
    Ttotal = Ttotal0
    rx = np.array(N0[:, 0])
    ry = np.array(N0[:, 1])
    set_beta = []
    for i in np.linspace(0, endtime, times + 1):
        v, w, o, beta = iter_state(Ttotal, Ftotal, v, w, deltaT, o, beta)
        N = point_move(N, rx, ry, o, beta)
        F = []
        T = []
        for i in np.arange(0, num):
            FT = force_point(nflist[i], FF[i], o, N[i])
            F.append(FT[0])
            T.append(FT[1])

        Ttotal = np.sum(T, 0)
        Ftotal = np.sum(F, 0)
        Ftotal[2] -= G
        Ftotal[1] = Ttotal[0] = Ttotal[2] = 0
        set_beta.append(beta)
    return beta * 180 / np.pi, v[0], Ttotal0[0], Ttotal0[2]

def minfun(FF):
    F0 = FF[:, [0]]
    F1 = FF[:, [1]]
    F2 = FF[:, [2]]
    F3 = FF[:, [3]]
    F4 = FF[:, [4]]
    F5 = FF[:, [5]]
    F6 = FF[:, [6]]
    F7 = FF[:, [7]]
    F8 = FF[:, [8]]
    F9 = FF[:, [9]]

    res = np.array([func(np.array([F0[0,0], F1[0,0], F2[0,0], F3[0,0], F4[0,0], F5[0,0], F6[0,0], F7[0,0], F8[0,0], F9[0,0]]))])
    len = np.size(F0, 0)
    for i in np.arange(1, len):
        tmp = np.array([func(np.array([F0[i, 0], F1[i, 0], F2[i, 0], F3[i, 0], F4[i, 0], F5[i, 0], F6[i, 0], F7[i, 0], F8[i, 0], F9[i, 0]])) ])
        res = np.vstack((res, tmp))

    """============================目标函数============================"""
    print(res)
    aim = np.fabs(res[:, [0]] - 0.524 * np.ones((len, 1))) * 500 + np.fabs(res[:, [2]]) * 300 + np.fabs(res[:, [3]]) * 20
    return aim


def main():

    """============================变量设置============================"""
    x1 = [50, 100]  # 第一个决策变量范围
    x2 = [50, 100]  # 第二个决策变量范围
    x3 = [50, 100]
    x4 = [50, 100]
    x5 = [50, 100]
    x6 = [50, 100]
    x7 = [50, 100]
    x8 = [50, 100]
    x9 = [50, 100]
    x10 = [50, 100]


    b1 = [1, 1]  # 第一个决策变量边界，1表示包含范围的边界，0表示不包含
    b2 = [1, 1]  # 第二个决策变量边界，1表示包含范围的边界，0表示不包含
    b3 = [1, 1]
    b4 = [1, 1]
    b5 = [1, 1]
    b6 = [1, 1]
    b7 = [1, 1]
    b8 = [1, 1]
    b9 = [1, 1]
    b10 = [1, 1]


    ranges = np.vstack([x1, x2, x3, x4, x5, x6, x7, x8, x9, x10]).T  # 生成自变量的范围矩阵，使得第一行为所有决策变量的下界，第二行为上界
    borders = np.vstack([b1, b2, b3, b4, b5, b6, b7, b8, b9, b10]).T  # 生成自变量的边界矩阵
    varTypes = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0])  # 决策变量的类型，0表示连续，1表示离散
    """==========================染色体编码设置========================="""
    Encoding = 'BG'  # 'BG'表示采用二进制/格雷编码
    codes = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]  # 决策变量的编码方式，设置两个0表示两个决策变量均使用二进制编码
    precisions = [4, 4, 4, 4, 4, 4, 4, 4, 4, 4]  # 决策变量的编码精度，表示二进制编码串解码后能表示的决策变量的精度可达到小数点后6位
    scales = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]  # 0表示采用算术刻度，1表示采用对数刻度
    FieldD = ea.crtfld(Encoding, varTypes, ranges, borders, precisions, codes, scales)  # 调用函数创建译码矩阵
    """=========================遗传算法参数设置========================"""
    NIND = 30  # 种群个体数目
    MAXGEN = 30  # 最大遗传代数
    maxormins = [1]  # 列表元素为1则表示对应的目标函数是最小化，元素为-1则表示对应的目标函数是最大化
    selectStyle = 'rws'  # 采用轮盘赌选择
    recStyle = 'xovdp'  # 采用两点交叉
    mutStyle = 'mutbin'  # 采用二进制染色体的变异算子
    pc = 0.7  # 交叉概率
    pm = 1  # 整条染色体的变异概率（每一位的变异概率=pm/染色体长度）
    Lind = int(np.sum(FieldD[0, :]))  # 计算染色体长度
    obj_trace = np.zeros((MAXGEN, 2))  # 定义目标函数值记录器
    var_trace = np.zeros((MAXGEN, Lind))  # 染色体记录器，记录历代最优个体的染色体
    """=========================开始遗传算法进化========================"""
    start_time = time.time()  # 开始计时
    Chrom = ea.crtpc(Encoding, NIND, FieldD)  # 生成种群染色体矩阵
    variable = ea.bs2real(Chrom, FieldD)  # 对初始种群进行解码
    ObjV = minfun(variable)  # 计算初始种群个体的目标函数值
    FitnV = ea.ranking(maxormins * ObjV)  # 根据目标函数大小分配适应度值
    best_ind = np.argmax(FitnV)  # 计算当代最优个体的序号
    # 开始进化
    for gen in range(MAXGEN):
        print(gen)
        SelCh = Chrom[ea.selecting(selectStyle, FitnV, NIND - 1), :]  # 选择
        SelCh = ea.recombin(recStyle, SelCh, pc)  # 重组
        SelCh = ea.mutate(mutStyle, Encoding, SelCh, pm)  # 变异
        # 把父代精英个体与子代的染色体进行合并，得到新一代种群
        Chrom = np.vstack([Chrom[best_ind, :], SelCh])
        Phen = ea.bs2real(Chrom, FieldD)  # 对种群进行解码(二进制转十进制)
        print(Phen)
        ObjV = minfun(Phen)  # 求种群个体的目标函数值
        FitnV = ea.ranking(maxormins * ObjV)  # 根据目标函数大小分配适应度值
        # 记录
        best_ind = np.argmax(FitnV)  # 计算当代最优个体的序号
        obj_trace[gen, 0] = np.sum(ObjV) / ObjV.shape[0]  # 记录当代种群的目标函数均值
        obj_trace[gen, 1] = ObjV[best_ind]  # 记录当代种群最优个体目标函数值
        var_trace[gen, :] = Chrom[best_ind, :]  # 记录当代种群最优个体的染色体
        print(best_ind)
    # 进化完成
    end_time = time.time()  # 结束计时
    ea.trcplot(obj_trace, [['种群个体平均目标函数值', '种群最优个体目标函数值']])  # 绘制图像
    """============================输出结果============================"""
    best_gen = np.argmin(obj_trace[:, [1]])
    print('最优解的目标函数值：', obj_trace[best_gen, 1])
    variable = ea.bs2real(var_trace[[best_gen], :], FieldD)  # 解码得到表现型（即对应的决策变量值）
    print('最优解的决策变量值为：')
    print(variable)
    print('用时：', end_time - start_time, '秒')

main()


