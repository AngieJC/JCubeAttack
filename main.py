from gurobipy import *

def KeySchedule(m, r): # 利用主密钥k生成第r轮的轮密钥并添加到模型中
    m.update()
    if r / 2 <= 4:
        for i in range(16):
            RK_i_2 = m.getVarByName("RK^" + str(r - 2) + "_" + str(i))
            RK_i_1 = m.getVarByName("RK^" + str(r - 1) + "_" + str(i))
            K_i_1_2 = m.getVarByName("K_" + str(int(((r - 1) / 2)) * 16 + i))
            flag = m.getVarByName("flag_" + str(int(((r - 1) / 2)) * 16 + i))
            m.addConstr(RK_i_2 - RK_i_1 - 3 * K_i_1_2 + 2 * flag == 0)
            m.addConstr(RK_i_1 + K_i_1_2 - flag >= 0)
            m.addConstr(2 * K_i_1_2 - flag >= 0)
            m.addConstr(- K_i_1_2 + flag >= 0)
            m.addConstr(- RK_i_1 - 2 * K_i_1_2 + flag + 1 >= 0)

def KeyScheduleSingle(m, r):
    m.update()
    if r / 2 <= 4:
        for i in range(16):
            RK_i_1 = m.getVarByName("RK^" + str(r - 1) + "_" + str(i))
            K_i_1_2 = m.getVarByName("K_" + str(int(((r - 1) / 2)) * 16 + i))
            flag = m.getVarByName("flag_" + str(int(((r - 1) / 2)) * 16 + i))
            m.addConstr(RK_i_1 == K_i_1_2)
            m.addConstr(RK_i_1 == flag)

# AND
def L1(m, a, b, c, copys):
    if len(copys) == 0:
        m.addConstr(c == a)
        m.addConstr(c == b)
    elif len(copys) == 2:
        m.addConstr(copys[0] <= a)
        m.addConstr(copys[1] <= a)
        m.addConstr(a <= copys[0] + copys[1])
        m.addConstr(c == copys[1])
        m.addConstr(c == b)
    elif len(copys) == 4:
        m.addConstr(copys[0] <= a)
        m.addConstr(copys[1] <= a)
        m.addConstr(a <= copys[0] + copys[1])
        m.addConstr(copys[2] <= b)
        m.addConstr(copys[3] <= b)
        m.addConstr(b <= copys[2] + copys[3])
        m.addConstr(c == copys[1])
        m.addConstr(c == copys[3])

# XOR
def L2(m, a, b, c, copys):
    if len(copys) == 0:
        m.addConstr(c >= b)
        m.addConstr(c == a + b)
    elif len(copys) == 2:
        m.addConstr(copys[0] <= a)
        m.addConstr(copys[1] <= a)
        m.addConstr(a <= copys[0] + copys[1])
        m.addConstr(c >= b)
        m.addConstr(c == copys[1] + b)
    elif len(copys) == 4:
        m.addConstr(copys[2] <= b)
        m.addConstr(copys[3] <= b)
        m.addConstr(b <= copys[2] + copys[3])
        m.addConstr(c >= copys[3])
        m.addConstr(c == copys[1] + copys[3])

# 3比特XOR
def L3(m, a, b, c, d, copys):
    m.addConstr(copys[0] <= a)
    m.addConstr(copys[1] <= a)
    m.addConstr(a <= copys[0] + copys[1])
    m.addConstr(copys[2] <= b)
    m.addConstr(copys[3] <= b)
    m.addConstr(b <= copys[2] + copys[3])
    m.addConstr(copys[4] <= c)
    m.addConstr(copys[5] <= c)
    m.addConstr(c <= copys[4] + copys[5])
    m.addConstr(d == copys[1] + copys[3] + copys[5])
    m.addConstr(d >= copys[1] + copys[3])

def GenerateJModel(r):
    m = Model("JSuperPoly")
    L0 = []
    R0 = []
    K = []
    flag = []
    for i in range(16):
        L0.append(m.addVar(lb=0.0, ub=1.0, vtype=GRB.BINARY, name="L^0_" + str(i)))
        R0.append(m.addVar(lb=0.0, ub=1.0, vtype=GRB.BINARY, name="R^0_" + str(i)))
    for i in range(64):
        K.append(m.addVar(lb=0.0, ub=1.0, vtype=GRB.BINARY, name="K_" + str(i)))
        flag.append(m.addVar(lb=0.0, ub=2.0, vtype=GRB.INTEGER, name="flag_" + str(i)))
    u = L0 + R0

    for i in range(1, r + 1):
        RK = []
        a = []
        b = []
        constant_i_var = []
        c = []
        d = []
        e = []
        f = []
        g = []
        L = []
        R = []
        con_i = [0 for k in range(16)]
        i_bin_str = bin(i - 1).replace("0b", "")
        for k in range(len(i_bin_str)):
            con_i[k + 16 - len(i_bin_str)] = int(i_bin_str)
        for j in range(16):
            RK.append(m.addVar(lb=0.0, ub=1.0, vtype=GRB.BINARY, name="RK^" + str(i - 1) + "_" + str(j)))
            a.append(m.addVar(lb=0.0, ub=1.0, vtype=GRB.BINARY, name="a^" + str(i) + "_" + str(j)))
            b.append(m.addVar(lb=0.0, ub=1.0, vtype=GRB.BINARY, name="b^" + str(i) + "_" + str(j)))
            if con_i[j] == 1:
                constant_i_var.append(m.addVar(lb=0.0, ub=1.0, vtype=GRB.BINARY, name="i^" + str(i) + "_" + str(j))) # 这里的i需要处理
            c.append(m.addVar(lb=0.0, ub=1.0, vtype=GRB.BINARY, name="c^" + str(i) + "_" + str(j)))
            e.append(m.addVar(lb=0.0, ub=1.0, vtype=GRB.BINARY, name="e^" + str(i) + "_" + str(j)))
            f.append(m.addVar(lb=0.0, ub=1.0, vtype=GRB.BINARY, name="f^" + str(i) + "_" + str(j)))
            g.append(m.addVar(lb=0.0, ub=1.0, vtype=GRB.BINARY, name="g^" + str(i) + "_" + str(j)))
            L.append(m.addVar(lb=0.0, ub=1.0, vtype=GRB.BINARY, name="L^" + str(i) + "_" + str(j)))
            R.append(m.addVar(lb=0.0, ub=1.0, vtype=GRB.BINARY, name="R^" + str(i) + "_" + str(j)))
        for j in range(14):
            d.append(m.addVar(lb=0.0, ub=1.0, vtype=GRB.BINARY, name="d^" + str(i) + "_" + str(j)))
        for j in range(16):
            RK_copy1 = m.addVar(lb=0.0, ub=1.0, vtype=GRB.BINARY, name="RK^" + str(i - 1) + "_" + str(j) + "_copy1")
            RK_copy2 = m.addVar(lb=0.0, ub=1.0, vtype=GRB.BINARY, name="RK^" + str(i - 1) + "_" + str(j) + "_copy2")
            L_copy1 = m.addVar(lb=0.0, ub=1.0, vtype=GRB.BINARY, name="L^" + str(i - 1) + "_" + str(j) + "_copy1")
            L_copy2 = m.addVar(lb=0.0, ub=1.0, vtype=GRB.BINARY, name="L^" + str(i - 1) + "_" + str(j) + "_copy2")
            L1(m, RK[j], u[j], a[j], [RK_copy1, RK_copy2, L_copy1, L_copy2])
            RK[j] = RK_copy1
            u[j] = L_copy1

            R_copy1 = m.addVar(lb=0.0, ub=1.0, vtype=GRB.BINARY, name="R^" + str(i - 1) + "_" + str(j) + "_copy1")
            R_copy2 = m.addVar(lb=0.0, ub=1.0, vtype=GRB.BINARY, name="R^" + str(i - 1) + "_" + str(j) + "_copy2")
            if con_i[j] == 1:
                L2(m, u[j + 16], constant_i_var[0], b[j], [R_copy1, R_copy2])
                del(constant_i_var[0])
            else:
                m.addConstr(R_copy2 == b[j])
            u[j + 16] = R_copy1

            L2(m, a[j], b[j], c[j], [])

            if j < 14:
                c_1_copy1 = m.addVar(lb=0.0, ub=1.0, vtype=GRB.BINARY, name="c^" + str(i) + "_1_" + str(j) + "_copy1")
                c_1_copy2 = m.addVar(lb=0.0, ub=1.0, vtype=GRB.BINARY, name="c^" + str(i) + "_1_" + str(j) + "_copy2")
                c_2_copy1 = m.addVar(lb=0.0, ub=1.0, vtype=GRB.BINARY, name="c^" + str(i) + "_2_" + str(j) + "_copy1")
                c_2_copy2 = m.addVar(lb=0.0, ub=1.0, vtype=GRB.BINARY, name="c^" + str(i) + "_2_" + str(j) + "_copy2")
                L1(m, c[j + 1], c[j + 2], d[j], [c_1_copy1, c_1_copy2, c_2_copy1, c_2_copy2])
                c[j + 1] = c_1_copy1
                c[j + 2] = c_2_copy1
            if j < 14:
                L2(m, c[j], d[j], e[j], [])
            else:
                m.addConstr(e[j] == c[j])

        for j in range(16):
            e_3_copy1 = m.addVar(lb=0.0, ub=1.0, vtype=GRB.BINARY, name="e^" + str(i) + "_3_" + str(j) + "_copy1")
            e_3_copy2 = m.addVar(lb=0.0, ub=1.0, vtype=GRB.BINARY, name="e^" + str(i) + "_3_" + str(j) + "_copy2")
            e_9_copy1 = m.addVar(lb=0.0, ub=1.0, vtype=GRB.BINARY, name="e^" + str(i) + "_9_" + str(j) + "_copy1")
            e_9_copy2 = m.addVar(lb=0.0, ub=1.0, vtype=GRB.BINARY, name="e^" + str(i) + "_9_" + str(j) + "_copy2")
            e_14_copy1 = m.addVar(lb=0.0, ub=1.0, vtype=GRB.BINARY, name="e^" + str(i) + "_14_" + str(j) + "_copy1")
            e_14_copy2 = m.addVar(lb=0.0, ub=1.0, vtype=GRB.BINARY, name="e^" + str(i) + "_14_" + str(j) + "_copy2")
            L3(m, e[(j + 3) % 16], e[(j + 9) % 16], e[(j + 14) % 16], f[j], [e_3_copy1, e_3_copy2, e_9_copy1, e_9_copy2, e_14_copy1, e_14_copy2])
            e[(j + 3) % 16] = e_3_copy1
            e[(j + 9) % 16] = e_9_copy1
            e[(j + 14) % 16] = e_14_copy1

        for j in range(16):
            f_copy1 = m.addVar(lb=0.0, ub=1.0, vtype=GRB.BINARY, name="f^" + str(i) + "_" + str(j) + "_copy1")
            f_copy2 = m.addVar(lb=0.0, ub=1.0, vtype=GRB.BINARY, name="f^" + str(i) + "_" + str(j) + "_copy2")
            L1(m, f[j], RK[j], g[j], [f_copy1, f_copy2])
            f[j] = f_copy1

            L2(m, g[j], u[16 + j], L[j], [])
            L2(m, f[j], u[j], R[j], [])

        if i % 2 == 1 and i == r: # 奇数轮次的最后一轮，例如7轮加密中的第7轮，没有取反的第8轮
            KeyScheduleSingle(m, i)
        elif i % 2 == 0:
            KeySchedule(m, i)
        u = L + R

    return m


def SearchDegree(m, r, I):
    for i in range(32):
        if i < 16:
            u0 = m.getVarByName("L^0_" + str(i))
        else:
            u0 = m.getVarByName("R^0_" + str(i - 16))
        if i in I:
            # 公开变量(明文)中立方元素全为1
            m.addConstr(u0 == 1)
        else:
            # 非立方元素全为0
            m.addConstr(u0 == 0)

    # L^{r}[i] = g^{r}[i] ^ R^{r-1}[i]
    for i in range(16):
        L = m.getVarByName("L^" + str(r) + "_" + str(i))
        R = m.getVarByName("R^" + str(r) + "_" + str(i))
        if i == 0:
            m.addConstr(L == 1)
            m.addConstr(R == 0)
        else:
            m.addConstr(L == 0)
            m.addConstr(R == 0)


    # 目标函数：令表示密钥的变量求和最大，最大值为超级多项式的次数

    return m


def main(r, I):
    m = GenerateJModel(r)
    m = SearchDegree(m, r, I)
    m.setParam("PoolSearchMode", 2)
    m.setParam("PoolSolutions", 2000000000)
    m.write("angiejc.lp")
    m.optimize()
    '''
    nSolutions = m.SolCount
    for solution in range(nSolutions):
        m.setParam(GRB.Param.SolutionNumber, solution)
        m.write("angiejc_" + str(solution) + ".sol")
        '''

if __name__ == '__main__':
    I = [i for i in range(31)]
    main(1, I)