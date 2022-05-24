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

def L1(m, a, b, c):
    m.addConstr(c == a)
    m.addConstr(c == b)

def L2(m, a, b, c):
    m.addConstr(c >= b)
    m.addConstr(c == a + b)

def L3(m, a, b, c, d):
    m.addConstr(d == a + b + c)
    m.addConstr(d >= a + b)

def GenerateJModel(r):
    # New MILP model M
    m = Model("JSuperPoly")

    # M.var ← u^0[0..31]
    # M.var ← K[0..63]
    u0 = []
    K = []
    flag = []
    for i in range(32):
        u0.append(m.addVar(lb=0.0, ub=1.0, vtype=GRB.BINARY, name="u^0_" + str(i)))
        K.append(m.addVar(lb=0.0, ub=1.0, vtype=GRB.BINARY, name="K_" + str(2 * i)))
        flag.append(m.addVar(lb=0.0, ub=2.0, vtype=GRB.INTEGER, name="flag_" + str(2 * i)))
        K.append(m.addVar(lb=0.0, ub=1.0, vtype=GRB.BINARY, name="K_" + str(2 * i + 1)))
        flag.append(m.addVar(lb=0.0, ub=2.0, vtype=GRB.INTEGER, name="flag_" + str(2 * i + 1)))
    u = u0

    # for i = 1; i <= r; i ← i + 1 do
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
            L1(m, RK[j], u[j], a[j])
            if con_i[j] == 1:
                L2(m, u[j + 16], constant_i_var[0], b[j])
                del(constant_i_var[0])
            else:
                m.addConstr(u[j + 16] == b[j])
            L2(m, a[j], b[j], c[j])
            if j < 14:
                L1(m, c[j + 1], c[j + 2], d[j])
            if j < 14:
                L2(m, c[j], d[j], e[j])
            else:
                m.addConstr(e[j] == c[j])
            L3(m, e[(j + 3) % 16], e[(j + 9) % 16], e[(j + 14) % 16], f[j])
            L1(m, f[j], RK[j], g[j])
            L2(m, g[j], u[16 + j], L[j])
            L2(m, f[j], u[j], R[j])
        if i % 2 == 1 and i == r: # 奇数轮次的最后一轮，例如7轮加密中的第7轮，没有取反的第8轮
            KeyScheduleSingle(m, i)
        elif i % 2 == 0:
            KeySchedule(m, i)
        u = L + R

    return m


def main():
    m = GenerateJModel(8)
    m.write("angiejc.lp")
    m.optimize()

if __name__ == '__main__':
    main()