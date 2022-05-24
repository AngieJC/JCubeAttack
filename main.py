from gurobipy import *

def KeySchedule(m, r): # 利用主密钥k生成第r轮的轮密钥并添加到模型中
    if 0 <= (r - 1) / 2 < 4:
        if (r - 1) % 2 == 0:
            return

def L1(m, a, b, c):
    m.addConstr(c == a)
    m.addConstr(c == b)

def L2(m, a, b, c):
    m.addConstr(c >= a)
    m.addConstr(c >= b)
    m.addConstr(c >= a + b)

def L3(m, a, b, c, d):
    m.addConstr(d >= a)
    m.addConstr(d >= b)
    m.addConstr(d >= c)
    m.addConstr(d >= a + b + c)

def GenerateJModel(R):
    # New MILP model M
    m = Model("JSuperPoly")

    # M.var ← u^0[0..31]
    # M.var ← K[0..63]
    u0 = []
    for i in range(32):
        u0.append(m.addVar(lb=0.0, ub=1.0, vtype=GRB.BINARY, name="u^0_" + str(i)))
        # k = m.addVar(lb=0.0, ub=1.0, obj=0.0, vtype=GRB.BINARY, name="K_" + str(2 * i))
        # k = m.addVar(lb=0.0, ub=1.0, obj=0.0, vtype=GRB.BINARY, name="K_" + str(2 * i + 1))
    u = u0

    # for i = 1; i <= r; i ← i + 1 do
    for i in range(1, R + 1):
        RK = []
        a = []
        b = []
        con_i = []
        c = []
        d = []
        e = []
        f = []
        g = []
        L = []
        R = []
        for j in range(16):
            RK.append(m.addVar(lb=0.0, ub=1.0, vtype=GRB.BINARY, name="RK^" + str(i - 1) + "_" + str(j)))
            a.append(m.addVar(lb=0.0, ub=1.0, vtype=GRB.BINARY, name="a^" + str(i) + "_" + str(j)))
            b.append(m.addVar(lb=0.0, ub=1.0, vtype=GRB.BINARY, name="b^" + str(i) + "_" + str(j)))
            con_i.append(m.addVar(lb=0.0, ub=1.0, vtype=GRB.BINARY, name="i^" + str(i) + "_" + str(j))) # 这里的i需要处理
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
            L2(m, u[j + 16], con_i[j], b[j])
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
        u = L + R
        # KeySchedule(m, i)

    return m


def main():
    m = GenerateJModel(1)
    m.write("angiejc.lp")
    m.optimize()

if __name__ == '__main__':
    main()