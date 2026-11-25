# cbf_sim.py 的正确性检查
# 运行：D:\python\python.exe test_cbf_sim.py
import itertools
import numpy as np
import cbf_sim as c

rng = np.random.default_rng(0)

# 1) λ2 解析梯度 vs 有限差分（选 λ2 是单特征值的样本）
X = rng.uniform(0, 5, (40, 2))
N_save = c.N
A, D, diff = c.weights(X)
w, V = c.spectrum(A)
assert w[2] - w[1] > 1e-3 and w[1] > 1e-6, "样本 λ2 不是单特征值或图不连通，换种子"
g = c.lambda2_grad(X, A, D, diff, V[:, 1])
eps, worst = 1e-6, 0.0
for i in range(40):
    for ax in range(2):
        Xp, Xm = X.copy(), X.copy()
        Xp[i, ax] += eps; Xm[i, ax] -= eps
        fd = (c.spectrum(c.weights(Xp)[0])[0][1] - c.spectrum(c.weights(Xm)[0])[0][1]) / (2 * eps)
        worst = max(worst, abs(fd - g[i, ax]))
print(f"梯度最大误差 {worst:.2e}")
assert worst < 1e-5

# 2) 局部 QP 精确解 vs 网格暴力搜索
grid = np.linspace(-c.UMAX, c.UMAX, 401)
G = np.array(list(itertools.product(grid, grid)))
worst_obj = 0.0
for _ in range(200):
    u_nom = rng.normal(0, 0.6, (1, 2)); gi = rng.normal(0, 1, (1, 2)); b = rng.normal(0, 0.3, 1)
    u, act, bad = c.cbf_filter(u_nom, gi, b)
    if bad[0]:
        assert c.UMAX * np.abs(gi).sum() < b[0]
        continue
    assert (gi[0] @ u[0]) >= b[0] - 1e-7 and np.all(np.abs(u) <= c.UMAX + 1e-12)
    feas = G[G @ gi[0] >= b[0]]
    best = ((feas - u_nom[0]) ** 2).sum(1).min()
    obj = ((u[0] - u_nom[0]) ** 2).sum()
    worst_obj = max(worst_obj, obj - best)
print(f"QP 目标值比网格最优多出 ≤ {worst_obj:.2e}（网格步长 {grid[1]-grid[0]:.4f}）")
assert worst_obj < 1e-5

# 3) 谱上界在随机加权图上 ≥ 支撑图直径
import networkx as nx
fails = 0
for s in range(300):
    Xs = rng.uniform(0, rng.uniform(3, 12), (60, 2))
    A, _, _ = c.weights(Xs)
    G2 = nx.from_numpy_array((A > 0).astype(int))
    if not nx.is_connected(G2):
        continue
    lam2 = c.spectrum(A)[0][1]
    old = c.N; c.N = 60
    B, _ = c.bound_now(A, lam2)
    c.N = old
    fails += nx.diameter(G2) > B
print("上界失败次数：", fails)
assert fails == 0
print("全部检查通过")
