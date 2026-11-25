# 冻结拓扑（机器人不动、没有控制器）下的三组实验，给海报用
#   E1 通信半径扫描：半径 R 变化时，连通比例、指令到达率、λ2、实测轮数 / 直径 / 谱上界
#   E2 flag 与 consensus 对比：同一张图上，两种传播方式需要多少轮（检验论文 A1 / §6）
#   E3 节点损失：随机删 15% 与区域性删 15%（模拟局部塌方），连通性与到达率
# 运行：D:\python\python.exe static_experiments.py   → results/static_experiments.json
# 注意：全部是静态快照检查，不是论文 §5 的评估（那需要运动 + CBF 控制器）。

import json
import math
import random
import numpy as np
import networkx as nx
import sanity_check_200 as s          # 复用 make_graph / spectral_bound / flood_rounds

N, D_MAX, SEEDS = s.N, s.D_MAX, list(s.SEEDS)
BETA = 1.0 / (D_MAX + 1)              # consensus 步长；β < 1/最大度 保证稳定且数值不越界
THRESH = 0.95                         # consensus "响应" 阈值：所有机器人 z ≥ 0.95


def source_of(pos):
    """离操作终端 (0,0) 最近的机器人"""
    return min(range(len(pos)), key=lambda i: math.dist(pos[i], (0, 0)))


def laplacian_eigs(G):
    L = nx.laplacian_matrix(G, nodelist=sorted(G.nodes)).toarray().astype(float)
    return np.linalg.eigvalsh(L)


def consensus_rounds(G, src, thresh=THRESH, beta=BETA, kmax=10**7):
    """单个 stubborn 源 z_s = 1，其他机器人从 0 开始：
         z_i ← z_i + β Σ_j a_ij (z_j − z_i)
       误差 e_k = (I − β L_g)^k e_0（L_g = 去掉源节点的 grounded Laplacian）。
       用特征分解直接算第 k 轮的误差，二分查找"所有机器人 z ≥ thresh"的最小 k。"""
    nodes = sorted(G.nodes)
    L = nx.laplacian_matrix(G, nodelist=nodes).toarray().astype(float)
    keep = [k for k, v in enumerate(nodes) if v != src]
    Lg = L[np.ix_(keep, keep)]
    mu, V = np.linalg.eigh(Lg)
    c = V.T @ (-np.ones(len(keep)))           # e_0 = 0 − 1
    lam = 1 - beta * mu                        # 每个模态每轮的衰减系数

    def worst_err(k):
        return np.max(np.abs(V @ (c * lam ** k)))

    if worst_err(kmax) > 1 - thresh:
        return None, float(mu[0])
    lo, hi = 0, 1
    while worst_err(hi) > 1 - thresh:
        lo, hi = hi, hi * 2
    while hi - lo > 1:
        mid = (lo + hi) // 2
        lo, hi = (lo, mid) if worst_err(mid) <= 1 - thresh else (mid, hi)
    return hi, float(mu[0])


def describe(G, src):
    """连通性、到达率；连通时再算 λ2、λN、实测轮数、离心率、直径、谱上界"""
    reach = len(nx.node_connected_component(G, src)) / G.number_of_nodes()
    row = {"connected": nx.is_connected(G), "reach": reach}
    if row["connected"]:
        e = laplacian_eigs(G)
        row.update(lam2=float(e[1]), lamN=float(e[-1]),
                   flood=s.flood_rounds(G, src), ecc=nx.eccentricity(G, src),
                   diam=nx.diameter(G), bound=s.spectral_bound(e[1], e[-1], G.number_of_nodes()))
    return row


out = {"N": N, "d_max": D_MAX, "seeds": SEEDS, "beta": BETA, "threshold": THRESH}

# ---------------- E1：通信半径扫描 ----------------
RADII = [1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0]
e1 = []
for R in RADII:
    s.R = R
    for seed in SEEDS:
        G, pos = s.make_graph(seed)
        e1.append({"R": R, "seed": seed, **describe(G, source_of(pos))})
    rows = [r for r in e1 if r["R"] == R]
    ok = [r for r in rows if r["connected"]]
    print(f"E1 R={R}: 连通 {len(ok)}/20, 平均到达 {np.mean([r['reach'] for r in rows]):.0%}"
          + (f", λ2 {np.median([r['lam2'] for r in ok]):.3f}, 实测 {np.median([r['flood'] for r in ok]):.0f}"
             f", 上界 {np.median([r['bound'] for r in ok]):.0f}" if ok else ""))
out["E1_radius_sweep"] = e1

# ---------------- E2：flag 与 consensus ----------------
# R ≥ 3.5 m 时每个机器人的 8 条链路已经"用满"，生成的图与 3.5 m 完全相同，所以只取 2.5/3.0/3.5，并按边集去重
e2, seen = [], set()
for R in [2.5, 3.0, 3.5]:
    s.R = R
    for seed in SEEDS:
        G, pos = s.make_graph(seed)
        key = frozenset(map(frozenset, G.edges))
        if not nx.is_connected(G) or key in seen:
            continue
        seen.add(key)
        src = source_of(pos)
        k_cons, mu_g = consensus_rounds(G, src)
        e2.append({"R": R, "seed": seed, "lam2": float(laplacian_eigs(G)[1]),
                   "flag": s.flood_rounds(G, src), "consensus": k_cons,
                   "lam_grounded": mu_g, "d_s": G.degree(src),
                   "grounded_upper": G.degree(src) / (N - 1)})
fl = np.array([r["flag"] for r in e2]); co = np.array([r["consensus"] for r in e2], dtype=float)
print(f"E2: {len(e2)} 张连通图; flag 中位数 {np.median(fl):.0f} 轮, consensus 中位数 {np.median(co):.0f} 轮, "
      f"比值中位数 {np.median(co / fl):.0f}×")
print("E2 检查 λ(L_g) ≤ d_s/(N−1) 全部成立：", all(r["lam_grounded"] <= r["grounded_upper"] + 1e-12 for r in e2))
lam2 = np.array([r["lam2"] for r in e2])
def rank_corr(a, b):
    ra, rb = np.argsort(np.argsort(a)), np.argsort(np.argsort(b))
    return float(np.corrcoef(ra, rb)[0, 1])
lg = np.array([r["lam_grounded"] for r in e2]); ds = np.array([r["d_s"] for r in e2])
print(f"E2 秩相关（越接近 ±1 越相关）：flag~λ2 {rank_corr(lam2, fl):+.2f}, consensus~λ2 {rank_corr(lam2, co):+.2f}, "
      f"consensus~λ(L_g) {rank_corr(lg, co):+.2f}, consensus~d_s {rank_corr(ds, co):+.2f}")
pred = np.log(1 / (1 - THRESH)) / (BETA * lg)
print(f"E2 实测 consensus 轮数 / 预测 ln(20)/(β·λ(L_g))：{np.min(co / pred):.2f}–{np.max(co / pred):.2f}")
out["E2_flag_vs_consensus"] = e2

# ---------------- E3：节点损失 ----------------
s.R = 3.0
LOSS = 0.15
e3 = []
for seed in SEEDS:
    G, pos = s.make_graph(seed)
    src = source_of(pos)
    if not nx.is_connected(G):
        continue                                     # 只看损失前连通的图
    k = int(round(LOSS * N))
    rng = random.Random(1000 + seed)
    rand_removed = rng.sample([i for i in G.nodes if i != src], k)
    # 区域性损失：东北区域（离 (15,15) 最近的 30 个机器人）——对应论文场景里的二次塌方
    region_removed = sorted((i for i in G.nodes if i != src), key=lambda i: math.dist(pos[i], (15, 15)))[:k]
    for mode, removed in (("none", []), ("random", rand_removed), ("region", region_removed)):
        H = G.copy()
        H.remove_nodes_from(removed)                # 剩下的机器人不重新建链（冻结拓扑）
        e3.append({"seed": seed, "mode": mode, "removed": len(removed), **describe(H, src)})
for mode in ("none", "random", "region"):
    rows = [r for r in e3 if r["mode"] == mode]
    ok = [r for r in rows if r["connected"]]
    print(f"E3 {mode:6s}: 连通 {len(ok)}/{len(rows)}, 平均到达 {np.mean([r['reach'] for r in rows]):.1%}"
          + (f", 实测轮数中位数 {np.median([r['flood'] for r in ok]):.0f}, λ2 中位数 {np.median([r['lam2'] for r in ok]):.3f}" if ok else ""))
out["E3_node_loss"] = e3

with open("results/static_experiments.json", "w", encoding="utf-8") as f:
    json.dump(out, f, indent=1)
print("已保存 results/static_experiments.json")
