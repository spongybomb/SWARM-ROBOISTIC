# 小验证：在随机生成的 200 机器人通信图上，检查
#     实测传播轮数 = 注入点离心率 ≤ 图直径 ≤ 谱上界 B(G)
# 注意：这是"冻结拓扑"的静态检查，没有运动、没有 CBF 控制器，
#       不是论文 §5 的评估实验，只是一个示例。
# 运行：D:\python\python.exe sanity_check_200.py

import math
import random
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

# ===== 1. 参数（来自 Extended Abstract §1）=====
N = 200          # 机器人数量
SIZE = 20.0      # 场地 20 m × 20 m
R = 3.0          # 通信半径（论文场景是 2 m；改成 2.0 可复现"全部不连通"）
D_MAX = 8        # 每个机器人最多 8 条链路
SEEDS = range(20)  # 20 个随机种子


def make_graph(seed):
    """随机撒点；距离 < R 的才能连；每个机器人只保留最近的 8 个邻居，
    且两边都选了对方才算一条边（这样保证度数 ≤ 8）。"""
    rng = random.Random(seed)
    pos = [(rng.uniform(0, SIZE), rng.uniform(0, SIZE)) for _ in range(N)]
    nearest = []
    for i in range(N):
        cand = []
        for j in range(N):
            if i != j:
                d = math.dist(pos[i], pos[j])
                if d < R:
                    cand.append((d, j))
        cand.sort()
        nearest.append({j for _, j in cand[:D_MAX]})
    G = nx.Graph()
    G.add_nodes_from(range(N))
    for i in range(N):
        for j in nearest[i]:
            if i in nearest[j]:
                G.add_edge(i, j)
    return G, pos


def spectral_bound(lam2, lamN, n):
    """Chung–Faber–Manteuffel 直径上界"""
    return math.floor(math.acosh(n - 1) / math.acosh((lamN + lam2) / (lamN - lam2))) + 1


def flood_rounds(G, source):
    """逐轮模拟 flag 传播：每轮只能往外走一跳（用"本轮/下一轮"两个集合）。"""
    have = {source}
    rounds = 0
    while len(have) < G.number_of_nodes():
        new = set(have)
        for i in have:
            new.update(G.neighbors(i))
        if new == have:          # 传不动了 = 图不连通
            return None
        have = new
        rounds += 1
    return rounds


# ===== 2. 对每个种子做检查 =====
rows = []
for seed in SEEDS:
    G, pos = make_graph(seed)
    connected = nx.is_connected(G)
    # 指令从离操作终端 (0,0) 最近的机器人注入
    source = min(range(N), key=lambda i: math.dist(pos[i], (0, 0)))
    row = {"seed": seed, "connected": connected,
           "max_deg": max(d for _, d in G.degree())}
    if connected:
        L = nx.laplacian_matrix(G).toarray().astype(float)
        eig = np.linalg.eigvalsh(L)
        lam2, lamN = eig[1], eig[-1]
        row.update(lam2=lam2, lamN=lamN,
                   flood=flood_rounds(G, source),
                   ecc=nx.eccentricity(G, source),
                   diam=nx.diameter(G),
                   bound=spectral_bound(lam2, lamN, N))
    else:
        row.update(components=nx.number_connected_components(G),
                   reach=len(nx.node_connected_component(G, source)) / N)
    rows.append(row)

# ===== 3. 打印结果 =====
print(f"{'seed':>4} {'连通':>4} {'最大度':>5} {'λ2':>7} {'λN':>6} "
      f"{'实测轮数':>6} {'离心率':>5} {'直径':>4} {'谱上界':>6}")
for r in rows:
    if r["connected"]:
        print(f"{r['seed']:>4} {'是':>4} {r['max_deg']:>5} {r['lam2']:>7.4f} {r['lamN']:>6.2f} "
              f"{r['flood']:>8} {r['ecc']:>6} {r['diam']:>5} {r['bound']:>7}")
    else:
        print(f"{r['seed']:>4} {'否':>4} {r['max_deg']:>5}   不连通，{r['components']} 个分量，指令只能传到 {round(r['reach']*N)}/{N} 个机器人")

ok = [r for r in rows if r["connected"]]
holds = all(r["flood"] == r["ecc"] <= r["diam"] <= r["bound"] for r in ok)
print(f"\n连通的图：{len(ok)}/{len(rows)}")
if ok:
    print("不等式 实测 = 离心率 ≤ 直径 ≤ 谱上界 是否全部成立：", holds)
else:
    print("没有连通的图，无法检查不等式（随机撒点不保证连通——这正是需要 λ2 下限的原因）")

# ===== 4. 画图：每个连通种子的 实测轮数 / 直径 / 谱上界 =====
if ok:
    plt.rcParams["font.size"] = 22
    fig, ax = plt.subplots(figsize=(11, 6.5))
    x = range(len(ok))
    ax.plot(x, [r["bound"] for r in ok], "s", color="gray", label="spectral bound B(G)", ms=9)
    ax.plot(x, [r["diam"] for r in ok], "o", color="black", mfc="none", mew=2, label="diameter", ms=9)
    ax.plot(x, [r["flood"] for r in ok], "^", color="#1f4e79", label="measured flag rounds", ms=9)
    ax.set_xticks(list(x), [str(r["seed"]) for r in ok], fontsize=18)
    ax.set_xlabel("random seed (connected graphs only)")
    ax.set_ylabel("rounds")
    ax.set_ylim(0, 60)
    ax.legend(fontsize=18)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig("sanity_check_200.png", dpi=300)
    fig.savefig("sanity_check_200.pdf")
    print("已保存 sanity_check_200.png / .pdf")
