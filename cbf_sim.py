# 运动模拟：名义搜索控制器  vs  名义控制器 + 连通性 CBF（λ2 下限）
# 对应论文 §4 的机制与 §5 的评估方案（缩小版）。
#
# 模型（与论文一致的部分）
#   单积分器 ẋ_i = u_i；通信半径 R = 2 m；平滑边权 a_ij(d)：d < r0 时为 1，r0 ≤ d < R 余弦衰减到 0
#   h(x) = λ2(L(x)) − λ_min；每个机器人解一个局部 QP：
#       min ½‖u_i − u_i^nom‖²   s.t.  g_i·u_i ≥ b_i,  |u_i,x|, |u_i,y| ≤ u_max
#       b_i = −γ h ‖g_i‖² / Σ_j ‖g_j‖²（Σ b_i = −γ h；分配需要全局量，与集中式 λ2 一样属参考实现）
#   其中 g_i = ∂λ2/∂x_i = Σ_j ∂a_ij/∂x_i (v_i − v_j)²（v = 归一化 Fiedler 向量）
#   指令以 flag 方式逐轮传播，每轮一跳（A1），每个控制周期一轮通信
#
# 与论文不同 / 简化的部分（海报上要写明）
#   - λ2 和 Fiedler 向量是集中式精确计算的（参考实现，不是去中心化估计，A3 视为成立）
#   - 没有 8 条链路上限（上限会让边权不光滑）；谱上界因此用"支撑图"的 λN 计算（见 bound_now）
#   - 不做避碰保证，只在名义控制器里加了近距离排斥
#
# 运行：D:\python\python.exe cbf_sim.py   → results/cbf_sim.json（约 10–20 分钟）
import json
import math
import time
import numpy as np

N, R = 200, 2.0
R0 = 0.6 * R
UMAX, DT, T_END = 0.5, 0.1, 60.0          # 单轴速度上限 0.5 m/s，控制周期 0.1 s，60 s 后发撤离指令
GAMMA = 1.0
MARGIN = 0.15                              # 实际按 (1 + 15%)·λ_min 执行，吸收离散步长造成的越界
K_NOM, REPULSE = 0.8, 0.35                # 名义控制器：奔向目标点的增益；排斥距离 (m)
MAX_ROUNDS = 400                          # 指令传播最多看 400 轮（= 40 s）

# 三种环境：目标点（被搜索的位置）分布不同；所有机器人都从入口 (0..6, 0..6) 的一团出发
ENVS = {
    "open":     {"size": (20.0, 20.0), "desc": "open 20 × 20 m"},
    "corridor": {"size": (60.0, 6.67), "desc": "corridor 60 × 6.7 m"},
    "rooms":    {"size": (30.0, 30.0), "desc": "three rooms in 30 × 30 m",
                 "rooms": [(25.0, 5.0), (5.0, 25.0), (25.0, 25.0)], "room_r": 3.5},
}


# ---------------- 图与谱 ----------------
def weights(X):
    """平滑边权矩阵 A 与距离矩阵 D"""
    diff = X[:, None, :] - X[None, :, :]
    D = np.sqrt((diff ** 2).sum(-1))
    A = np.where(D < R0, 1.0, 0.5 * (1 + np.cos(np.pi * (D - R0) / (R - R0))))
    A[D >= R] = 0.0
    np.fill_diagonal(A, 0.0)
    return A, D, diff


def spectrum(A):
    L = np.diag(A.sum(1)) - A
    w, V = np.linalg.eigh(L)
    return w, V


def lambda2_grad(X, A, D, diff, v):
    """g_i = Σ_j ∂a_ij/∂x_i (v_i − v_j)²；过渡区以外梯度为 0"""
    mid = (D >= R0) & (D < R)
    dadd = np.zeros_like(D)
    dadd[mid] = -0.5 * np.pi / (R - R0) * np.sin(np.pi * (D[mid] - R0) / (R - R0))
    with np.errstate(invalid="ignore", divide="ignore"):
        unit = np.where(D[..., None] > 1e-9, diff / D[..., None], 0.0)
    coef = dadd * (v[:, None] - v[None, :]) ** 2
    return (coef[..., None] * unit).sum(1)


# ---------------- 局部 QP（盒约束 + 一个半空间），精确解 ----------------
def cbf_filter(u_nom, g, b):
    """对每个机器人：min ½‖u − u_nom‖²  s.t.  g·u ≥ b,  |u| ≤ UMAX（逐轴）
       KKT ⇒ u(η) = clip(u_nom + η g)，η ≥ 0；g·u(η) 随 η 单调不减 → 对 η 二分。
       返回 u、是否激活、是否不可行。"""
    u0 = np.clip(u_nom, -UMAX, UMAX)
    phi0 = (g * u0).sum(1)
    active = phi0 < b
    feasible = UMAX * np.abs(g).sum(1) >= b
    u = u0.copy()
    idx = np.where(active & feasible)[0]
    if len(idx):
        gi, ui, bi = g[idx], u_nom[idx], b[idx]
        lo = np.zeros(len(idx))
        hi = np.ones(len(idx))
        for _ in range(60):                                   # grow the upper bound until feasible
            ok = (gi * np.clip(ui + hi[:, None] * gi, -UMAX, UMAX)).sum(1) >= bi
            if ok.all():
                break
            hi = np.where(ok, hi, hi * 4)
        for _ in range(50):                                   # bisection on eta
            mid = 0.5 * (lo + hi)
            ok = (gi * np.clip(ui + mid[:, None] * gi, -UMAX, UMAX)).sum(1) >= bi
            hi, lo = np.where(ok, mid, hi), np.where(ok, lo, mid)
        u[idx] = np.clip(ui + hi[:, None] * gi, -UMAX, UMAX)
    bad = np.where(active & ~feasible)[0]                     # infeasible: best effort, logged
    u[bad] = UMAX * np.sign(g[bad])
    return u, active, active & ~feasible


# ---------------- 环境与名义控制器 ----------------
def make_targets(env, rng):
    cfg = ENVS[env]
    W, H = cfg["size"]
    if env == "rooms":
        c = np.array(cfg["rooms"])[rng.integers(0, 3, N)]
        ang, rad = rng.uniform(0, 2 * np.pi, N), cfg["room_r"] * np.sqrt(rng.uniform(0, 1, N))
        return c + np.c_[rad * np.cos(ang), rad * np.sin(ang)]
    return np.c_[rng.uniform(0, W, N), rng.uniform(0, H, N)]


def initial_positions(rng):
    """入口处 6 × 6 m 的一团（保证连通且 λ2 足够大）"""
    return rng.uniform(0, 6, (N, 2))


def nominal(X, targets, D, diff):
    u = K_NOM * (targets - X)
    close = (D < REPULSE) & (D > 1e-9)
    with np.errstate(invalid="ignore", divide="ignore"):
        rep = np.where(close[..., None], diff / D[..., None] ** 2, 0.0).sum(1)
    return u + 0.05 * rep


# ---------------- 指标 ----------------
def components_from(A, src):
    seen, frontier = {src}, [src]
    nb = [np.nonzero(A[i] > 0)[0] for i in range(len(A))]
    while frontier:
        nxt = []
        for i in frontier:
            for j in nb[i]:
                if j not in seen:
                    seen.add(j); nxt.append(j)
        frontier = nxt
    return seen


def bound_now(A, lam2):
    """谱直径上界 B(λ2_weighted, λN_support)。
       边权 ≤ 1 ⇒ L_w ≼ L_support ⇒ λ2_w ≤ λ2_s；用 λ2_w（更小）与支撑图的 λN 仍是 diam 的有效上界。"""
    S = (A > 0).astype(float)
    lamN = np.linalg.eigvalsh(np.diag(S.sum(1)) - S)[-1]
    if lam2 <= 1e-12:
        return None, float(lamN)
    return math.floor(math.acosh(N - 1) / math.acosh((lamN + lam2) / (lamN - lam2))) + 1, float(lamN)


def coverage_grid(env):
    W, H = ENVS[env]["size"]
    gx, gy = np.meshgrid(np.arange(0.5, W, 1.0), np.arange(0.5, H, 1.0))
    return np.c_[gx.ravel(), gy.ravel()]


# ---------------- 一次运行 ----------------
def run(env, seed, lam_min, keep_positions=False):
    """lam_min = None → 只有名义控制器；否则加 CBF"""
    rng = np.random.default_rng(seed)
    X = initial_positions(rng)
    targets = make_targets(env, rng)
    cells = coverage_grid(env)
    searched = np.zeros(len(cells), bool)
    src = int(np.argmin(np.linalg.norm(X, axis=1)))          # 离操作终端 (0,0) 最近的机器人做注入点
    log = {"t": [], "lam2": [], "reach": [], "coverage": [], "active": []}
    viol, infeas, steps_active, corr = 0, 0, 0, 0.0
    backtracks, backtrack_steps, min_gap = 0, 0, float("inf")
    A, D, diff = weights(X)
    w, V = spectrum(A)
    h0 = w[1] - (lam_min or 0)
    steps = int(T_END / DT)
    rounds, received = None, None
    for k in range(steps + MAX_ROUNDS):
        t = k * DT
        u_nom = nominal(X, targets, D, diff)
        if lam_min is not None:
            lam2, v = w[1], V[:, 1]
            g = lambda2_grad(X, A, D, diff, v)
            # 约束总量 Σ g_i·u_i ≥ −γ h 按 ‖g_i‖² 分给各机器人（梯度为 0 的机器人分到 0，不会无解）
            h = lam2 - (1 + MARGIN) * lam_min
            gn2 = (g ** 2).sum(1)
            b = -GAMMA * h * gn2 / max(gn2.sum(), 1e-12)
            u, act, bad = cbf_filter(u_nom, g, b)
            infeas += int(bad.sum()); steps_active += int(act.sum())
            corr += float(np.linalg.norm(u - np.clip(u_nom, -UMAX, UMAX), axis=1).sum())
        else:
            u = np.clip(u_nom, -UMAX, UMAX)
            act = np.zeros(N, bool)
        # 离散时间安全检查：若这一步会让 λ2 跌破 λ_min，就把步长减半重试（最多 8 次，最后退回 u = 0）
        X_new = X + DT * u
        A_new, D_new, diff_new = weights(X_new)
        w_new, V_new = spectrum(A_new)
        if lam_min is not None:
            tries = 0
            while w_new[1] < lam_min and tries < 8:
                u = 0.5 * u if tries < 7 else 0.0 * u
                X_new = X + DT * u
                A_new, D_new, diff_new = weights(X_new)
                w_new, V_new = spectrum(A_new)
                tries += 1
            backtracks += tries
            backtrack_steps += tries > 0
        X, A, D, diff, w, V = X_new, A_new, D_new, diff_new, w_new, V_new
        if lam_min is not None and w[1] < lam_min - 1e-9:
            viol += 1
        if lam_min is not None:
            min_gap = min(min_gap, float(w[2] - w[1]))
        searched |= (np.linalg.norm(cells[:, None, :] - X[None, :, :], axis=2) < 1.0).any(1)
        if k % 10 == 0 or k == steps:
            log["t"].append(round(t, 2)); log["lam2"].append(float(w[1]))
            log["reach"].append(len(components_from(A, src)) / N)
            log["coverage"].append(float(searched.mean())); log["active"].append(int(act.sum()))
        # ---- 第 60 s：注入撤离指令，之后每个控制周期传播一跳（拓扑在移动） ----
        if k == steps:
            snap_lam2 = float(w[1])
            snap_bound, snap_lamN = bound_now(A, snap_lam2)
            snap_reach = len(components_from(A, src)) / N
            snap_cov = float(searched.mean())
            received = np.zeros(N, bool); received[src] = True
            X_inject = X.copy(); A_inject = A.copy()
        elif k > steps:
            new = received | (A[:, received] > 0).any(1)
            received = new
            if received.all():
                rounds = k - steps
                break
    extra = {}
    if keep_positions:                                         # 给空间快照图用：注入时的位置、边、最终收到指令的机器人
        ii, jj = np.nonzero(np.triu(A_inject) > 0)
        extra = {"X": X_inject.tolist(), "edges": np.c_[ii, jj].tolist(), "src": src,
                 "received": received.tolist(), "targets": targets.tolist()}
    return {**extra, "env": env, "seed": seed, "lam_min": lam_min,
            "lam2_end": snap_lam2, "bound_end": snap_bound, "lamN_support_end": snap_lamN,
            "reach_end": snap_reach, "coverage_end": snap_cov,
            "override_rounds": rounds, "override_received": float(received.mean()),
            "lam2_init": float(h0 + (lam_min or 0)),
            "violation_steps": viol, "backtrack_steps": backtrack_steps, "backtracks": backtracks,
            "min_spectral_gap": None if min_gap == float("inf") else min_gap, "infeasible_agent_steps": infeas,
            "active_agent_steps": steps_active, "correction_sum": corr, "log": log}


def _job(args):
    import os
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    env, seed, lam_min = args
    return run(env, seed, lam_min)


if __name__ == "__main__":
    # 每完成一次就追加写入 results/cbf_sim_runs.jsonl；中断后重跑会跳过已完成的组合
    import os
    import sys
    from multiprocessing import Pool
    seeds = range(int(sys.argv[1]) if len(sys.argv) > 1 else 20)
    workers = int(sys.argv[2]) if len(sys.argv) > 2 else 12
    floors = [None, 0.05, 0.2]
    path = "results/cbf_sim_runs.jsonl"
    done = set()
    if os.path.exists(path):
        for line in open(path, encoding="utf-8"):
            r = json.loads(line)
            done.add((r["env"], r["seed"], r["lam_min"]))
    jobs = [(e, s, f) for e in ENVS for f in floors for s in seeds if (e, s, f) not in done]
    print(f"已完成 {len(done)}，待运行 {len(jobs)}，并行 {workers} 个进程", flush=True)
    t0 = time.time()
    with Pool(workers) as pool, open(path, "a", encoding="utf-8") as fout:
        for k, r in enumerate(pool.imap_unordered(_job, jobs), 1):
            fout.write(json.dumps(r) + "\n"); fout.flush()
            print(f"[{k}/{len(jobs)}] {r['env']:8s} floor={r['lam_min']} seed={r['seed']:2d}  λ2_end={r['lam2_end']:.3f}  "
                  f"recv={r['override_received']:.0%}  cover={r['coverage_end']:.0%}  rounds={r['override_rounds']}  "
                  f"bound={r['bound_end']}  viol={r['violation_steps']}  backtrack={r['backtrack_steps']}  ({time.time() - t0:.0f}s)",
                  flush=True)
    runs = [json.loads(l) for l in open(path, encoding="utf-8")]
    runs.sort(key=lambda r: (list(ENVS).index(r["env"]), str(r["lam_min"]), r["seed"]))
    with open("results/cbf_sim.json", "w", encoding="utf-8") as f:
        json.dump({"params": {"N": N, "R": R, "r0": R0, "umax": UMAX, "dt": DT, "T_end": T_END, "gamma": GAMMA,
                              "margin": MARGIN, "floors": floors, "envs": {k: v["desc"] for k, v in ENVS.items()}},
                   "runs": runs}, f)
    print(f"已保存 results/cbf_sim.json（{len(runs)} 次运行）")
