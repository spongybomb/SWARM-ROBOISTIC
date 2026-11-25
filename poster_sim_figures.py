# 运动模拟的海报图（先运行 cbf_sim.py 生成 results/cbf_sim.json）
#   F6_sim_summary   三种环境 × 三种设置：指令到达率、已搜索面积、指令传遍所需轮数
#   F7_sim_snapshot  空间快照（种子 0）：上排不加下限，下排加 λ_min = 0.05；青绿 = 收到撤离指令
#   F8_code          cbf_filter 源码截图（来自 cbf_sim.py 原文）
# 运行：D:\python\python.exe poster_sim_figures.py
import inspect
import json
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import cbf_sim as c

NAVY, TEAL, SLATE, RUST, GREY = "#142C47", "#1F7A6D", "#5B7083", "#B5542B", "#9AA5B1"
MM = 1 / 25.4
plt.rcParams.update({
    "font.family": "serif", "font.serif": ["Times New Roman", "DejaVu Serif"],
    "mathtext.fontset": "stix", "font.size": 16, "axes.edgecolor": SLATE,
    "axes.labelcolor": NAVY, "xtick.color": SLATE, "ytick.color": SLATE,
    "axes.spines.top": False, "axes.spines.right": False, "savefig.transparent": True,
})
os.makedirs("results/poster", exist_ok=True)
ENV_LABEL = {"open": "Collapsed-building void", "corridor": "Collapsed tunnel", "rooms": "Separated survivor pockets"}
ENV_SHORT = {"open": "Building\nvoid", "corridor": "Collapsed\ntunnel", "rooms": "Survivor\npockets"}
SET_LABEL = {None: "no floor", 0.05: r"$\lambda_{\min}=0.05$", 0.2: r"$\lambda_{\min}=0.2$"}
SET_COLOR = {None: GREY, 0.05: TEAL, 0.2: NAVY}


def save(fig, name):
    fig.savefig(f"results/poster/{name}.png", dpi=300)
    fig.savefig(f"results/poster/{name}.pdf")
    plt.close(fig)
    print("已保存", name)


# ---------------- F6 汇总 ----------------
d = json.load(open("results/cbf_sim.json", encoding="utf-8"))
runs, envs, floors = d["runs"], list(d["params"]["envs"]), d["params"]["floors"]
summary = {}
for e in envs:
    for f in floors:
        rs = [r for r in runs if r["env"] == e and r["lam_min"] == f]
        done = [r["override_rounds"] for r in rs if r["override_rounds"] is not None]
        summary[(e, f)] = {
            "n": len(rs),
            "received": 100 * np.mean([r["override_received"] for r in rs]),
            "all_received": sum(r["override_rounds"] is not None for r in rs),
            "coverage": 100 * np.mean([r["coverage_end"] for r in rs]),
            "rounds_med": float(np.median(done)) if done else None,
            "bound_med": float(np.median([r["bound_end"] for r in rs if r["bound_end"]])) if f is not None else None,
            "viol": sum(r["violation_steps"] for r in rs),
            "backtrack_frac": float(np.mean([r["backtrack_steps"] for r in rs])) / (c.T_END / c.DT),
        }
json.dump({f"{e}|{f}": v for (e, f), v in summary.items()}, open("results/cbf_sim_summary.json", "w"), indent=1)

fig, axes = plt.subplots(1, 3, figsize=(262 * MM, 92 * MM))
xg = np.arange(len(envs))
wbar = 0.26
for k, f in enumerate(floors):
    xs = xg + (k - 1) * wbar
    rec = [summary[(e, f)]["received"] for e in envs]
    cov = [summary[(e, f)]["coverage"] for e in envs]
    axes[0].bar(xs, rec, wbar, color=SET_COLOR[f], label=SET_LABEL[f])
    axes[1].bar(xs, cov, wbar, color=SET_COLOR[f])
    for x, e in zip(xs, envs):
        s = summary[(e, f)]
        axes[0].text(x, 2, f"{s['all_received']}/{s['n']}", ha="center", va="bottom", fontsize=10,
                     color=NAVY if f is None else "white", rotation=90)
        if s["rounds_med"] is not None and f is not None:
            axes[2].bar(x, s["rounds_med"], wbar, color=SET_COLOR[f])
            axes[2].plot([x - wbar / 2, x + wbar / 2], [s["bound_med"]] * 2, color=RUST, lw=2.2)
for ax, title, ylab in ((axes[0], "A   Robots that got the order", "% of robots"),
                        (axes[1], "B   Area searched by 60 s", "% of cells"),
                        (axes[2], "C   Rounds to reach all", "rounds (median)")):
    ax.set_xticks(xg, [ENV_SHORT[e] for e in envs], fontsize=12)
    ax.set_title(title, loc="left", fontsize=15, color=NAVY, fontweight="bold")
    ax.set_ylabel(ylab)
axes[0].set_ylim(0, 105); axes[1].set_ylim(0, 105)
axes[0].plot([], [], color=RUST, lw=2.2, label="spectral bound (C, median)")
axes[0].legend(loc="lower left", bbox_to_anchor=(0, 1.13), ncol=4, fontsize=13, frameon=False)
axes[2].set_ylim(0, 135)
axes[2].text(0.99, 0.99, "no floor: order rarely" + "\n" + "reached all robots", transform=axes[2].transAxes, ha="right", va="top", fontsize=11, color=SLATE)
fig.subplots_adjust(left=0.07, right=0.99, bottom=0.14, top=0.78, wspace=0.38)
save(fig, "F6_sim_summary")
for (e, f), s in summary.items():
    print(f"{e:8s} {str(f):5s}: 收到 {s['received']:.0f}%（全部收到 {s['all_received']}/{s['n']}），搜索 {s['coverage']:.0f}%，"
          f"轮数中位数 {s['rounds_med']}，上界中位数 {s['bound_med']}，越界步数 {s['viol']}，回退步比例 {s['backtrack_frac']:.0%}")

# ---------------- F7 空间快照 ----------------
fig = plt.figure(figsize=(262 * MM, 104 * MM))
gs = fig.add_gridspec(2, 4, height_ratios=[3.05, 1.0], hspace=0.42, wspace=0.08,
                      left=0.01, right=0.99, top=0.89, bottom=0.09)
slots = {("open", None): gs[0, 0], ("open", 0.05): gs[0, 1], ("rooms", None): gs[0, 2], ("rooms", 0.05): gs[0, 3],
         ("corridor", None): gs[1, 0:2], ("corridor", 0.05): gs[1, 2:4]}
for (e, f), slot in slots.items():
    r = c.run(e, 0, f, keep_positions=True)
    X, rec = np.array(r["X"]), np.array(r["received"])
    ax = fig.add_subplot(slot)
    segs = [[X[i], X[j]] for i, j in r["edges"]]
    ax.add_collection(LineCollection(segs, colors=GREY, linewidths=0.3, alpha=0.6))
    ax.scatter(*np.array(r["targets"]).T, s=3, color=RUST, alpha=0.18, lw=0)
    ax.scatter(*X[~rec].T, s=8, color=GREY, lw=0)
    ax.scatter(*X[rec].T, s=8, color=TEAL, lw=0)
    ax.plot(0, 0, marker="s", color=RUST, ms=7)
    W, H = c.ENVS[e]["size"]
    ax.set_xlim(-1, W + 1); ax.set_ylim(-1, H + 1); ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(True); sp.set_color("#D5DCE3")
    ax.set_title(f"{ENV_LABEL[e]}" + "\n" + f"{'no floor' if f is None else 'floor 0.05'}: {rec.mean() * 100:.0f}% reached",
                 fontsize=10.5, color=NAVY, loc="left", linespacing=1.1)
fig.text(0.01, 0.015, "Seed 0 at t = 60 s.  Teal: received the evacuation order;  grey: never received;  "
         "red square: commander's terminal;  faint red: search targets.", fontsize=10.5, color=SLATE)
save(fig, "F7_sim_snapshot")

# ---------------- F8 源码截图（cbf_filter 节选，逐行取自 cbf_sim.py，不改写）----------------
src = inspect.getsource(c.cbf_filter).split("\n")
keep = ("def cbf_filter", "u0 = np.clip", "active = phi0", "phi0 = (g", "feasible = UMAX", "for _ in range(50)",
        "mid = 0.5", "ok = (gi * np.clip(ui + mid", "hi, lo = np.where", "u[idx] = np.clip", "return u, active")
lines = []
for l in src:
    if any(l.strip().startswith(k) for k in keep):
        if l.strip().startswith("for _ in range(50)") and lines:
            lines.append("    ...                                                   # bracket search omitted")
        lines.append(l.rstrip())
code = "\n".join(lines)
fig = plt.figure(figsize=(262 * MM, (len(lines) + 2.2) * 5.0 * MM))
fig.patch.set_facecolor("#F4F6F9")
fig.text(0.015, 0.96, code, family="Consolas", fontsize=12.5, va="top", color=NAVY, linespacing=1.25)
fig.savefig("results/poster/F8_code.png", dpi=300, facecolor=fig.get_facecolor())
plt.close(fig)
print("已保存 F8_code（", len(lines), "行 ）")
print(code)