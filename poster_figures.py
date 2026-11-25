# 海报用图（A1 竖版，按实际打印尺寸出图，风格参考 ADMA 样板：衬线字体、深蓝/青绿/锈红、直接标数值）
# 先运行 static_experiments.py，再运行：D:\python\python.exe poster_figures.py
# 输出：results/poster/*.png（300 DPI）和 *.pdf（矢量）
import json
import math
import os
import numpy as np
import matplotlib.pyplot as plt
import sanity_check_200 as s

NAVY, TEAL, SLATE, RUST, GREY = "#142C47", "#1F7A6D", "#5B7083", "#B5542B", "#9AA5B1"
MM = 1 / 25.4
plt.rcParams.update({
    "font.family": "serif", "font.serif": ["Times New Roman", "DejaVu Serif"],
    "mathtext.fontset": "stix", "font.size": 17, "axes.edgecolor": SLATE,
    "axes.labelcolor": NAVY, "xtick.color": SLATE, "ytick.color": SLATE,
    "axes.spines.top": False, "axes.spines.right": False, "savefig.transparent": True,
})
os.makedirs("results/poster", exist_ok=True)
d = json.load(open("results/static_experiments.json", encoding="utf-8"))
N = d["N"]


def save(fig, name):
    fig.savefig(f"results/poster/{name}.png", dpi=300)
    fig.savefig(f"results/poster/{name}.pdf")
    plt.close(fig)
    print("已保存", name)


# ---------- F1 通信半径扫描（半栏 128 × 100 mm）----------
e1 = d["E1_radius_sweep"]
radii = sorted({r["R"] for r in e1 if r["R"] <= 4.0})
conn = [sum(r["connected"] for r in e1 if r["R"] == R) for R in radii]
reach = [100 * np.mean([r["reach"] for r in e1 if r["R"] == R]) for R in radii]
fig, ax = plt.subplots(figsize=(128 * MM, 100 * MM))
ax.axvspan(1.76, 2.24, color=RUST, alpha=0.10, lw=0)
ax.text(2.0, 3, "paper\nscenario", ha="center", va="bottom", fontsize=12, color=RUST)
ax.bar(radii, [100 * c / 20 for c in conn], width=0.32, color=NAVY, label="networks connected")
ax.plot(radii, reach, "o-", color=TEAL, lw=2.5, ms=8, label="robots reached (mean)")
ax.text(1.92, reach[radii.index(2.0)], f"{reach[radii.index(2.0)]:.0f}%", ha="right", va="center",
        color=TEAL, fontsize=15)
# 横轴标签两行：半径 + 连通网络数，避免数字和折线重叠
ax.set_xticks(radii, [f"{R:g}\n{c}/20" for R, c in zip(radii, conn)], fontsize=14)
ax.set_xlabel("radio range R (m)  /  networks connected")
ax.set_ylabel("%")
ax.set_yticks(range(0, 101, 20))
ax.set_ylim(0, 132)
ax.legend(loc="upper left", fontsize=13, frameon=False, ncol=1)
fig.tight_layout()
save(fig, "F1_radius_sweep")

# ---------- F2 上界是否成立（整栏 262 × 100 mm）----------
uniq, seen = [], set()
for r in e1:
    if r["connected"] and r["R"] <= 3.5:
        key = (r["seed"], r["flood"], r["diam"], r["bound"], round(r["lam2"], 9))
        if key not in seen:
            seen.add(key)
            uniq.append(r)
uniq.sort(key=lambda r: r["lam2"])
holds = sum(r["flood"] == r["ecc"] <= r["diam"] <= r["bound"] for r in uniq)
ratio = [r["bound"] / r["flood"] for r in uniq]
fig, ax = plt.subplots(figsize=(262 * MM, 115 * MM))
x = np.arange(len(uniq))
ax.vlines(x, [r["flood"] for r in uniq], [r["bound"] for r in uniq], color=GREY, lw=1.2, alpha=0.6)
ax.plot(x, [r["bound"] for r in uniq], "s", color=SLATE, ms=7, label="spectral bound B(G)")
ax.plot(x - 0.18, [r["diam"] for r in uniq], "o", mfc="none", mec=NAVY, mew=1.8, ms=8, label="diameter")
ax.plot(x + 0.18, [r["flood"] for r in uniq], "^", color=TEAL, ms=8, label="measured flag rounds")
ax.set_xticks([])
ax.set_xlabel(r"connected random networks, sorted by $\lambda_2$ (low $\rightarrow$ high)")
ax.set_ylabel("rounds")
ax.set_ylim(0, 118)
ax.legend(loc="upper right", ncol=3, fontsize=14, frameon=False)
ax.text(0.01, 0.84, f"inequality held in {holds}/{len(uniq)} networks;  bound / measured = "
        f"{min(ratio):.1f}–{max(ratio):.1f}×", transform=ax.transAxes, va="top", fontsize=15, color=NAVY)
fig.tight_layout()
save(fig, "F2_bound_check")
print(f"F2：{holds}/{len(uniq)} 成立；倍数 {min(ratio):.2f}–{max(ratio):.2f}，中位数 {np.median(ratio):.2f}")

# ---------- F3 flag vs consensus（整栏 262 × 110 mm，两幅）----------
e2 = d["E2_flag_vs_consensus"]
l2 = np.array([r["lam2"] for r in e2])
fl = np.array([r["flag"] for r in e2])
co = np.array([r["consensus"] for r in e2], float)
lg = np.array([r["lam_grounded"] for r in e2])
pred = math.log(1 / (1 - d["threshold"])) / (d["beta"] * lg)
fig, (a1, a2) = plt.subplots(1, 2, figsize=(262 * MM, 125 * MM), gridspec_kw={"width_ratios": [1.25, 1], "wspace": 0.42})
a1.semilogy(l2, co, "o", color=RUST, ms=7, alpha=0.85, label="consensus (all robots ≥ 0.95)")
a1.semilogy(l2, fl, "^", color=TEAL, ms=7, alpha=0.9, label="flag forwarding")
a1.set_xlabel(r"$\lambda_2$ of the network")
a1.set_ylabel("rounds until every robot complies")
a1.set_ylim(8, 3e4)
a1.legend(loc="center right", fontsize=13, frameon=False)
a1.text(0.03, 0.97, f"median {np.median(fl):.0f} vs {np.median(co):,.0f} rounds  (≈{np.median(co / fl):.0f}×)",
        transform=a1.transAxes, va="top", fontsize=14, color=NAVY)
a1.set_title("A   Same networks, two protocols", loc="left", fontsize=16, color=NAVY, fontweight="bold")
lim = [2500, 12000]
a2.plot(lim, lim, "--", color=GREY, lw=1.5)
a2.plot(pred, co, "o", color=RUST, ms=7, alpha=0.85)
a2.set_xlim(lim); a2.set_ylim(lim)
a2.set_xlabel(r"predicted  $\ln 20\,/\,(\beta\,\lambda(L_g))$")
a2.set_ylabel("measured consensus rounds")
a2.text(0.05, 0.95, f"ratio {np.min(co / pred):.2f}–{np.max(co / pred):.2f}", transform=a2.transAxes,
        va="top", fontsize=14, color=NAVY)
a2.set_title(r"B   Grounded-Laplacian prediction", loc="left", fontsize=16, color=NAVY, fontweight="bold")
fig.subplots_adjust(left=0.08, right=0.98, bottom=0.15, top=0.9, wspace=0.42)
save(fig, "F3_flag_vs_consensus")

# ---------- F4 节点损失（整栏 262 × 90 mm）----------
e3 = d["E3_node_loss"]
modes = [("none", "no loss"), ("random", "random 15%"), ("region", "regional 15%\n(NE collapse)")]
fig, (b1, b2) = plt.subplots(1, 2, figsize=(262 * MM, 105 * MM))
cols = [NAVY, SLATE, RUST]
tot = sum(r["mode"] == "none" for r in e3)
for k, (m, lab) in enumerate(modes):
    rows = [r for r in e3 if r["mode"] == m]
    c = sum(r["connected"] for r in rows)
    b1.bar(k, c, color=cols[k], width=0.6)
    b1.text(k, c + 0.3, f"{c}/{tot}", ha="center", va="bottom", fontsize=15, color=cols[k])
    ok = [r["lam2"] for r in rows if r["connected"]]
    b2.bar(k, np.median(ok), color=cols[k], width=0.6)
    b2.text(k, np.median(ok) + 0.001, f"{np.median(ok):.3f}", ha="center", va="bottom", fontsize=15, color=cols[k])
for b, title, ylab in ((b1, "A   Still connected after loss", "networks"),
                       (b2, r"B   $\lambda_2$ of those still connected (median)", r"$\lambda_2$")):
    b.set_xticks(range(3), [lab for _, lab in modes], fontsize=14)
    b.set_title(title, loc="left", fontsize=16, color=NAVY, fontweight="bold")
    b.set_ylabel(ylab)
b1.set_ylim(0, tot + 3)
b2.set_ylim(0, 0.068)
fig.tight_layout()
save(fig, "F4_node_loss")

# ---------- F5 设计曲线（半栏 128 × 100 mm）----------
LN, REQ = 16, 20
xs = [0.10 + 0.005 * k for k in range(381)]
ys = [s.spectral_bound(x, LN, N) for x in xs]
fig, ax = plt.subplots(figsize=(128 * MM, 100 * MM))
ax.step(xs, ys, where="post", color=NAVY, lw=2.5)
ax.axhline(REQ, color=RUST, ls="--", lw=1.8)
ax.axvline(0.36, color=RUST, ls=":", lw=1.8)
ax.text(1.95, REQ + 1, "requirement: 20 rounds", ha="right", color=RUST, fontsize=14)
ax.text(0.40, 37, r"floor $\geq$ 0.36", color=RUST, fontsize=14)
ax.set_xlabel(r"connectivity floor $\lambda_{\min}$")
ax.set_ylabel("worst-case rounds (bound)")
ax.set_xlim(0.1, 2.0); ax.set_ylim(5, 42)
ax.set_title(r"N = 200, $\lambda_N \leq 16$ (analytic)", loc="right", fontsize=13, color=SLATE)
fig.tight_layout()
save(fig, "F5_design_curve")
