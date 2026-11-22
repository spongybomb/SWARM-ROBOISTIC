# 画海报第 5 板块的曲线图：连通性下限 λ_min  →  最坏情况传播轮数
# 运行方法：在命令行输入  D:\python\python.exe plot_design_curve.py
# 运行后会在同一个文件夹生成 design_curve.png 和 design_curve.pdf

import math
import matplotlib.pyplot as plt

# ===== 1. 参数（来自你的 Extended Abstract §4）=====
N = 200          # 机器人数量
LAMBDA_N = 16    # λ_N 的上界 = 2 × d_max = 2 × 8
REQUIREMENT = 20 # 监管要求：指令必须在 20 轮内传到所有机器人


# ===== 2. Chung–Faber–Manteuffel 直径上界公式 =====
def bound(lam2):
    top = math.acosh(N - 1)
    bottom = math.acosh((LAMBDA_N + lam2) / (LAMBDA_N - lam2))
    return math.floor(top / bottom) + 1


# ===== 3. 计算曲线上的点：λ_min 从 0.10 到 2.00 =====
xs = [0.10 + 0.01 * k for k in range(191)]
ys = [bound(x) for x in xs]

# 找到满足 20 轮要求所需的最小 λ_min
need = next(x for x, y in zip(xs, ys) if y <= REQUIREMENT)
print("满足", REQUIREMENT, "轮要求的最小 λ_min ≈", round(need, 2))


# ===== 4. 画图 =====
plt.rcParams["font.size"] = 24        # 海报上看得清，至少 24
fig, ax = plt.subplots(figsize=(10, 7))

ax.step(xs, ys, where="post", color="black", linewidth=3)
ax.axhline(REQUIREMENT, color="gray", linestyle="--", linewidth=2)
ax.axvline(need, color="gray", linestyle=":", linewidth=2)
ax.text(1.2, REQUIREMENT + 1, f"requirement: {REQUIREMENT} rounds", color="gray")

ax.set_xlabel(r"connectivity floor $\lambda_{\min}$")
ax.set_ylabel("worst-case rounds")
ax.set_xlim(0.1, 2.0)
ax.grid(alpha=0.3)

fig.tight_layout()
fig.savefig("design_curve.png", dpi=300)
fig.savefig("design_curve.pdf")
print("已保存 design_curve.png / design_curve.pdf")
