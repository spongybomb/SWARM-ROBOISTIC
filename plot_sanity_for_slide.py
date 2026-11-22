# 把 sanity_check_200.py 的结果画成适合 16:9 幻灯片的版本
# 运行：D:\python\python.exe plot_sanity_for_slide.py   →  sanity_check_slide.png
import matplotlib.pyplot as plt
import sanity_check_200 as s   # 会先重新跑一遍检查（R = 3.0 m）

ok = [r for r in s.rows if r["connected"]]
plt.rcParams["font.size"] = 15
fig, ax = plt.subplots(figsize=(8, 5.2))
x = range(len(ok))
ax.plot(x, [r["bound"] for r in ok], "s", color="#8193A6", ms=9, label="spectral bound B(G)")
ax.plot(x, [r["diam"] for r in ok], "o", color="#10243A", mfc="none", mew=2, ms=10, label="graph diameter")
ax.plot(x, [r["flood"] for r in ok], "^", color="#13B8A6", ms=10, label="measured flag rounds")
ax.set_xticks(list(x), [str(r["seed"]) for r in ok], fontsize=11)
ax.set_xlabel(f"random graph (seed), connected cases: {len(ok)}/{len(s.rows)}")
ax.set_ylabel("rounds until every robot has the order")
ax.set_ylim(0, 62)
ax.legend(loc="upper right", fontsize=12, ncol=3, frameon=False, bbox_to_anchor=(1, 1.1))
ax.grid(alpha=0.3)
for side in ("top", "right"):
    ax.spines[side].set_visible(False)
fig.tight_layout()
fig.savefig("sanity_check_slide.png", dpi=300)
print("已保存 sanity_check_slide.png")
