# 环境自检：运行后应在同目录生成 env_check.png
import matplotlib, numpy, scipy, networkx, qrcode
import matplotlib.pyplot as plt
print("matplotlib", matplotlib.__version__, "| numpy", numpy.__version__,
      "| scipy", scipy.__version__, "| networkx", networkx.__version__)
fig, ax = plt.subplots(figsize=(4, 3))
ax.plot([0, 1], [0, 1])
fig.savefig(__file__.replace(".py", ".png"), dpi=150)
print("OK")
