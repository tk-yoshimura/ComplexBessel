import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

methods = ["hankel", "millerbackward", "powerseries"]
colors = ["red", "green", "blue"]

for nu in np.arange(-64, 64 + 1) / 4:
    plt.clf()

    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(1, 1, 1)

    ax.set_xlabel("Re z")
    ax.set_ylabel("Im z")

    ax.set_xlim([0, 64])
    ax.set_ylim([0, 64])

    if nu == int(nu):
        nu = int(nu)

    for method, color in zip(methods, colors):
        filename = f"besselj_nu{nu}_{method}_convergence"

        data = pd.read_csv(f"../results/{filename}.csv")

        r, i, relerr = data["r"].to_numpy(), data["i"].to_numpy(), data["relerr"].to_numpy()

        r = np.concatenate([[0], r]).reshape((129, 129))
        i = np.concatenate([[0], i]).reshape((129, 129))
        relerr = np.concatenate([[1], relerr]).reshape((129, 129))

        logrelerr = np.log10(relerr)

        ax.contour(r, i, logrelerr, levels=[-36.5], colors=color, linestyles="-")
        ax.contour(r, i, logrelerr, levels=[-31.5], colors=color, linestyles="dashed")

    plt.legend([Line2D([], [], color=c) for c in colors], methods)
    plt.savefig(f"../results/besselj_nu{nu}_convergence.png")