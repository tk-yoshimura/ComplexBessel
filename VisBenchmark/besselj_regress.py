import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

methods = ["powerseries"]
colors = ["blue"]

for nu in np.arange(-64, 64 + 1) / 4:
    plt.clf()

    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(1, 1, 1)

    ax.set_xlabel("Re z")
    ax.set_ylabel("Im z")

    ax.set_xlim([0, 42])
    ax.set_ylim([0, 42])

    if nu == int(nu):
        nu = int(nu)

    for method, color in zip(methods, colors):
        filename = f"besselj_nu{nu}_{method}_convergence"

        data = pd.read_csv(f"../results_ddouble/{filename}.csv")

        r, i, relerr = data["r"].to_numpy(), data["i"].to_numpy(), data["relerr"].to_numpy()

        r = np.concatenate([[0], r]).reshape((85, 85))
        i = np.concatenate([[0], i]).reshape((85, 85))
        relerr = np.concatenate([[1], relerr]).reshape((85, 85))

        logrelerr = np.log10(relerr)

        ax.contour(r, i, logrelerr, levels=[-29.5], colors=color, linestyles="-")

    i = np.arange(0, 43)
    r = 7.5 + np.abs(nu) * (3.57e-1 + np.abs(nu) * 5.23e-3) + i * (4.67e-1 + np.abs(nu) * -1.51e-2)

    ax.plot(r, i, color="green")

    plt.legend([Line2D([], [], color=c) for c in colors], methods)
    plt.savefig(f"../results_ddouble/besselj_nu{nu}_convergence_regress.png")