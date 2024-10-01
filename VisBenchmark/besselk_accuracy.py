import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

for nu in np.arange(0, 64 + 1) / 4:
    plt.clf()

    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(1, 1, 1)

    ax.set_xlabel("Re z")
    ax.set_ylabel("Im z")

    ax.set_xlim([0, 42])
    ax.set_ylim([0, 42])

    if nu == int(nu):
        nu = int(nu)

    filename = f"besselk_nu{nu}_accuracy"

    data = pd.read_csv(f"../results_ddouble/{filename}.csv")

    r, i, relerr = data["r"].to_numpy(), data["i"].to_numpy(), data["relerr"].to_numpy()

    r = np.concatenate([[0], r]).reshape((85, 85))
    i = np.concatenate([[0], i]).reshape((85, 85))
    relerr = np.concatenate([[1], relerr]).reshape((85, 85))

    with np.errstate(divide='ignore'):
        relerr = np.clip(relerr, 1e-80, 1)
        logrelerr = np.log10(relerr)

    ax.contourf(r, i, logrelerr, levels=[-80, -29.5, -29, -28.5, -28, -27.5, -27, -26.5, -26, 0])
    
    plt.savefig(f"../results_ddouble/besselk_nu{nu}_accuracy.png")