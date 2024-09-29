import os, itertools
import numpy as np
import mpmath as mpf
import tqdm

dirpath = "../mpmath_expected/besselk/"

os.makedirs(dirpath, exist_ok=True)

mpf.mp.dps = 64

nus = [0, 1, 2, 3, 4, 6, 8, 12, 16, 0.5, 1.5, 0.25, 15.75, 15.5]
xs = np.arange(0, 128 + 1) / 2

for nu in tqdm.tqdm(nus, "Processing...", total=len(nus)):

    f = open(f"{dirpath}/besselk_nu{nu}.csv", "w")

    f.write("r,i,z\n")

    for r, i in itertools.product(xs, xs):
        if r == 0 and i == 0:
            continue

        c = r + i * 1j

        z = mpf.besselk(nu, c)
        s = str(z).replace("j", "i").replace(" ", "").replace("(", "").replace(")", "")
        f.write(f"{r},{i},{s}\n")

    f.close()