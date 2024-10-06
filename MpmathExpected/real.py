import os, itertools
import numpy as np
import mpmath as mpf
import tqdm

mpf.mp.dps = 128

fs = {"besselj": mpf.besselj, "bessely": mpf.bessely, "besseli": mpf.besseli, "besselk": mpf.besselk}
nus = np.arange(-64, 64+1, 1) / 4
xs = [1/256, 1/128, 1/64, 1/32, 1/16, 1/8, 1/4, 1/2, 1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 38.625, 38.75, 38.875, 64]

for fname in fs.keys():
    dirpath = f"../mpmath_expected/{fname}_real/"
    os.makedirs(dirpath, exist_ok=True)

    func = fs[fname]

    for x in tqdm.tqdm(xs, f"Processing... {fname}", total=len(xs)):

        f = open(f"{dirpath}/{fname}_real_x{x}.csv", "w")

        for nu in nus:
            y = func(nu, x)
            f.write(f"{mpf.nstr(y, 40, strip_zeros=False, min_fixed=2, max_fixed=1)}\n")

        f.close()