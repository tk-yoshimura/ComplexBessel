import os, itertools
import numpy as np
import mpmath as mpf
import tqdm

mpf.mp.dps = 128

fs = {"besselj": mpf.besselj, "bessely": mpf.bessely, "besseli": mpf.besseli, "besselk": mpf.besselk}
nus = [16.25, 16.5, 16.75, 17, 17.25, 17.5, 17.75, 18, 
       18.25, 18.5, 18.75, 19, 19.5, 20, 20.25, 20.5, 20.75, 
       21, 63.75, 64, 127.75, 128, 255.5, 255.75, 256,
       -16.25, -16.5, -16.75, -17, -17.25, -17.5, -17.75, -18, 
       -18.25, -18.5, -18.75, -19, -19.5, -20, -20.25, -20.5, -20.75, 
       -21, -63.75, -64, -127.75, -128, -255.5, -255.75, -256]
xs = [0, 1/8, 1/4, 1/2, 1, 2, 4, 8, 16, 32, 64, 128, 256]

for fname in fs.keys():
    dirpath = f"../mpmath_expected/{fname}_complex_nuplus16/"
    os.makedirs(dirpath, exist_ok=True)

    func = fs[fname]

    for nu in tqdm.tqdm(nus, f"Processing... {fname}", total=len(nus)):

        if fname == "besselk" and nu < 0:
            continue

        nu_str = str(nu).replace(".", "p").replace("-", "m")

        f = open(f"{dirpath}/{fname}_complex_nuplus16_nu{nu}.csv", "w")

        f.write("r,i,z\n")

        try:
            for r, i in itertools.product(xs, xs):
                if r == 0 and i == 0:
                    continue

                c = r + i * 1j

                y = func(nu, c)

                s = mpf.nstr(y, 40, strip_zeros=False, min_fixed=2, max_fixed=1)
                s = s.replace(" ", "").replace("j", "i").replace("(", "").replace(")", "")

                f.write(f"{r},{i},{s}\n")
        except:
            pass

        f.close()