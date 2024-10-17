import os, itertools
import numpy as np
import mpmath as mpf
import tqdm

mpf.mp.dps = 256

fs = {"besselj": mpf.besselj, "bessely": mpf.bessely, "besseli": mpf.besseli, "besselk": mpf.besselk}
nus = np.arange(0, 64+1, 1) / 4
xs = [0, 1/8, 1/4, 1/2, 1, 2, 4, 8, 16, 32, 64, 128, 256, 
      -1/8, -1/4, -1/2, -1, -2, -4, -8, -16, -32, -64, -128, -256]

for nu_abs in tqdm.tqdm(nus, "Processing...", total=len(nus)):
    for fname in fs.keys():
        dirpath = f"../mpmath_expected/{fname}_plusi/"
        os.makedirs(dirpath, exist_ok=True)

        func = fs[fname]

        for nu in [nu_abs, -nu_abs] if nu_abs > 0 else [0]:
            if fname == "besselk" and nu < 0:
                continue

            nu_str = str(nu).replace(".", "p").replace("-", "m")

            f = open(f"{dirpath}/{fname}_nu{nu}.csv", "w")

            f.write("r,i,z\n")

            try:
                for r, i in itertools.product(xs, xs):
                    if r <= 0 and i == 0:
                        continue
                    if i < 0:
                        continue
                    
                    c = r + i * 1j

                    y = func(nu, c)

                    s = mpf.nstr(y, 40, strip_zeros=False, min_fixed=2, max_fixed=1)
                    s = s.replace(" ", "").replace("j", "i").replace("(", "").replace(")", "")

                    f.write(f"{r},{i},{s}\n")
            except:
                pass

            f.close()