import os, itertools
import numpy as np
import mpmath as mpf
import tqdm

mpf.mp.dps = 512

fs = {"besselj": mpf.besselj, "bessely": mpf.bessely, "besseli": mpf.besseli, "besselk": mpf.besselk}
nus = [16.25, 16.5, 16.75, 17, 17.25, 17.5, 17.75, 18, 
       18.25, 18.5, 18.75, 19, 19.5, 20, 20.25, 20.5, 20.75, 
       21, 63.75, 64, 127.75, 128, 255.5, 255.75, 256,
       -16.25, -16.5, -16.75, -17, -17.25, -17.5, -17.75, -18, 
       -18.25, -18.5, -18.75, -19, -19.5, -20, -20.25, -20.5, -20.75, 
       -21, -63.75, -64, -127.75, -128, -255.5, -255.75, -256]
xs = [1/8, 1/4, 1/2, 1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 64, 128, 256]

for fname in fs.keys():
    dirpath = f"../mpmath_expected/{fname}_real_nuplus16/"
    os.makedirs(dirpath, exist_ok=True)

    func = fs[fname]

    f_com = open(f"{dirpath}/{fname}_real.txt", "w")

    for nu in tqdm.tqdm(nus, f"Processing... {fname}", total=len(nus)):

        if fname == "besselk" and nu < 0:
            continue

        nu_str = str(nu).replace(".", "p").replace("-", "m")

        f_com.write(f"ddouble[] nu{nu_str}_expecteds = {{\n")

        f = open(f"{dirpath}/{fname}_real_nuplus16_nu{nu}.csv", "w")

        try:
            for x in xs:
                y = func(nu, x)
                f.write(f"{mpf.nstr(y, 40, strip_zeros=False, min_fixed=2, max_fixed=1)}\n")
                f_com.write(f"    \"{mpf.nstr(y, 40, strip_zeros=False, min_fixed=2, max_fixed=1)}\",\n")
        except:
            pass

        f_com.write("};\n")

        f.close()

    f_com.close()