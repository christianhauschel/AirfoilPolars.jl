"""Convert QBlade / CCBlade polar ot Dust polar files"""

# %%

from pathlib import Path
import numpy as np
import proplot as pplt

# Parameters

name_airfoil = "0009"
dir_in = Path("data/polars/djimatrice300rtk/")
pattern_fname = f"{name_airfoil}_{name_airfoil}_*.plr"
dir_out = dir_in / "dust"

dummy_mach = np.linspace(0, 0.5, 2)

# %%

fnames = list(dir_in.glob(pattern_fname))
dir_out.mkdir(parents=True, exist_ok=True)

fname_out = dir_out / f"{name_airfoil}.dat"

def parse_file(filename):
    # get file extension from Path
    ext = filename.suffix[1:]

    if ext == "dat":
        filetype = "ccblade"
    elif ext == "plr":
        filetype = "qblade"
    else:
        raise ValueError("Unknown filetype")

    alpha = []
    cl = []
    cd = []
    cm = []
    info = ""
    Re = 1.0
    Mach = 0.0

    with open(filename, "r") as f:
        if filetype == "ccblade":
            info = next(f).strip()
            Re = float(next(f).strip())
            Mach = float(next(f).strip())
            raise NotImplementedError("CCBlade polar not working as it does not extrapolate cm!!!")
        elif filetype == "qblade":
            for _ in range(9):
                next(f)
            info = next(f).strip()
            for _ in range(3):
                next(f)
            Re = float(next(f).split()[1])
            for _ in range(3):
                next(f)

        for line in f:
            parts = line.split()
            alpha.append(float(parts[0]))
            cl.append(float(parts[1]))
            cd.append(float(parts[2]))
            cm.append(float(parts[3]))

    alpha = np.array(alpha)
    cl = np.array(cl)
    cd = np.array(cd)
    cm = np.array(cm)


    if max(abs(alpha)) > 4:  # Raw values of alpha are in [deg]
        alpha = [np.radians(a) for a in alpha]

    return info, Re, Mach, alpha, cl, cd, cm


# %%

from airfoil import sampling

n_alpha = 99
# alpha = np.linspace(-180, 180, n_alpha)
order = 10
alpha1 = 180 * (sampling("polynomial", -1, 1, n_alpha//2+1, m=np.pi/2, order=order) - 1)
alpha2 = -alpha1
alpha = np.concatenate((alpha1[:-1], np.flip(alpha2)))
alpha_rad = np.radians(alpha)

delta_alpha = np.diff(alpha)
print(delta_alpha)

n_files = len(fnames)
n_mach = len(dummy_mach)
# TODO
# nonlinear spacing
# spacing should be closer at -15 to 15 degrees range
cl = np.zeros((n_files, n_alpha, n_mach))
cd = np.zeros((n_files, n_alpha, n_mach))
cm = np.zeros((n_files, n_alpha, n_mach))
Re = np.zeros(n_files)
Mach = np.zeros(n_files)

for i, fname in enumerate(fnames):
    info, Re[i], Mach[i], _alpha, _cl, _cd, _cm = parse_file(fname)
    for j in range(n_mach):
        cl[i,:,j] = np.interp(alpha_rad, _alpha, _cl)
        cd[i,:,j] = np.interp(alpha_rad, _alpha, _cd)
        cm[i,:,j] = np.interp(alpha_rad, _alpha, _cm)

n_Re = len(Re)

# fig, ax = pplt.subplots(figsize=(7, 1))
# ax.plot(alpha, np.zeros_like(alpha), marker="o", ms=1, ls="none", label="cosine")
# ax.legend()


# fig, ax = pplt.subplots(figsize=(5, 4))
# ax.plot(_alpha, _cl)
# ax.plot(alpha_rad, cl[-1,:], ".-", lw=0.5)    


# %%

# write file
fname_out = dir_out / f"{name_airfoil}.c81"

with open(fname_out, "w") as f:
    f.write(f"{n_Re} 0 0\n")
    f.write(f"0 1\n")
    f.write(f"0.158 0.158\n")
    f.write(f"COMMENT#1\n")
    

    for j in range(n_Re):
        f.write(f"  {Re[j]:0.2f} 0.2\n")
        f.write(f"Vahana1        {n_mach:02d}{n_alpha:02d}{n_mach:02d}{n_alpha:02d}{n_mach:02d}{n_alpha:02d}\n")
        
        for ma in dummy_mach:
            f.write(f"\t{ma:6.3f}")
        f.write("\n")
        for i, a in enumerate(alpha):
            f.write(f"{a}")
            for k in range(n_mach):
                f.write(f"\t{cl[j,i,k]:6.3f}")
            f.write("\n")

        for ma in dummy_mach:
            f.write(f"\t{ma:6.3f}")
        f.write("\n")
        for i, a in enumerate(alpha):
            f.write(f"{a}")
            for k in range(n_mach):
                f.write(f"\t{cd[j,i,k]:6.3f}")
            f.write("\n")

        for ma in dummy_mach:
            f.write(f"\t{ma:6.3f}")
        f.write("\n")
        for i, a in enumerate(alpha):
            f.write(f"{a}")
            for k in range(n_mach):
                f.write(f"\t{cm[j,i,k]:6.3f}")
            f.write("\n")
  
        f.write("\n")
   


# %%