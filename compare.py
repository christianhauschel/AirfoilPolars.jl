# %%

import airfoilprep as ap
import pandas as pd
from copy import copy
import proplot as pplt

df = pd.read_csv('/mnt/a/Code/10_drones/Propeller/data/polars/NACA4412/NACA4412_Re0.050_M0.00_N5.0.csv')
cd_max = 1.5

p = ap.Polar(50000, df.alpha.values, df.cl.values, df.cd.values, df.cm.values)

p_ext = p.extrapolate(cd_max)

fig, ax = pplt.subplots(figsize=(9,3), nrows=1, ncols=3, sharex=True, sharey=False)
ax[0].plot(p_ext.alpha, p_ext.cl, label='Extrapolated')
ax[0].plot(p.alpha, p.cl, label='Original')
ax[1].plot(p_ext.alpha, p_ext.cd, label='Extrapolated')
ax[1].plot(p.alpha, p.cd, label='Original')
ax[2].plot(p_ext.alpha, p_ext.cm, label='Extrapolated')
ax[2].plot(p.alpha, p.cm, label='Original')

# p_ext_3D = p_ext.correction3D(0.22, 0.53, 4.24)

p_ext_3D = p_ext.correction3D( 0.12037037037037036, 1.1626923076923077, 8.4823)

fig, ax = pplt.subplots(figsize=(9,3), nrows=1, ncols=3, sharex=True, sharey=False)
ax[0].plot(p_ext_3D.alpha, p_ext_3D.cl, label='3D Corrected')
ax[0].plot(p_ext.alpha, p_ext.cl, label='Extrapolated')
ax[1].plot(p_ext_3D.alpha, p_ext_3D.cd, label='3D Corrected')
ax[1].plot(p_ext.alpha, p_ext.cd, label='Extrapolated')
ax[2].plot(p_ext_3D.alpha, p_ext_3D.cm, label='3D Corrected')
ax[2].plot(p_ext.alpha, p_ext.cm, label='Extrapolated')


# %%