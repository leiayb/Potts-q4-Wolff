import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import nuTc
import mainrun as s

Ls = s.Ls
nu = nuTc.nu
# nu = 2 / 3

# read measurements from mainrun.py (m, C, chi VS L)
meas = pd.read_csv('main-meas-2fac-3.csv', index_col=0, header=0)
meas.columns = [int(col) for col in pd.to_numeric(meas.columns)]
# print(meas)

colors = ['olivedrab', 'darkorange', 'purple']
expnames = ['beta', 'alpha', 'gamma']
anlytc = [-1/12, 2/3, 7/6]
fitnames = [str(i)+'fit' for i in s.Ls]

# linear fit log-log plot and find slope, which is exponent/nu
for n, i in enumerate(meas.index):
    plt.figure()
    p = np.polyfit(np.log(meas.columns), np.log(meas.loc[i]), 1)  # linear fit
    print(expnames[n],'=', p[0] * nu, (anlytc[n] - p[0] * nu) / anlytc[n])
    plt.plot(np.log(meas.columns), p[0] * np.log(meas.columns) + p[1], color=colors[n])
    plt.scatter(np.log(meas.columns), np.log(meas.loc[i]), label=i, color=colors[n])
    plt.xlabel(r'$\ln\ L$', fontsize=15)
    plt.ylabel(r'$\ln$ '+i, fontsize=15)
    plt.savefig('exp-'+i+'.png', bbox_inches='tight')
    plt.savefig('exp-'+i+'.pdf', bbox_inches='tight')
    plt.show()

