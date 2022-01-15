import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import setup as s
from scipy.optimize import curve_fit

# reading, plotting, and fitting data
BC = pd.read_csv('bc-data-2fac-2.csv', index_col=0, header=0)
BC.columns = pd.to_numeric(BC.columns)
# print(BC)

colors = ['olivedrab', 'darkorange', 'purple']
fitnames = [str(i)+'fit' for i in s.Ls]
allp = []

for n, i in enumerate(BC.index):
    params, _ = curve_fit(s.logistic, BC.columns, BC.loc[i])
    allp.append(params)
    BC.loc[fitnames[n]] = s.logistic(BC.columns, *params)
    plt.plot(BC.columns, s.logistic(BC.columns, *params), color=colors[n])
    plt.plot(BC.loc[i], label='L = '+str(i), alpha=0.5, color=colors[n])
plt.legend()
plt.xlabel(r'$k_B T/J$', fontsize=15)
plt.ylabel(r'Binder cumulant $= 1 - \langle m^4 \rangle / 3 \langle m^2 \rangle^2$',
           fontsize=15)
plt.savefig('BC-originl.pdf', bbox_inches='tight')
plt.savefig('BC-originl.png', bbox_inches='tight')
# plt.show()


# finding intersection point of fits to find Tc
BC.loc['12diff'] = np.abs(BC.loc[fitnames[0]] - BC.loc[fitnames[1]])
BC.loc['13diff'] = np.abs(BC.loc[fitnames[0]] - BC.loc[fitnames[2]])
BC.loc['23diff'] = np.abs(BC.loc[fitnames[1]] - BC.loc[fitnames[2]])

# make sure to exclude near-intersection at low end of axis
Tcs = BC.iloc[:, 5:].idxmin(axis='columns')[-3:]
print(Tcs)
# Tc = np.average(Tcs)
Tc = 0.905
print(Tc)
print((Tc - 1 / np.log(3)) * np.log(3))

# guessing different nu's to make the three curves line up
alist = [1.45, 1.5, 1.55]
for a in alist:
    plt.figure()
    for n, i in enumerate(BC.index[:3]):
        plt.plot(i**a * (BC.columns - Tc), s.logistic(BC.columns, *allp[n]), color=colors[n], label='L = '+str(i))
    plt.title(r'$1/\nu = {}$'.format(a), fontsize=15)
    plt.xlabel(r'$L^{1/\nu} ( T - T_c )$', fontsize=15)
    plt.ylabel('Binder cumulant fit', fontsize=15)
    plt.legend()
    plt.savefig('BC-scaling'+str(int(100*a))+'.pdf', bbox_inches='tight')
    plt.savefig('BC-scaling'+str(int(100*a))+'.png', bbox_inches='tight')
    # plt.show()

# manually put in denominator:
# nu = 1 / 0.67
nu = 1 / 1.45
print((nu - 2 / 3) * 1.5)

