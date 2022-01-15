import numpy as np
import pandas as pd
import nuTc
import setup as s

q = 4
Ls = np.logspace(0.8, 1.25, 6)
# Ls = [6, 7, 8, 9, 10]
Ls = [int(L) for L in Ls]
# Tc = nuTc.Tc
Tc = 0.905
print('Tc =', Tc)
filename = 'main-meas-2fac-3.csv'

meas = np.zeros((3, np.size(Ls)))

def mainrun():
    for Lspot in range(np.size(Ls)):
        L = Ls[Lspot]
        print(L)
        m = 0
        m2 = 0
        E = 0
        E2 = 0
        grid = np.random.choice(range(q), L**2).reshape((L, L))
        nflips = int(L**2) * 10**3
        frac_count = 1 / 20
        for i in range(nflips):
            # if (L==7 or L==8) and (i==100 or i==101 or i==102 or i==103 or i==104 or i==105):
            #     s.spin_map(grid)
            s.Wolff(grid, Tc)
            # consider system flips to be thermalized after certain number of flips
            if i > nflips * frac_count:
                m += s.m(grid) / (nflips * (1 - frac_count))
                m2 += s.m(grid)**2 / (nflips * (1 - frac_count))
                E += s.H(grid) / (nflips * (1 - frac_count))
                E2 += s.H(grid)**2 / (nflips * (1 - frac_count))
        meas[0, Lspot] = m
        meas[1, Lspot] = (E2 - E**2) / (L**2 * Tc**2)
        meas[2, Lspot] = (m2 - m**2) * L**2 / Tc

    a = pd.DataFrame(data=meas, columns=Ls, index=['magnetization', 'heat capacity', 'susceptibility'])
    a.to_csv(filename)

if __name__ == '__main__':
    mainrun()