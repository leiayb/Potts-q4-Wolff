import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import setup as s

m_wolff = []
m_singl = []

L = 5
# Ts = np.linspace(0.8, 1.4, 61)
Ts = [0.8, 1.1, 1.4]
q = 4

frac_count = 1 / 3
nflips = 10**4
for T in Ts:
    T *= 2 
    mag = 0
    grid = np.random.choice(range(q), L**2).reshape((L, L))
    for i in range(nflips):
        s.Wolff(grid, T)
        if i > nflips * frac_count:
            mag += s.m(grid) / (nflips * (1 - frac_count))
        # if T == 1.4:
        #     s.spin_map(grid, 'rand')
    m_wolff.append(mag)

print('w done')
nflips = 10**6

for T in Ts:
    print(T)
    mag = 0
    grid = np.random.choice(range(q), L**2).reshape((L, L))
    for i in range(nflips):
        s.single_flip(grid, T)
        if i > nflips * frac_count:
            mag += s.m(grid) / (nflips * (1 - frac_count))
    m_singl.append(mag)
    print(mag)
    # s.spin_map(grid, str(int(10*T)))

# a = pd.DataFrame(data=data, columns=Ts, index=['wolff-singl'])
# a.to_csv('wolff-vs-single.csv')

plt.plot(Ts, m_wolff, label='Wolff')
plt.plot(Ts, m_singl, label='single')
plt.ylabel(r'order parameter $m$')
plt.xlabel(r'temperature $T$')
plt.legend()
plt.savefig('single-vs-wolff.png')