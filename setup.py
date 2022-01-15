import numpy as np
import pandas as pd
from matplotlib import colors
import matplotlib.pyplot as plt

q = 4
Ts = np.linspace(0.8, 1.0, 21)
Ls = [6, 7, 8]
filename = 'bc-data-2fac-3.csv'

# np.random.seed(9)
# L = 4


def spin_map(array, label='randlabel'):
    L = np.shape(array)[0]
    # create discrete colormap
    cmap = colors.ListedColormap(['moccasin', 'orange', 'peru', 'saddlebrown'])
    bounds = [0, 1, 2, 3, 4]
    norm = colors.BoundaryNorm(bounds, cmap.N)

    fig, ax = plt.subplots()
    ax.imshow(array, cmap=cmap, norm=norm)

    # draw gridlines
    ax.grid(which='major', axis='both', linestyle='-', color='k', linewidth=2)
    ax.set_xticks(np.arange(0.5, L + .4, 1))
    ax.set_yticks(np.arange(0.5, L + .4, 1))

    plt.tick_params(bottom=False, left=False, labelbottom=False, labelleft=False)

    # plt.savefig('spin-map-'+label+'.png', bbox_inches='tight')
    plt.show()


def H(array):
    L = np.shape(array)[0]
    H = 0
    for i in range(L):
        for j in range(L):
            if array[i - 1, j] == array[i, j]:
                H -= 1
            if array[i, j - 1] == array[i, j]:
                H -= 1
    return H


def localH(array, i, j):
    L = np.shape(array)[0]
    # requires fewer calculations than H,
    # so it's useful for doing the spin flips more quickly
    H = 0
    # periodic boundary conditions
    # minx = i - 1 if i > 0 else L - 1
    # miny = j - 1 if i > 0 else L - 1
    maxx = i + 1 if i < L - 1 else 0
    maxy = j + 1 if j < L - 1 else 0
    for x in [i - 1, maxx]:
        if array[x, j] == array[i, j]:
            H -= 1
    for y in [j - 1, maxy]:
        if array[i, y] == array[i, j]:
            H -= 1
    return H


def m(array):
    """
    treating states like vectors
    """
    L = np.shape(array)[0]
    m = 0
    for qi in range(q):
        m += np.count_nonzero(array==qi) * np.exp(2j * np.pi * qi / q)
    return np.absolute(m) / L**2 / q


def m_altdef(array):
    L = np.shape(array)[0]
    mNp = max([np.count_nonzero(array==p) for p in range(q)])
    return (q * mNp / L**2 - 1) / (q - 1)


def Wolff(grid, T):
    """
    one iteration of the Wolff cluster flip
    """
    L = np.shape(grid)[0]
    sites_in_cluster = []
    i, j = np.random.choice(range(L), 2)
    sites_in_cluster.append([i, j])
    old_val = grid[i, j]

    p = 1 - np.exp(- 1 / T)
    for [i, j] in sites_in_cluster:
        # periodic boundary conditions
        minx = i - 1 if i > 0 else L - 1
        miny = j - 1 if j > 0 else L - 1
        maxx = i + 1 if i < L - 1 else 0
        maxy = j + 1 if j < L - 1 else 0
        # left and right
        for x in [minx, maxx]:
            # don't want to repeat with sites already in cluster
            if [x, j] not in sites_in_cluster and grid[i, j] == grid[x, j]:
                rand = np.random.rand()
                if rand < p:
                    sites_in_cluster.append([x, j])
        # up and down
        for y in [miny, maxy]:
            # don't want to repeat with sites already in cluster
            if [i, y] not in sites_in_cluster and grid[i, j] == grid[i, y]:
                rand = np.random.rand()
                if rand < p:
                    sites_in_cluster.append([i, y])

    # want to flip cluster to different value
    l = list(range(q))
    l.remove(old_val)
    new_val = np.random.choice(l)
    for [i, j] in sites_in_cluster:
        grid[i, j] = new_val.copy()

    # return np.shape(sites_in_cluster)[0]


def single_flip(grid, T):
    L = np.shape(grid)[0]
    flipsite = np.random.choice(range(L), 2)
    x, y = flipsite

    newgrid = grid.copy()

    old_val = grid[x, y].copy()
    l = list(range(q))
    l.remove(old_val)
    new_val = np.random.choice(l)
    newgrid[x, y] = new_val

    E = localH(grid, x, y)
    newE = localH(newgrid, x, y)
    p = min(1, np.exp((E - newE) / T))
    rand = np.random.rand()
    if rand < p:
        grid[x, y] = new_val.copy()


def binder_cumulant(Ls, Ts, m, filename): 
    """
    finding curves of binder cumulant for different Ls
    returns csv file
    """
    BC = np.zeros((np.size(Ls), np.size(Ts)))
    frac_count = 1 / 3
    # ntrials = 3000
    for Lspot in range(np.size(Ls)):
        L = Ls[Lspot]
        print(L)
        nflips = int(L**2 / 2) * 10**4
        for Tspot in range(np.size(Ts)):
            m4 = 0
            m2 = 0
            # for _ in range(ntrials):
            grid = np.random.choice(range(q), L**2).reshape((L, L))
            for i in range(nflips):
                Wolff(grid, Ts[Tspot])
                if i > nflips * frac_count:
                    m4 += m(grid)**4 / (nflips * (1 - frac_count))
                    m2 += m(grid)**2 / (nflips * (1 - frac_count))
            # if round(Ts[Tspot], 2) in [0.80, 0.90, 1.00, 1.10, 1.20, 1.30, 1.40]:
            #     spin_map(grid, str(L) + '-' + str(int(10 * Ts[Tspot])))
            BC[Lspot, Tspot] = 1 - m4 / m2**2 / 3

    a = pd.DataFrame(data=BC, columns=Ts, index=Ls)
    a.to_csv(filename)


def logistic(x, A, C, x0, k, N):
    """
    Assuming the BC curves follow logistic curves (therefore max is 2/3),
    two parameters are k and x0, data is x
    """
    return A + (C - A) / (1 + np.exp(k * (x - x0)))**(1 / N)


if __name__ == '__main__':
    binder_cumulant(Ls, Ts, m, filename)