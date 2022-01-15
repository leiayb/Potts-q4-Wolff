    L = np.shape(array)[0]
    sites_in_cluster = []
    i, j = np.random.choice(range(L), 2)
    sites_in_cluster.append([i, j])

    p = 1 - np.exp(-2 / T)
    for [i, j] in sites_in_cluster:
        # periodic boundary conditions
        minx = i - 1 if i > 0 else L - 1
        miny = i - 1 if i > 0 else L - 1
        maxx = i + 1 if i < L - 1 else 0
        maxy = j + 1 if j < L - 1 else 0
        # left and right
        for x in [minx, maxx]:
            # don't want to repeat with sites already in cluster
            if [x, j] not in sites_in_cluster and array[i, j] == array[x, j]:
                rand = np.random.rand()
                if rand < p:
                    sites_in_cluster.append([x, j])
        # up and down
        for y in [miny, maxy]:
            # don't want to repeat with sites already in cluster
            if [i, y] not in sites_in_cluster and array[i, j] == array[i, y]:
                rand = np.random.rand()
                if rand < p:
                    sites_in_cluster.append([i, y])


    # want to flip cluster to different value
    old_val = array[sites_in_cluster[0][0], sites_in_cluster[0][1]]
    l = list(range(q))
    l.remove(old_val)
    new_val = np.random.choice(l)
    for [i, j] in sites_in_cluster:
        array[i, j] = new_val.copy()