import numpy as np

coord = np.loadtxt("tables/cluster-bica.txt", usecols=[1, 2], dtype=float)
param = np.loadtxt("tables/cluster-bica.txt", usecols=[3, 4, 5, 6], dtype=float)
cluster = np.loadtxt("tables/cluster-bica.txt", usecols=[0], dtype=str)
field_name = np.loadtxt("tables/cluster-bica.txt", usecols=[7], dtype=str)
field_param = np.loadtxt("tables/cluster-bica.txt", usecols=[8, 9, 10, 11])
ebv_table = np.loadtxt("tables/smc_res5_radius7.txt", usecols=[0, 1, 2])

k = 0
for i in range(len(cluster)):
    d = np.zeros(len(ebv_table))
    for j in range(len(ebv_table)):
        d[j] = np.sqrt((coord[i, 0] - ebv_table[j, 0]) ** 2 + (coord[i, 1] - ebv_table[j, 1]) ** 2)
    d_min = np.min(d)
    print(d_min)
    if (d_min < 0.11):  # equivalente a um raio de 7 arcmin
        if k == 0:
            with open('results/parameter-cat.txt', 'w') as f:
                f.write('{} {} {} {} {} {} {} {} {} {} {} {} {}'.format(cluster[i], coord[i, 0], coord[i, 1],
                                                                        param[i, 0], param[i, 1], param[i, 2],
                                                                        param[i, 3],
                                                                        ebv_table[np.argmin(d), 2],
                                                                        field_name[i], field_param[i, 0],
                                                                        field_param[i, 1],
                                                                        field_param[i, 2], field_param[i, 3]))
                k = 1
        else:
            with open('results/parameter-cat.txt', 'a') as f:
                f.write('\n{} {} {} {} {} {} {} {} {} {} {} {} {}'.format(cluster[i], coord[i, 0], coord[i, 1],
                                                                          param[i, 0], param[i, 1], param[i, 2],
                                                                          param[i, 3],
                                                                          ebv_table[np.argmin(d), 2],
                                                                          field_name[i], field_param[i, 0],
                                                                          field_param[i, 1],
                                                                          field_param[i, 2], field_param[i, 3]))
