import os

import matplotlib.pyplot as plt
import numpy as np

file = 'hl-photometry'  # nome do arquivo de fotometria gerado
param_table = "results/parameter-cat.txt"
cluster_list = np.loadtxt(param_table, usecols=(0), dtype=str)
field_list = np.loadtxt(param_table, usecols=(8), dtype=str)


def split_file(field, clustername):
    filename = "fits-splus/" + field + '/' + clustername
    filters = ['u', 'f378', 'f395', 'f410', 'f430', 'g', 'f515', 'r', 'f660', 'i', 'f861', 'z']
    if not os.path.exists("cluster-phot/" + field):
        os.makedirs("cluster-phot/" + field)
    if not os.path.exists("cluster-phot/" + field + "/" + clustername.upper()):
        os.makedirs("cluster-phot/" + field + "/" + clustername.upper())

    with open(filename, 'r') as f:
        # Reads the file with each line as an element of a list
        data = f.readlines()

        # put everything together in a single string
        data = "".join(data)

    # Get the header
    header = data.split(filters[0] + ".fits")[0]
    data = filters[0] + ".fits" + data.split(filters[0] + ".fits")[1]

    # Split and save table for each filter, last filter is done latter
    for i in range(len(filters) - 1):
        # Create new file and save table data
        with open("cluster-phot/" + field + "/" + clustername.upper() + "/" + clustername.upper() + "-" + filters[i],
                  'w') as f:
            f.write(header)
            f.write(data.split(filters[i + 1] + ".fits")[0])

        # Remove already saved table from data
        data = filters[i + 1] + ".fits" + data.split(filters[i + 1] + ".fits")[1]

    # Save table for the last filter
    with open("cluster-phot/" + field + "/" + clustername.upper() + "/" + clustername.upper() + "-" + filters[-1],
              'w') as f:
        f.write(header)
        f.write(data)


for i in range(len(cluster_list)):
    path = "fits-splus/" + field_list[i] + "/" + cluster_list[i].lower()
    if (os.path.isfile(path) == True):
        split_file(field_list[i], cluster_list[i].lower())


# APÓS FAZER A FOTOMETRIA, COLOCAR OS OUTPUTS DO DAOPHOT NA PASTA cluster-phot e rodar este programa, em seguida
# rodar o test.py, que realiza o ajuste dos parâmetros
def load_photometrytables(field, clustername):
    cluster = []
    # cluster=np.loadtxt(field + '/' + clustername, usecols=[0, 1, 2, 3, 4, 5, 6], skiprows=79)

    for i in colors:
        cluster.append(np.loadtxt('cluster-phot/' + field + '/' + clustername + '/'
                                  + clustername + '-' + i, usecols=[0, 1, 2, 3, 4, 5, 6], skiprows=79))
    return cluster


def magaper_graph(field, cluster, clustername):
    cmap = plt.get_cmap('nipy_spectral')
    for i in range(len(colors)):
        plt.plot(cluster[i][:, 0], cluster[i][:, 4], marker='o', label=colors[i], color=cmap(25 * i))
    plt.legend(bbox_to_anchor=(1.0, 1.05), loc="upper left", frameon=False)
    plt.xlabel('radius(pixel)', fontsize=15)
    plt.ylabel('mag', fontsize=15)
    # plt.xlim(0,np.max(cluster[0][:,0])+40)
    plt.title('{}_{}'.format(field, clustername))  # field+'_'+clustername)
    plt.tick_params(labelsize=15)
    plt.savefig('results/photometry/' + field + "_" + clustername + '.png', dpi=300, bbox_inches='tight')
    plt.close()
    # plt.show()


def cluster_photometrictable(cluster):
    mag = np.zeros(len(colors))
    mag_err = np.zeros(len(colors))
    min_aper_index = np.zeros(len(colors))

    for i in range(len(colors)):
        delta_maga = cluster[i][0, 4] - cluster[i][1, 4]

        for j in range(2, len(cluster[1])):
            delta_magb = cluster[i][j - 1, 4] - cluster[i][j, 4]
            if (delta_magb < delta_maga):
                delta_maga = delta_magb
                min_aper_index[i] = j - 1

            else:
                min_aper_index[i] = j - 1

    for i in range(len(colors)):
        mag[i] = cluster[i][int(np.min(min_aper_index) / 2), 4]
        mag_err[i] = cluster[i][int(np.min(min_aper_index / 2)), 5]
        aper = cluster[i][int(np.min(min_aper_index) / 2), 0]  # /0.55
    mag = list(mag)
    mag_err[mag_err < 0.0001] = np.mean(mag_err[:])
    mag_err = list(mag_err)
    cluster_phot = mag + mag_err
    cluster_phot = np.array(cluster_phot)
    cluster_phot = (cluster_phot)
    return cluster_phot, aper


# u.fits,f378.fits,f395.fits,f410.fits,f430.fits,g.fits,f515.fits,r.fits,f660.fits,i.fits,f861.fits,z.fits
param_table = "results/parameter-cat.txt"  # cat.txt
colors = ['u', 'f378', 'f395', 'f410', 'f430', 'g', 'f515', 'r', 'f660', 'i', 'f861', 'z']
cluster_list = np.loadtxt(param_table, usecols=(0), dtype=str)
field_list = np.loadtxt(param_table, usecols=(8), dtype=str)
# colors=['u','f378','f395','f410','f430','g','f515','r','f660','i','f861','z']
# cluster_list=['pal3']
kx = [4.916, 4.637, 4.467, 4.289, 4.091, 3.629, 3.325, 2.527, 2.317, 1.825, 1.470, 1.363]
# kx=[0,0,0,0,0,0,0,0,0,0,0,0] #E(B-V) como output. Mudar tabela para mag.dat em chisquare_runtest.py
a = np.zeros((len(cluster_list), 3 + 2 * len(colors)))
param = np.loadtxt(param_table, usecols=(1, 2, 4, 6, 7))  # pego o e(b-v) para desavermelhar o aglomerado
ZP_field = np.loadtxt('tables/MC_ZPs.cat', usecols=0, dtype=str)
ZP_table = np.loadtxt('tables/MC_ZPs.cat', usecols=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], dtype=float)
for e in range(0, len(cluster_list)):
    a[e, 0:2] = param[e, 0:2]
    # a[e,2]=rh[e]

n_cluster = 0
for i in range(0, len(cluster_list)):
    cluster = cluster_list[i]
    field = field_list[i]
    field_locator = 0
    for x in range(len(ZP_field)):
        if (ZP_field[x] == field):
            field_locator = x
    cluster_path = "cluster-phot/" + field_list[i] + "/" + cluster_list[i].upper()
    # print(cluster_path)
    if (os.path.exists(cluster_path) == True):
        print(cluster_path)
        magaper_graph(field, load_photometrytables(field, cluster), cluster)
        a[i, 3:], a[i, 2] = cluster_photometrictable(load_photometrytables(field, cluster))
        for u in range(len(cluster_list)):
            for v in range(3, len(colors) + 3):
                a[u, v] = a[u, v] - kx[v - 3] * param[u, 4]  # desavermelhamento do aglomerado
                a[u, v] = a[u, v] - (25 - ZP_table[field_locator, v - 3])
        # print(ZP_table[field_locator,:])
        if n_cluster == 0:
            n_cluster += 1
            # escrevendo catálogo de fotometrias
            with open('results/' + file, 'w') as f:
                f.write(
                    '#object field RA DEC radius_aper(arcsecs) uJava F378 F395 F410 F430 gSDSS F515 rSDSS F660 iSDSS F861 zSDSS e_uJava e_F378 e_F395 e_F410 e_F430 e_gSDSS e_F515 e_rSDSS e_F660 e_iSDSS e_F861 e_zSDSS\n')
                f.write('{} '.format(cluster))
                f.write('{} '.format(field))
                for j in range(0, 2 * len(colors) + 3):
                    f.write('{:.3f} '.format(a[i, j]))

        else:
            with open('results/' + file, 'a') as f:
                f.write('\n{} '.format(cluster))
                f.write('{} '.format(field))
                for j in range(0, 2 * len(colors) + 3):
                    f.write('{:.3f} '.format(a[i, j]))
