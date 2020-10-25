# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np

import chisquare_v4 as cs

# Calculando os parâmetros
literature_param = np.loadtxt('tables/mag2.dat', usecols=(0, 1, 2))
literature_param[:, 0] = np.log10(
    (literature_param[:, 0] / (1 - literature_param[:, 0] - (0.2485 + (1.78 * literature_param[:, 0])))) / (0.0207))
literature_param[:, 1] = np.log10(literature_param[:, 1])

# literature_name=np.loadtxt('results/parameter-cat.txt',usecols=(0),dtype=str)
print('Object [Fe/H](Lit, Calc, -, +) log(Age)(Lit, Calc, -, +)')

for i in range(0, len(literature_param)):

    cluster_phot = np.loadtxt('tables/mag2.dat',
                              usecols=[3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14])[i]
    cluster_err = 0.02 * np.ones((12))
    cluster_zeros = np.zeros((3))
    cluster = np.hstack((cluster_zeros, cluster_phot, cluster_err))
    # print(cluster)
    Nsint = 1000
    sint_cluster = np.zeros((Nsint, 27))
    sint_cluster[0, :] = cluster
    sint_cluster_param = np.zeros((Nsint, 3))
    for k in range(1, Nsint):
        sint_cluster[k, 0:3] = cluster[0:3]
        for m in range(0, 24):
            if m < 12:
                sint_cluster[k, 3 + m] = np.random.normal(cluster[3 + m], 1 * cluster[15 + m], 1)
            else:
                sint_cluster[k, 3 + m] = cluster[3 + m]
            if sint_cluster[k, 3 + m] < 0:
                sint_cluster[k, 3 + m] = -sint_cluster[k, 3 + m]
    for j in range(0, len(sint_cluster)):
        chi_ebv, chi_Z, chi_age = cs.full_fit(cs.load_table('tables/mag1.dat'),
                                              sint_cluster, j,
                                              (['uJava', 'uJava', 'uJava', 'uJava', 'F378', 'F378', 'F430', 'gSDSS',
                                                'F660', 'F861', 'uJava']),
                                              (['F660', 'rSDSS', 'iSDSS', 'zSDSS', 'F515', 'gSDSS', 'gSDSS', 'zSDSS',
                                                'gSDSS', 'uJava', 'F515']),
                                              # cs.load_table('best5_100age.dat')[-5,0:11],
                                              np.ones(11),
                                              (['iSDSS', 'rSDSS', 'uJava', 'F378', 'F410', 'F430', 'F515', 'F660',
                                                'F660', 'gSDSS', 'iSDSS']),
                                              (['zSDSS', 'zSDSS', 'zSDSS', 'zSDSS', 'zSDSS', 'iSDSS', 'zSDSS', 'F861',
                                                'iSDSS', 'zSDSS', 'F861']),
                                              np.ones(11),
                                              # cs.load_table('best5_100Z_test6.dat')[-5,0:11])1:29
                                              (
                                                  ['F660', 'gSDSS', 'rSDSS', 'F378', 'F378', 'F410', 'F410', 'F660',
                                                   'F410',
                                                   'F378', 'iSDSS']),
                                              (['iSDSS', 'F861', 'F861', 'F395', 'gSDSS', 'F430', 'F515', 'F861',
                                                'gSDSS', 'F660', 'F861']),
                                              np.ones(11))
        # cs.load_table('best5_100Z_test6.dat')[-5,0:11])tables

        np.set_printoptions(precision=2)

        # Ajuste com aglomerados desavermelhados, usar tabela mag1.dat (mudar no cabeçalho da função)
        chi_Z = np.log10((chi_Z[0] / (1 - chi_Z[0] - (0.2485 + (1.78 * chi_Z[0])))) / (0.0207)), np.log10(chi_Z[1]), \
                chi_Z[2]
        sint_cluster_param[j, :] = chi_Z

        # Ajuste com aglomerados avermelhados, usar tabela mag.dat (mudar no cabeçalho da função)
        # chi_ebv=np.log10((chi_ebv[0]/(1-chi_ebv[0]-(0.2485+(1.78*chi_ebv[0]))))/(0.0207)),np.log10(chi_ebv[1]),chi_ebv[2]
        # sint_cluster_param[j,:]=chi_ebv

    np.set_printoptions(precision=2)
    # percentils de idade e metalicidade
    Z_percentiles = [np.percentile(sint_cluster_param[:, 0], 50, interpolation='nearest'),
                     np.percentile(sint_cluster_param[:, 0], 16, interpolation='nearest'),
                     np.percentile(sint_cluster_param[:, 0], 84, interpolation='nearest')]
    age_percentiles = [np.percentile(sint_cluster_param[:, 1], 50, interpolation='nearest'),
                       np.percentile(sint_cluster_param[:, 1], 16, interpolation='nearest'),
                       np.percentile(sint_cluster_param[:, 1], 84, interpolation='nearest')]

    print("{} ({:.2f} {:.2f} {:.2f} {:.2f}) ({:.2f} {:.2f} {:.2f} {:.2f})".format(i,
                                                                                  literature_param[i, 0],
                                                                                  Z_percentiles[0], Z_percentiles[1],
                                                                                  Z_percentiles[2],
                                                                                  literature_param[i, 1],
                                                                                  age_percentiles[0],
                                                                                  age_percentiles[1],
                                                                                  age_percentiles[2]))

    if i == 0:
        with open('results/sints/clusters_output.txt', 'w') as f:
            f.write('#Lit_[Fe/H] Calc_[Fe/H] -_[Fe/H] +_[Fe/H] Lit_log(Age) Calc_log(Age) -_log(Age) +_log(Age) \n')
            f.write("{} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f}\n".format(i,
                                                                                          literature_param[i, 0],
                                                                                          Z_percentiles[0],
                                                                                          Z_percentiles[1],
                                                                                          Z_percentiles[2],
                                                                                          literature_param[i, 1],
                                                                                          age_percentiles[0],
                                                                                          age_percentiles[1],
                                                                                          age_percentiles[2]))

    else:
        with open('results/sints/clusters_output.txt', 'a') as f:
            f.write("{} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f}\n".format(i,
                                                                                          literature_param[i, 0],
                                                                                          Z_percentiles[0],
                                                                                          Z_percentiles[1],
                                                                                          Z_percentiles[2],
                                                                                          literature_param[i, 1],
                                                                                          age_percentiles[0],
                                                                                          age_percentiles[1],
                                                                                          age_percentiles[2]))

    plt.subplots()

    plt.subplot(1, 2, 1)
    plt.hist(sint_cluster_param[:, 0], bins=np.arange((Z_percentiles[1] - 1), (Z_percentiles[2] + 1), 0.2))
    plt.vlines(literature_param[i, 0], 0, Nsint / 2, label='literature', color='black')
    plt.vlines(Z_percentiles[0], 0, Nsint / 3, label='median', color='orange')
    plt.vlines(Z_percentiles[1], 0, Nsint / 4, label='16%', color='red')
    plt.vlines(Z_percentiles[2], 0, Nsint / 4, label='84%', color='green')
    plt.legend()
    plt.xlabel('[Fe/H]', fontsize=15)
    plt.ylabel('N', fontsize=15)
    plt.tick_params(labelsize=15)

    plt.subplot(1, 2, 2)
    plt.hist(sint_cluster_param[:, 1], bins=np.arange((age_percentiles[1] - 0.15), (age_percentiles[2] + 0.15), 0.1))
    plt.vlines(literature_param[i, 1], 0, Nsint / 2, label='literature', color='black')
    plt.vlines(age_percentiles[0], 0, Nsint / 3, label='median', color='orange')
    plt.vlines(age_percentiles[1], 0, Nsint / 4, label='16%', color='red')
    plt.vlines(age_percentiles[2], 0, Nsint / 4, label='84%', color='green')
    plt.legend()
    plt.xlabel('log(Age)', fontsize=15)
    plt.ylabel('N', fontsize=15)
    plt.tick_params(labelsize=15)
    plt.tight_layout()
    plt.savefig('results/sints/mcmc_hists/cluster' + str(i) + '.png', dpi=300, bbox_inches='tight')
    plt.close()

results_floats = np.loadtxt('results/sints/clusters_output.txt', usecols=[0, 1, 2, 3, 4, 5, 6, 7, 8], dtype=float)

fig0, ax0 = plt.subplots(ncols=1)
# x=np.arange(0,len(results_strings))
ax0.set_xlabel('Clusters')
ax0.set_ylabel(r'|$\Delta$[Fe/H]|', fontsize=15)

for i in range(len(results_floats)):
    ax0.scatter(results_floats[i, 0], (results_floats[i, 1] - results_floats[i, 2]), marker='o')

plt.tight_layout()
plt.savefig('results/sints/clusters_Z.png', dpi=300, bbox_inches='tight')

fig1, ax1 = plt.subplots(ncols=1)
ax1.set_xlabel('Clusters')
ax1.set_ylabel(r'|$\Delta$log(age)|', fontsize=15)

for i in range(len(results_floats)):
    ax1.scatter(results_floats[i, 0], (results_floats[i, 5] - results_floats[i, 6]), marker='o')

plt.tight_layout()
plt.savefig('results/sints/clusters_age.png', dpi=300, bbox_inches='tight')
