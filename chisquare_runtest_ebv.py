# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
from tabulate import tabulate

import chisquare_v4 as cs

# Calculando os parâmetros
nomes = np.loadtxt('results/photometry-cat.txt', usecols=[0], dtype=str)
fields = np.loadtxt('results/photometry-cat.txt', usecols=[1], dtype=str)
literature_param = np.loadtxt('results/parameter-cat.txt', usecols=(5, 3, 7))
literature_name = np.loadtxt('results/parameter-cat.txt', usecols=(0), dtype=str)
print('Object field [Fe/H](Lit, Calc, -, +) log(Age)(Lit, Calc, -, +) E(B-V)(Calc, -, +)')
for i in range(0, len(nomes)):
    cluster_locator = 0  # locating the cluster in the literature table
    while (literature_name[cluster_locator] != nomes[i]):
        cluster_locator = cluster_locator + 1

    cluster = np.loadtxt('results/photometry-cat.txt',
                         usecols=[2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,
                                  25, 26, 27, 28])[i]
    # print(cluster)
    Nsint = 2
    sint_cluster = np.zeros((Nsint, 27))
    sint_cluster[0, :] = cluster
    sint_cluster_param = np.zeros((Nsint, 3))
    for k in range(1, Nsint):
        sint_cluster[k, 0:3] = cluster[0:3]
        for m in range(0, 24):
            if m < 12:
                sint_cluster[k, 3 + m] = np.random.normal(cluster[3 + m], 3 * cluster[15 + m], 1)
                # print(cluster[15+m])
            else:
                sint_cluster[k, 3 + m] = cluster[3 + m]
            if sint_cluster[k, 3 + m] < 0:
                sint_cluster[k, 3 + m] = -sint_cluster[k, 3 + m]
    # print(sint_cluster[:,3:15])
    # break
    for j in range(0, len(sint_cluster)):
        chi_ebv, chi_Z, chi_age = cs.full_fit(cs.load_table('tables/mag_new.dat'),
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
        # print (nomes[i], chi_age[0:3],chi_Z[0:3],chi_ebv[0:3])

        # Ajuste com aglomerados desavermelhados, usar tabela mag1.dat (mudar no cabeçalho da função)
        # chi_Z=np.log10((chi_Z[0]/(1-chi_Z[0]-(0.2485+(1.78*chi_Z[0]))))/(0.0207)),np.log10(chi_Z[1]),chi_Z[2]
        # sint_cluster_param[j,:]=chi_Z
        # print(nomes[i],chi_Z[0:3],np.abs(chi_Z[0]-literature[cluster_locator,0]),
        # np.abs(chi_Z[1]-literature[cluster_locator,1]))

        # Ajuste com aglomerados avermelhados, usar tabela mag.dat (mudar no cabeçalho da função)
        chi_ebv = np.log10((chi_ebv[0] / (1 - chi_ebv[0] - (0.2485 + (1.78 * chi_ebv[0])))) / (0.0207)), np.log10(
            chi_ebv[1]), chi_ebv[2]
        sint_cluster_param[j, :] = chi_ebv
        # print(nomes[i], chi_ebv[0:3], np.abs(chi_ebv[0] - literature[cluster_locator, 0]),
        # np.abs(chi_ebv[1] - literature[cluster_locator, 1]))

    np.set_printoptions(precision=2)

    Z_percentiles = [np.percentile(sint_cluster_param[:, 0], 50, interpolation='nearest'),
                     np.percentile(sint_cluster_param[:, 0], 16, interpolation='nearest'),
                     np.percentile(sint_cluster_param[:, 0], 84, interpolation='nearest')]
    age_percentiles = [np.percentile(sint_cluster_param[:, 1], 50, interpolation='nearest'),
                       np.percentile(sint_cluster_param[:, 1], 16, interpolation='nearest'),
                       np.percentile(sint_cluster_param[:, 1], 84, interpolation='nearest')]
    ebv_percentiles = [np.percentile(sint_cluster_param[:, 2], 50, interpolation='nearest'),
                       np.percentile(sint_cluster_param[:, 2], 16, interpolation='nearest'),
                       np.percentile(sint_cluster_param[:, 2], 84, interpolation='nearest')]
    # print(Z_percentiles, age_percentiles)
    print("{} {} ({:.2f} {:.2f} {:.2f} {:.2f}) ({:.2f} {:.2f} {:.2f} {:.2f}) ({:.2f} {:.2f} {:.2f} {:.2f})".format(
        nomes[i], fields[i],
        literature_param[cluster_locator, 0],
        Z_percentiles[0], Z_percentiles[1], Z_percentiles[2],
        literature_param[cluster_locator, 1],
        age_percentiles[0], age_percentiles[1], age_percentiles[2],
        literature_param[cluster_locator, 2],
        ebv_percentiles[0], ebv_percentiles[1], ebv_percentiles[2]))

    if i == 0:
        with open('results/Z-age-ebv/clusters_output.txt', 'w') as f:
            f.write(
                '#Object field Lit_[Fe/H] Calc_[Fe/H] -_[Fe/H] +_[Fe/H] Lit_log(Age) Calc_log(Age) -_log(Age) +_log(Age) Lit_E(B-V) Calc_E(B-V) -_E(B-V) +_E(B-V)\n')
            f.write("{} {} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f}\n"
                    .format(nomes[i], fields[i],
                            literature_param[cluster_locator, 0],
                            Z_percentiles[0], Z_percentiles[1], Z_percentiles[2],
                            literature_param[cluster_locator, 1],
                            age_percentiles[0], age_percentiles[1], age_percentiles[2],
                            literature_param[cluster_locator, 2],
                            ebv_percentiles[0], ebv_percentiles[1], ebv_percentiles[2]))
            # f.write('{} {:.2f} {:.2f} {:.2f} {:.4f} {:.2f}\n'.format(nomes[i], chi_Z[0], chi_Z[1], chi_Z[2],
            # np.abs(chi_Z[0] - literature[i, 0]),
            # np.abs(chi_Z[1] - literature[i, 0])))
    else:
        with open('results/Z-age-ebv/clusters_output.txt', 'a') as f:
            f.write("{} {} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f}\n"
                    .format(nomes[i], fields[i],
                            literature_param[cluster_locator, 0],
                            Z_percentiles[0], Z_percentiles[1], Z_percentiles[2],
                            literature_param[cluster_locator, 1],
                            age_percentiles[0], age_percentiles[1], age_percentiles[2],
                            literature_param[cluster_locator, 2],
                            ebv_percentiles[0], ebv_percentiles[1], ebv_percentiles[2]))
            # f.write('{} {:.2f} {:.2f} {:.2f} {:.4f} {:.2f}\n'.format(nomes[i], chi_Z[0], chi_Z[1], chi_Z[2],
            #                                                        np.abs(chi_Z[0] - literature[i, 0]),
            #                                                        np.abs(chi_Z[1] - literature[i, 1])))

    #######################################################################################################################
    #######################################################################################################################
    #######################################################################################################################
    #######################################################################################################################

    plt.subplots()  # (figsize=(12, 12))
    # plt.suptitle(fields[i]+'_'+nomes[i])
    plt.subplot(3, 1, 1)
    plt.hist(sint_cluster_param[:, 0], bins=np.arange((Z_percentiles[1] - 1), (Z_percentiles[2] + 1), 0.2))
    plt.vlines(literature_param[cluster_locator, 0], 0, Nsint / 2, label='literature', color='black')
    plt.vlines(Z_percentiles[0], 0, Nsint / 3, label='median', color='orange')
    plt.vlines(Z_percentiles[1], 0, Nsint / 4, label='16%', color='red')
    plt.vlines(Z_percentiles[2], 0, Nsint / 4, label='84%', color='green')
    # plt.legend(bbox_to_anchor=(1.0,1.05), loc="upper left", frameon=False)
    plt.xlabel('[Fe/H]', fontsize=15)
    plt.ylabel('N', fontsize=15)
    plt.tick_params(labelsize=15)

    plt.subplot(3, 1, 2)
    plt.hist(sint_cluster_param[:, 1], bins=np.arange((age_percentiles[1] - 0.15), (age_percentiles[2] + 0.15), 0.1))
    plt.vlines(literature_param[cluster_locator, 1], 0, Nsint / 2, label='literature', color='black')
    plt.vlines(age_percentiles[0], 0, Nsint / 3, label='median', color='orange')
    plt.vlines(age_percentiles[1], 0, Nsint / 4, label='16%', color='red')
    plt.vlines(age_percentiles[2], 0, Nsint / 4, label='84%', color='green')
    # plt.legend(bbox_to_anchor=(1.0,1.05), loc="upper left", frameon=False)
    plt.xlabel('log(Age)', fontsize=15)
    plt.ylabel('N', fontsize=15)
    plt.tick_params(labelsize=15)
    plt.tight_layout()

    plt.subplot(3, 1, 3)
    plt.hist(sint_cluster_param[:, 2], bins=np.arange((ebv_percentiles[1] - 0.15), (ebv_percentiles[2] + 0.15), 0.1))
    plt.vlines(literature_param[cluster_locator, 1], 0, Nsint / 2, label='literature', color='black')
    plt.vlines(ebv_percentiles[0], 0, Nsint / 3, label='median', color='orange')
    plt.vlines(ebv_percentiles[1], 0, Nsint / 4, label='16%', color='red')
    plt.vlines(ebv_percentiles[2], 0, Nsint / 4, label='84%', color='green')
    plt.legend(bbox_to_anchor=(1.0, 1.05), loc="upper left", frameon=False)
    plt.xlabel('E(B-V)', fontsize=15)
    plt.ylabel('N', fontsize=15)
    plt.tick_params(labelsize=15)

    plt.tight_layout()
    plt.savefig('results/Z-age-ebv/mcmc_hists/' + fields[i] + '_' + nomes[i] + '.png', dpi=300, bbox_inches='tight')
    plt.close()
    # plt.show()

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################


results_strings = np.loadtxt('results/Z-age-ebv/clusters_output.txt', usecols=[0, 1], dtype=str)
results_floats = np.loadtxt('results/Z-age-ebv/clusters_output.txt', usecols=[2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13],
                            dtype=float)
results = np.hstack((results_strings, results_floats))

results_latex = tabulate(results, tablefmt="latex",
                         headers=["Object", "field", "Lit_[Fe/H]", "Calc_[Fe/H]", "-_[Fe/H]",
                                  "+_[Fe/H]", "Lit_log(Age)", "Calc_log(Age)", "-_log(Age)", "+_log(Age)",
                                  "Lit_E(B-V)", "Calc_E(B-V)", "-_E(B-V)", "+_E(B-V)"])

with open('results/Z-age-ebv/clusters_output_latex.txt', 'w') as f:
    f.write("{}".format(results_latex))

fig0, ax0 = plt.subplots(ncols=1)
x = np.arange(0, len(results_strings))
ax0.set_xlabel('Clusters')
ax0.set_ylabel(r'|$\Delta$[Fe/H]|', fontsize=15)
ax0.set_xticks(x)
ax0.set_xticklabels((results_strings[:, 0]), rotation='vertical', fontsize=8)
# plt.tick_params(labelsize=15)
for i in range(len(results_strings)):
    # print(i,np.abs(results_floats[i, 0] - results_floats[i, 1]))

    ax0.scatter(results_strings[i, 0], (results_floats[i, 0] - results_floats[i, 1]), marker='o',
                label=results_strings[i, 1] + '_' + results_strings[i, 0])
    # if(results_floats[i,2]!=results_floats[i,3]):
    #    ax0.vlines(results_strings[i,0],results_floats[i,2],results_floats[i,3])
plt.tight_layout()
plt.savefig('results/Z-age-ebv/clusters_Z.png', dpi=300, bbox_inches='tight')

fig1, ax1 = plt.subplots(ncols=1)
# plt.subplots(121)
x = np.arange(0, len(results_strings))
ax1.set_xlabel('Clusters')
ax1.set_ylabel(r'|$\Delta$log(age)|', fontsize=15)
ax1.set_xticks(x)
ax1.set_xticklabels(results_strings[:, 0], rotation='vertical', fontsize=8)  # adicionar campo no tick também?
for i in range(len(results_strings)):
    ax1.scatter(i, (results_floats[i, 4] - results_floats[i, 5]), marker='o',
                label=results_strings[i, 1] + '_' + results_strings[i, 0])
    # plt.xticks(results_strings[i,1],rotation=55)
    # if(results_floats[i,2]!=results_floats[i,3]):
    #    ax1.vlines(results_strings[i,0],results_floats[i,2],results_floats[i,3])
plt.tight_layout()
plt.savefig('results/Z-age-ebv/clusters_age.png', dpi=300, bbox_inches='tight')

fig2, ax2 = plt.subplots(ncols=1)
# plt.subplots(121)
x = np.arange(0, len(results_strings))
ax2.set_xlabel('Clusters')
ax2.set_ylabel(r'|$\Delta$E(B-V)|', fontsize=15)
ax2.set_xticks(x)
ax2.set_xticklabels(results_strings[:, 0], rotation='vertical', fontsize=8)  # adicionar campo no tick também?
for i in range(len(results_strings)):
    ax2.scatter(i, (results_floats[i, 8] - results_floats[i, 9]), marker='o',
                label=results_strings[i, 1] + '_' + results_strings[i, 0])
    # plt.xticks(results_strings[i,1],rotation=55)
    # if(results_floats[i,2]!=results_floats[i,3]):
    #    ax1.vlines(results_strings[i,0],results_floats[i,2],results_floats[i,3])
plt.tight_layout()
plt.savefig('results/Z-age-ebv/clusters_ebv.png', dpi=300, bbox_inches='tight')
