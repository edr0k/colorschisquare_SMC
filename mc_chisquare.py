# -*- coding: utf-8 -*-
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn import linear_model
from tabulate import tabulate

import chisquare_v4 as cs


def mc_chisquare(pasta, age_delta, Nsint, Ns, phot='photometry_cat.txt', models='mag.dat'):
    if not (os.path.exists('results/' + pasta)):
        os.makedirs('results/' + pasta)
        os.makedirs('results/' + pasta + '/mcmc_hists')

    nomes = np.loadtxt('results/' + phot, usecols=[0], dtype=str)
    fields = np.loadtxt('results/' + phot, usecols=[1], dtype=str)
    literature_param = np.loadtxt('results/parameter-cat.txt', usecols=(5, 3, 7))
    literature_name = np.loadtxt('results/parameter-cat.txt', usecols=0, dtype=str)
    print(pasta)
    print('Object field [Fe/H](Lit,Mode,Median,-,+) log(Age)(Lit,Mode,Median,-,+)')
    for i in range(0, len(nomes)):
        cluster_locator = 0  # locating the cluster in the literature table
        while (literature_name[cluster_locator] != nomes[i]):
            cluster_locator = cluster_locator + 1

        cluster = np.loadtxt('results/' + phot,
                             usecols=[2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
                                      24,
                                      25, 26, 27, 28])[i]

        sint_cluster = np.zeros((Nsint, 27))
        sint_cluster[0, :] = cluster
        sint_cluster_param = np.zeros((Nsint, 3))
        for k in range(1, Nsint):
            sint_cluster[k, 0:3] = cluster[0:3]
            for m in range(0, 24):
                if m < 12:
                    sint_cluster[k, 3 + m] = np.random.normal(cluster[3 + m], Ns * cluster[15 + m], 1)
                    # print(cluster[15+m])
                else:
                    sint_cluster[k, 3 + m] = cluster[3 + m]
                if sint_cluster[k, 3 + m] < 0:
                    sint_cluster[k, 3 + m] = -sint_cluster[k, 3 + m]
        # print(sint_cluster[:,3:15])
        # break
        for j in range(0, len(sint_cluster)):
            chi_ebv, chi_Z, chi_age = cs.full_fit(cs.load_table('tables/' + models),
                                                  sint_cluster, j,
                                                  (['uJava', 'uJava', 'uJava', 'uJava', 'F378', 'F378', 'F430', 'gSDSS',
                                                    'F660', 'F861', 'uJava']),
                                                  (
                                                  ['F660', 'rSDSS', 'iSDSS', 'zSDSS', 'F515', 'gSDSS', 'gSDSS', 'zSDSS',
                                                   'gSDSS', 'uJava', 'F515']),
                                                  # cs.load_table('best5_100age.dat')[-5,0:11],
                                                  np.ones(11),
                                                  (['iSDSS', 'rSDSS', 'uJava', 'F378', 'F410', 'F430', 'F515', 'F660',
                                                    'F660', 'gSDSS', 'iSDSS']),
                                                  (['zSDSS', 'zSDSS', 'zSDSS', 'zSDSS', 'zSDSS', 'iSDSS', 'zSDSS',
                                                    'F861',
                                                    'iSDSS', 'zSDSS', 'F861']),
                                                  np.ones(11),
                                                  # cs.load_table('best5_100Z_test6.dat')[-5,0:11])1:29
                                                  (['F660', 'gSDSS', 'rSDSS', 'F378', 'F378', 'F410', 'F410', 'F660',
                                                    'F410', 'F378', 'iSDSS']),
                                                  (['iSDSS', 'F861', 'F861', 'F395', 'gSDSS', 'F430', 'F515', 'F861',
                                                    'gSDSS', 'F660', 'F861']),
                                                  np.ones(11),
                                                  age_delta)
            # cs.load_table('best5_100Z_test6.dat')[-5,0:11])tables
            np.set_printoptions(precision=2)
            # print (nomes[i], chi_age[0:3],chi_Z[0:3],chi_ebv[0:3])

            # Ajuste com aglomerados desavermelhados, usar tabela mag1.dat (mudar no cabeçalho da função)
            # chi_Z=np.log10((chi_Z[0]/(1-chi_Z[0]-(0.2485+(1.78*chi_Z[0]))))/(0.0207)),np.log10(chi_Z[1]),chi_Z[2]
            chi_Z = np.log10(chi_Z[0] / 0.019), np.log10(chi_Z[1]), chi_Z[2]
            sint_cluster_param[j, :] = chi_Z
            # print(nomes[i],chi_Z[0:3],np.abs(chi_Z[0]-literature[cluster_locator,0]),
            # np.abs(chi_Z[1]-literature[cluster_locator,1]))

            # Ajuste com aglomerados avermelhados, usar tabela mag.dat (mudar no cabeçalho da função)
            # chi_ebv=np.log10((chi_ebv[0]/(1-chi_ebv[0]-(0.2485+(1.78*chi_ebv[0]))))/(0.0207)),np.log10(chi_ebv[1]),chi_ebv[2]
            # sint_cluster_param[j,:]=chi_ebv
            # print(nomes[i], chi_ebv[0:3], np.abs(chi_ebv[0] - literature[cluster_locator, 0]),
            # np.abs(chi_ebv[1] - literature[cluster_locator, 1]))

        np.set_printoptions(precision=2)
        # percentils de idade e metalicidade
        Z_percentiles = [np.percentile(sint_cluster_param[:, 0], 50, interpolation='nearest'),
                         np.percentile(sint_cluster_param[:, 0], 16, interpolation='nearest'),
                         np.percentile(sint_cluster_param[:, 0], 84, interpolation='nearest'),
                         float(stats.mode(sint_cluster_param[:, 0])[0])]
        age_percentiles = [np.percentile(sint_cluster_param[:, 1], 50, interpolation='nearest'),
                           np.percentile(sint_cluster_param[:, 1], 16, interpolation='nearest'),
                           np.percentile(sint_cluster_param[:, 1], 84, interpolation='nearest'),
                           float(stats.mode(sint_cluster_param[:, 1])[0])]
        # print(Z_percentiles, age_percentiles)
        print("{} {} ({:.2f} {:.2f} {:.2f} {:.2f} {:.2f}) ({:.2f} {:.2f} {:.2f} {:.2f} {:.2f})"
              .format(nomes[i], fields[i],
                      literature_param[cluster_locator, 0], Z_percentiles[3], Z_percentiles[0], Z_percentiles[1],
                      Z_percentiles[2],
                      literature_param[cluster_locator, 1], age_percentiles[3], age_percentiles[0], age_percentiles[1],
                      age_percentiles[2]))

        if i == 0:
            with open('results/' + pasta + '/clusters_output.txt', 'w') as f:
                f.write(
                    '#Object field Lit_[Fe/H] Mode_[Fe/H] Median_[Fe/H] 16_[Fe/H] 84_[Fe/H] Lit_log(Age) Mode_log(Age) Median_log(Age) 16_log(Age) 84_log(Age) \n')
                f.write("{} {} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f}\n"
                        .format(nomes[i], fields[i],
                                literature_param[cluster_locator, 0], Z_percentiles[3], Z_percentiles[0],
                                Z_percentiles[1],
                                Z_percentiles[2],
                                literature_param[cluster_locator, 1], age_percentiles[3], age_percentiles[0],
                                age_percentiles[1], age_percentiles[2]))
                # f.write('{} {:.2f} {:.2f} {:.2f} {:.4f} {:.2f}\n'.format(nomes[i], chi_Z[0], chi_Z[1], chi_Z[2],
                # np.abs(chi_Z[0] - literature[i, 0]),
                # np.abs(chi_Z[1] - literature[i, 0])))
        else:
            with open('results/' + pasta + '/clusters_output.txt', 'a') as f:
                f.write("{} {} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f}\n"
                        .format(nomes[i], fields[i],
                                literature_param[cluster_locator, 0], Z_percentiles[3], Z_percentiles[0],
                                Z_percentiles[1],
                                Z_percentiles[2],
                                literature_param[cluster_locator, 1], age_percentiles[3], age_percentiles[0],
                                age_percentiles[1], age_percentiles[2]))
                # f.write('{} {:.2f} {:.2f} {:.2f} {:.4f} {:.2f}\n'.format(nomes[i], chi_Z[0], chi_Z[1], chi_Z[2],
                #                                                        np.abs(chi_Z[0] - literature[i, 0]),
                #                                                        np.abs(chi_Z[1] - literature[i, 1])))

        #######################################################################################################################
        #######################################################################################################################
        #######################################################################################################################
        #######################################################################################################################

        # plt.subplots(3,3, constrained_layout=True)#(figsize=(12, 12))
        # plt.suptitle(fields[i]+'_'+nomes[i])
        mmm_Z = np.mean((Z_percentiles[0], Z_percentiles[3]))
        mmm_age = np.mean((age_percentiles[0], age_percentiles[3]))
        fig = plt.figure(constrained_layout=False)
        gs = fig.add_gridspec(nrows=2, ncols=2, width_ratios=[2, 1], height_ratios=[1, 2])

        # plt.subplot2grid((3, 3), (1, 0), colspan=2, rowspan=2)

        # plt.subplot2grid((3, 3), (1, 2), rowspan=2)
        # plt.subplot(221)
        a = fig.add_subplot(gs[1, 1])
        n, bins_Z, patches = plt.hist(sint_cluster_param[:, 0], orientation='horizontal', color='khaki')
        # bins=np.arange((literature_param[cluster_locator,0]-1),
        #               (literature_param[cluster_locator,0]+1),0.1),
        # orientation='horizontal')
        # print(n, bins_Z)
        a.axes.get_yaxis().set_visible(False)
        plt.hlines(literature_param[cluster_locator, 0], 0, np.max(n), label='literature', color='green', zorder=5)
        plt.hlines(Z_percentiles[3], 0, np.max(n), label='mode', color='purple', zorder=3)
        plt.hlines(Z_percentiles[0], 0, np.max(n), label='median', color='red', zorder=2)
        # plt.ylim(np.min(bins_Z) - 0.1, np.max(bins_Z) + 0.1)
        plt.ylim(-2.5, 1.5)
        plt.hlines(mmm_Z, 0, np.max(n), label='mmm', color='blue', zorder=4)
        plt.hlines(Z_percentiles[1], 0, np.max(n), label='16%', color='black', linestyles='dashed', zorder=5)
        plt.hlines(Z_percentiles[2], 0, np.max(n), label='84%', color='black', linestyles='dashed', zorder=5)
        # plt.legend(bbox_to_anchor=(1.0,1.05), loc="upper left", frameon=False)
        # plt.ylabel('[Fe/H]',fontsize=15)
        plt.xlabel('N', fontsize=15)
        # plt.tick_p    arams(labelsize=15)

        # plt.subplot2grid((3, 3), (0, 0), colspan=2)
        # plt.subplot(222)
        a = fig.add_subplot(gs[0, 0])
        n, bins_age, patches = plt.hist(sint_cluster_param[:, 1], orientation='vertical', color='khaki')
        # bins=np.arange((literature_param[cluster_locator,1]-1),
        #                                    (literature_param[cluster_locator,1]+1),0.1),
        #                  orientation='vertical')
        # print(n,bins_age)

        plt.vlines(literature_param[cluster_locator, 1], 0, np.max(n), label='literature', color='green', zorder=5)
        plt.vlines(age_percentiles[3], 0, np.max(n), label='mode', color='purple', zorder=3)
        plt.vlines(age_percentiles[0], 0, np.max(n), label='median', color='red', zorder=2)
        a.axes.get_xaxis().set_visible(False)
        # plt.xlim(np.min(bins_age) - 0.05, np.max(bins_age) + 0.05)
        plt.xlim(6.9, 10.2)
        plt.vlines(mmm_age, 0, np.max(n), label='mmm', color='blue', zorder=4)
        plt.vlines(age_percentiles[1], 0, np.max(n), label='16%', color='black', linestyles='dashed', zorder=5)
        plt.vlines(age_percentiles[2], 0, np.max(n), label='84%', color='black', linestyles='dashed', zorder=5)
        plt.legend(bbox_to_anchor=(1.0, 1.05), loc="upper left", frameon=False)
        # plt.xlabel('log(Age)', fontsize=15)
        plt.ylabel('N', fontsize=15)
        # plt.tick_params(labelsize=15)

        fig.add_subplot(gs[1, 0])
        plt.ylabel('[Fe/H]', fontsize=15)
        plt.xlabel('log(Age)', fontsize=15)
        # plt.ylim(np.min(bins_Z) - 0.1, np.max(bins_Z) + 0.1)
        # plt.xlim(np.min(bins_age) - 0.05, np.max(bins_age) + 0.05)
        plt.xlim(6.9, 10.2)
        plt.ylim(-2.5, 1.5)
        d = np.vstack([sint_cluster_param[:, 1], sint_cluster_param[:, 0]])
        density = stats.gaussian_kde(d)(d)
        sc = plt.scatter(sint_cluster_param[:, 1], sint_cluster_param[:, 0], marker='.', c=density, s=100, edgecolor='',
                         cmap='inferno_r', zorder=4)
        plt.scatter(mmm_age, mmm_Z, marker='o', color='blue', zorder=5)
        # plt.hlines(mmm_Z, np.min(bins_age) - 0.05, np.max(bins_age) + 0.05, color='blue', zorder=3)
        # plt.vlines(mmm_age, np.min(bins_Z) - 0.1, np.max(bins_Z) + 0.1, color='blue', zorder=3)
        plt.scatter(literature_param[cluster_locator, 1], literature_param[cluster_locator, 0], marker='x',
                    color='green',
                    zorder=5)

        # fig.add_subplot(gs[:, 2])
        # cb = mpl.colorbar.ColorbarBase(cax, orientation='vertical')
        position = fig.add_axes([0.9, 0.02, 0.01, 0.9])
        cb = plt.colorbar(sc, cax=position, orientation='vertical', label='density of models fitted (gaussian_kde)')
        # plt.axis('off')
        fig.subplots_adjust(wspace=0, hspace=0)

        # plt.tight_layout()
        plt.savefig('results/' + pasta + '/mcmc_hists/' + fields[i] + '_' + nomes[i] + '.png',
                    dpi=300)  # ,bbox_inches='tight')
        plt.close()
        # plt.show()

    #######################################################################################################################
    #######################################################################################################################
    #######################################################################################################################
    #######################################################################################################################

    results_strings = np.loadtxt('results/' + pasta + '/clusters_output.txt', usecols=[0, 1], dtype=str)
    results_floats = np.loadtxt('results/' + pasta + '/clusters_output.txt', usecols=[2, 3, 4, 5, 6, 7, 8, 9, 10, 11],
                                dtype=float)
    results = np.hstack((results_strings, results_floats))

    results_latex = tabulate(results, tablefmt="latex",
                             headers=["Object", "field", "Lit_[Fe/H]", "Mode_[Fe/H]", "Mean_[Fe/H]" "-_[Fe/H]",
                                      "+_[Fe/H]", "Lit_log(Age)", "Mode_log(Age)", "Mean_log(Age)", "-_log(Age)",
                                      "+_log(Age)"])

    with open('results/' + pasta + '/clusters_output_latex.txt', 'w') as f:
        f.write("{}".format(results_latex))


# Calculando os parâmetros

mc_chisquare(pasta='all-noreg', age_delta=((0, 0)), Nsint=2000, Ns=3,
             phot='photometry_cat.txt', models='mag_new.dat')
pasta = 'all-noreg'
results_strings = np.loadtxt('results/' + pasta + '/clusters_output.txt', usecols=[0, 1], dtype=str)
age_lit = np.loadtxt('results/' + pasta + '/clusters_output.txt', usecols=[7], dtype=float)
age_mode = np.loadtxt('results/' + pasta + '/clusters_output.txt', usecols=[8], dtype=float)
age_median = np.loadtxt('results/' + pasta + '/clusters_output.txt', usecols=[9], dtype=float)

df = pd.read_table('results/' + pasta + '/clusters_output.txt', sep=' ')
df9 = df[df['Lit_log(Age)'] > 9]

reg_mmm = linear_model.LinearRegression()
reg_mmm.fit(df.loc[:, "Lit_log(Age)"].values.reshape((-1, 1)),
            df.loc[:, "Lit_log(Age)"] - (df.loc[:, "Mode_log(Age)"] + df.loc[:, "Median_log(Age)"]) / 2)
reg_mmm9 = linear_model.LinearRegression()
reg_mmm9.fit(df9.loc[:, "Lit_log(Age)"].values.reshape((-1, 1)),
             df9.loc[:, "Lit_log(Age)"] - (df9.loc[:, "Mode_log(Age)"] + df9.loc[:, "Median_log(Age)"]) / 2)
print(reg_mmm.coef_[0], reg_mmm.intercept_)
print(reg_mmm9.coef_[0], reg_mmm9.intercept_)

reg_median = linear_model.LinearRegression()
reg_median.fit(df.loc[:, "Lit_log(Age)"].values.reshape((-1, 1)),
               df.loc[:, "Lit_log(Age)"] - df.loc[:, "Median_log(Age)"])
reg_median9 = linear_model.LinearRegression()
reg_median9.fit(df9.loc[:, "Lit_log(Age)"].values.reshape((-1, 1)),
                df9.loc[:, "Lit_log(Age)"] - df9.loc[:, "Median_log(Age)"])
print(reg_median.coef_[0], reg_median.intercept_)
print(reg_median9.coef_[0], reg_median9.intercept_)

reg_mode = linear_model.LinearRegression()
reg_mode.fit(df.loc[:, "Lit_log(Age)"].values.reshape((-1, 1)),
             df.loc[:, "Lit_log(Age)"] - df.loc[:, "Mode_log(Age)"])
reg_mode9 = linear_model.LinearRegression()
reg_mode9.fit(df9.loc[:, "Lit_log(Age)"].values.reshape((-1, 1)),
              df9.loc[:, "Lit_log(Age)"] - df9.loc[:, "Mode_log(Age)"])
print(reg_mode.coef_[0], reg_mode.intercept_)
print(reg_mode9.coef_[0], reg_mode9.intercept_)

mc_chisquare(pasta='all-fullreg-mmm', age_delta=((reg_mmm.coef_[0], reg_mmm.intercept_)), Nsint=2000, Ns=1,
             phot='photometry_cat.txt', models='mag_new.dat')
mc_chisquare(pasta='all-fullreg-median', age_delta=((reg_median.coef_[0], reg_median.intercept_)), Nsint=2000, Ns=1,
             phot='photometry_cat.txt', models='mag_new.dat')
mc_chisquare(pasta='all-fullreg-mode', age_delta=((reg_mode.coef_[0], reg_mode.intercept_)), Nsint=2000, Ns=3,
             phot='photometry_cat.txt', models='mag_new.dat')

mc_chisquare(pasta='all-9reg-mmm', age_delta=((reg_mmm9.coef_[0], reg_mmm9.intercept_)), Nsint=2000, Ns=1,
             phot='photometry_cat.txt', models='mag_new.dat')
mc_chisquare(pasta='all-9reg-median', age_delta=((reg_median9.coef_[0], reg_median9.intercept_)), Nsint=2000, Ns=1,
             phot='photometry_cat.txt', models='mag_new.dat')
mc_chisquare(pasta='all-9reg-mode', age_delta=((reg_mode9.coef_[0], reg_mode9.intercept_)), Nsint=2000, Ns=3,
             phot='photometry_cat.txt', models='mag_new.dat')
