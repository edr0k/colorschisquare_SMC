# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
from astropy.stats import biweight_midvariance

import chisquare_v4 as cs


def histOutline(dataIn, *args, **kwargs):
    (histIn, binsIn) = np.histogram(dataIn, *args, **kwargs)

    stepSize = binsIn[1] - binsIn[0]

    bins = np.zeros(len(binsIn) * 2 + 2, dtype=np.float)
    data = np.zeros(len(binsIn) * 2 + 2, dtype=np.float)
    for bb in range(len(binsIn)):
        bins[2 * bb + 1] = binsIn[bb]
        bins[2 * bb + 2] = binsIn[bb] + stepSize
        if bb < len(histIn):
            data[2 * bb + 1] = histIn[bb]
            data[2 * bb + 2] = histIn[bb]

    bins[0] = bins[1]
    bins[-1] = bins[-2]
    data[0] = 0
    data[-1] = 0

    return (bins, data)


def weight_gen(band, band_ref, file):
    color = ['#000000', '#FF0000', '#00FF00', '#0000FF', '#FFFF00', '#00FFFF', '#FF00FF', '#C0C0C0', '#800000',
             '#808000',
             '#008000']
    data = file
    plt.subplots(figsize=(13, 10))
    k = 1
    for j in range(0, len(band)):
        plt.subplot(2, 1 + len(band) / 2, k)
        k = k + 1
        for i in range(0, 5):
            gen = 0
            plt.plot(gen + 1, data[5 * gen + i, j], marker='.', color=color[i], label=i + 1, zorder=5 - i)
            for gen in range(1, len(data) / 5):
                plt.plot(gen + 1, data[5 * gen + i, j], marker='.', color=color[i], zorder=5 - i)
        plt.title(band[j] + '-' + band_ref[j])
        plt.xlabel('Generation')
        plt.ylabel('Weight value')

    plt.legend(loc='lower right', title='Rank')
    plt.tight_layout()
    # print ('image generated')
    plt.savefig('weitgh_gen.png')


# def bargraph_weight():

def histogram(band, band_ref, best_weights, training_data, aeg_models, chisquare_models):
    weights = best_weights[-5, 0:len(band)]
    print(weights)
    Z_bef = np.zeros((len(training_data)))
    age_bef = np.zeros((len(training_data)))
    ebv_bef = np.zeros((len(training_data)))
    Z_aft = np.zeros((len(training_data)))
    age_aft = np.zeros((len(training_data)))
    ebv_aft = np.zeros((len(training_data)))

    plt.subplots(figsize=(15, 7))
    for j in range(0, len(training_data)):
        chi = cs.fit(chisquare_models, training_data, j, band, band_ref, np.ones((len(band))))
        Z_bef[j] = chi[0] - aeg_models[j, 0]
        age_bef[j] = (chi[1] - aeg_models[j, 1]) / 10 ** 9
        ebv_bef[j] = chi[2] - aeg_models[j, 2]

    for j in range(0, len(training_data)):
        chi = cs.fit(chisquare_models, training_data, j, band, band_ref, np.array([weights]))
        Z_aft[j] = chi[0] - aeg_models[j, 0]
        age_aft[j] = (chi[1] - aeg_models[j, 1]) / 10 ** 9
        ebv_aft[j] = chi[2] - aeg_models[j, 2]

    plt.subplot(1, 3, 1)
    plt.hist(Z_bef, bins=np.arange(-0.055, 0.055, 0.01), color='#aaaa22', rwidth=0.95)
    bins_Z_aft, Z_aft_data = histOutline(Z_aft, bins=np.arange(-0.055, 0.055, 0.01))
    plt.plot(bins_Z_aft, Z_aft_data, '-', color="#1B6AC6", lw=4)
    plt.xlabel(r'$\Delta$Z' + '\n' + 'Before training: ' + r'$\mu={0:1.6f},\ \zeta={1:1.6f}$'.format(np.mean(Z_bef),
                                                                                                     biweight_midvariance(
                                                                                                         Z_bef))
               + '\n' + 'After training: ' + r'$\mu={0:1.6f},\ \zeta={1:1.6f}$'.format(np.mean(Z_aft),
                                                                                       biweight_midvariance(Z_aft)))
    plt.xlim(-0.05, 0.05)

    plt.subplot(1, 3, 2)
    plt.hist(age_bef, bins=np.arange(-1.05, 1.05, 0.1), color='#aaaa22', rwidth=0.95)
    bins_age_aft, age_aft_data = histOutline(age_aft, bins=np.arange(-1.05, 1.05, 0.1))
    plt.plot(bins_age_aft, age_aft_data, '-', color="#1B6AC6", lw=4)
    plt.xlabel(r'$\Delta$age/$10^{9}$' + '\n' + 'Before training: ' + r'$\mu={0:1.6f},\ \zeta={1:1.6f}$'.format(
        np.mean(age_bef), biweight_midvariance(age_bef))
               + '\n' + 'After training: ' + r'$\mu={0:1.6f},\ \zeta={1:1.6f}$'.format(np.mean(age_aft),
                                                                                       biweight_midvariance(age_aft)))
    plt.xlim(-0.5, 0.5)

    plt.subplot(1, 3, 3)
    plt.hist(ebv_bef, bins=np.arange(-0.45, 0.45, 0.05), label='Before training', color='#aaaa22', rwidth=0.95)
    bins_ebv_aft, ebv_aft_data = histOutline(ebv_aft, bins=np.arange(-0.45, 0.45, 0.05))
    plt.plot(bins_ebv_aft, ebv_aft_data, '-', color="#1B6AC6", lw=4, label='After training')
    plt.xlabel(
        '$\Delta$E(B-V)' + '\n' + 'Before training: ' + r'$\mu={0:1.6f},\ \zeta={1:1.6f}$'.format(np.mean(ebv_bef),
                                                                                                  biweight_midvariance(
                                                                                                      ebv_bef))
        + '\n' + 'After training: ' + r'$\mu={0:1.6f}, \zeta={1:1.6f}$'.format(np.mean(ebv_aft),
                                                                               biweight_midvariance(ebv_aft)))
    # plt.xlim(-5,2)
    plt.legend()

    # plt.suptitle('Histograms')
    plt.tight_layout()
    # plt.savefig('histogram_Z&Age.png')
    plt.show()


def histogram2(band, band_ref, best_weights, training_data, aeg_models, chisquare_models):
    weights = best_weights[-5, 0:len(band)]
    Z_bef = np.zeros((len(training_data)))
    age_bef = np.zeros((len(training_data)))
    ebv_bef = np.zeros((len(training_data)))
    Z_aft = np.zeros((len(training_data)))
    age_aft = np.zeros((len(training_data)))
    ebv_aft = np.zeros((len(training_data)))

    plt.subplots(figsize=(15, 7))
    for j in range(0, len(training_data)):
        # chi = cs.fit(chisquare_models, training_data, j, band, band_ref, np.ones((len(band))))
        chi2, chi1 = cs.full_fit(cs.load_table('chisquare_models_test.dat'), cs.load_table('training_data.dat'), j,
                                 (['uJava', 'uJava', 'uJava', 'uJava', 'F378', 'F378', 'F430', 'gSDSS', 'F660', 'F861',
                                   'uJava']),
                                 (['F660', 'rSDSS', 'iSDSS', 'zSDSS', 'F515', 'gSDSS', 'gSDSS', 'zSDSS', 'gSDSS',
                                   'uJava', 'F515']),
                                 cs.load_table('best5_100age.dat')[-5, 0:11],
                                 # fixar os pesos referentes as cores da idade antes e depois do treino para só ver a
                                 # influência dos pesos das cores da metalicidade
                                 # np.ones(11),
                                 (['iSDSS', 'rSDSS', 'uJava', 'F378', 'F410', 'F430', 'F515', 'F660', 'F660', 'gSDSS',
                                   'iSDSS']),
                                 (['zSDSS', 'zSDSS', 'zSDSS', 'zSDSS', 'zSDSS', 'iSDSS', 'zSDSS', 'F861', 'iSDSS',
                                   'zSDSS', 'F861']),
                                 np.ones(11))
        # cs.load_table('best5_100Z_test2.dat')[-5, 0:11])
        Z_bef[j] = chi2[0] - aeg_models[j, 0]
        age_bef[j] = (chi2[1] - aeg_models[j, 1]) / 10 ** 9
        ebv_bef[j] = chi2[2] - aeg_models[j, 2]

        # chi = cs.fit(chisquare_models, training_data, j, band, band_ref, np.array([weights]))
        chi2, chi1 = cs.full_fit(cs.load_table('chisquare_models_test.dat'), cs.load_table('training_data.dat'), j,
                                 (['uJava', 'uJava', 'uJava', 'uJava', 'F378', 'F378', 'F430', 'gSDSS', 'F660', 'F861',
                                   'uJava']),
                                 (['F660', 'rSDSS', 'iSDSS', 'zSDSS', 'F515', 'gSDSS', 'gSDSS', 'zSDSS', 'gSDSS',
                                   'uJava', 'F515']),
                                 cs.load_table('best5_100age.dat')[-5, 0:11],
                                 # np.ones(11),
                                 (['iSDSS', 'rSDSS', 'uJava', 'F378', 'F410', 'F430', 'F515', 'F660', 'F660', 'gSDSS',
                                   'iSDSS']),
                                 (['zSDSS', 'zSDSS', 'zSDSS', 'zSDSS', 'zSDSS', 'iSDSS', 'zSDSS', 'F861', 'iSDSS',
                                   'zSDSS', 'F861']),
                                 # np.ones(11))
                                 cs.load_table('best5_100Z_test4.dat')[-5, 0:11])
        Z_aft[j] = chi2[0] - aeg_models[j, 0]
        age_aft[j] = (chi2[1] - aeg_models[j, 1]) / 10 ** 9
        ebv_aft[j] = chi2[2] - aeg_models[j, 2]

    # N=20
    plt.subplot(1, 3, 1)
    plt.hist(Z_bef, bins=np.arange(-0.055, 0.055, 0.01), color='#aaaa22', rwidth=0.95)
    bins_Z_aft, Z_aft_data = histOutline(Z_aft, bins=np.arange(-0.055, 0.055, 0.01))
    plt.plot(bins_Z_aft, Z_aft_data, '-', color="#1B6AC6", lw=4)
    # plt.hist(Z_aft, bins=np.arange(-0.055, 0.055, 0.01), color='#aaab22', rwidth=0.95)
    plt.xlabel(r'$\Delta$Z' + '\n' + 'Before training: ' + r'$\mu={0:1.6f},\ \zeta={1:1.6f}$'.format(np.mean(Z_bef),
                                                                                                     biweight_midvariance(
                                                                                                         Z_bef))
               + '\n' + 'After training: ' + r'$\mu={0:1.6f},\ \zeta={1:1.6f}$'.format(np.mean(Z_aft),
                                                                                       biweight_midvariance(Z_aft)))
    plt.xlim(-0.05, 0.05)

    plt.subplot(1, 3, 2)
    plt.hist(age_bef, bins=np.arange(-1.05, 1.05, 0.1), color='#aaaa22', rwidth=0.95)
    bins_age_aft, age_aft_data = histOutline(age_aft, bins=np.arange(-1.05, 1.05, 0.1))
    plt.plot(bins_age_aft, age_aft_data, '-', color="#1B6AC6", lw=4)
    # plt.hist(age_aft, bins=np.arange(-1.05, 1.05, 0.1), color='#aaaa22', rwidth=0.95)
    plt.xlabel(r'$\Delta$age/$10^{9}$' + '\n' + 'Before training: ' + r'$\mu={0:1.6f},\ \zeta={1:1.6f}$'.format(
        np.mean(age_bef), biweight_midvariance(age_bef))
               + '\n' + 'After training: ' + r'$\mu={0:1.6f},\ \zeta={1:1.6f}$'.format(np.mean(age_aft),
                                                                                       biweight_midvariance(age_aft)))
    plt.xlim(-0.5, 0.5)

    plt.subplot(1, 3, 3)
    plt.hist(ebv_bef, bins=np.arange(-0.45, 0.45, 0.05), label='Before training', color='#aaaa22', rwidth=0.95)
    bins_ebv_aft, ebv_aft_data = histOutline(ebv_aft, bins=np.arange(-0.45, 0.45, 0.05))
    plt.plot(bins_ebv_aft, ebv_aft_data, '-', color="#1B6AC6", lw=4, label='After training')
    # plt.hist(ebv_aft, bins=np.arange(-0.4, 0.4, 0.05), label='After training', color='#aaaa22', rwidth=0.95)
    plt.xlabel(
        '$\Delta$E(B-V)' + '\n' + 'Before training: ' + r'$\mu={0:1.6f},\ \zeta={1:1.6f}$'.format(np.mean(ebv_bef),
                                                                                                  biweight_midvariance(
                                                                                                      ebv_bef))
        + '\n' + 'After training: ' + r'$\mu={0:1.6f}, \zeta={1:1.6f}$'.format(np.mean(ebv_aft),
                                                                               biweight_midvariance(ebv_aft)))
    # plt.xlim(-5,2)
    plt.legend()

    # plt.suptitle('Histograms')
    plt.tight_layout()
    # plt.savefig('histogram_Z&Age.png')
    plt.show()
