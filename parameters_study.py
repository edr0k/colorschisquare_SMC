# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np


def Z_study():
    table = np.loadtxt('tables/mag.dat')
    print(table)
    array_band = (['uJava', 'F378', 'F395', 'F410', 'F430', 'gSDSS', 'F515', 'rSDSS', 'F660', 'iSDSS', 'F861', 'zSDSS'])
    array_Z = np.array([0.0001, 0.0002, 0.0005, 0.0010, 0.0015, 0.0020, 0.0030, 0.0050,
                        0.0080, 0.0100, 0.0152, 0.0200, 0.0300, 0.0400, 0.0500])
    # calcular variação de cor com a metalicidade para ebv=0
    # fazer um plot para idades diferentes, de g-f515, u-z, z-f861, u-f378, u-z em função da metalicidade.
    for band1 in range(0, 11):
        for band2 in range(band1 + 1, 12):
            color = np.zeros((10, 16))
            for j in range(0, 15):  # Loop metalicidades
                for i in range(0, 70, 7):  # loop idade
                    if j == 0:
                        color[i // 7, 0] = table[21 * i + 1491 * j, 1]  # idade
                    color[i // 7, j + 1] = table[21 * i + 1491 * j, 3 + band1] - table[21 * i + 1491 * j, 3 + band2]
            print(color)
            cmap = plt.get_cmap('nipy_spectral')
            for i in range(0, 10):
                plt.plot(np.log10((array_Z / (1 - array_Z - (0.2485 + (1.78 * array_Z)))) / (0.0207)),
                         color[i, 1:], label='{:.2f}'.format(np.log10(color[i, 0])), color=cmap(25 * i), marker='.')
            plt.ylabel(array_band[band1] + '-' + array_band[band2], fontsize=15)
            plt.xlabel('[Fe/H]', fontsize=15)
            # plt.title('('+array_band[band1]+'-'+array_band[band2]+') x [Fe/H]')
            plt.legend(title='log(age)', bbox_to_anchor=(1.0, 1.05), loc="upper left", frameon=False)
            plt.tick_params(labelsize=15)
            plt.savefig('parameters_study/Z/' + array_band[band1] + '-' + array_band[band2] + '.jpg', dpi=300,
                        bbox_inches='tight')
            plt.close()


def age_study():
    table = np.loadtxt('tables/mag.dat')
    array_band = (['uJava', 'F378', 'F395', 'F410', 'F430', 'gSDSS', 'F515', 'rSDSS', 'F660', 'iSDSS', 'F861', 'zSDSS'])
    array_age = np.arange(6.6, 10.13, 0.05)
    # calcular variação de cor com a metalicidade para ebv=0
    # fazer um plot para idades diferentes, de g-f515, u-z, z-f861, u-f378, u-z em função da metalicidade.
    for band1 in range(0, 11):  # loop idade
        for band2 in range(band1 + 1, 12):
            color = np.zeros((15, 72))
            for i in range(0, 71):
                for j in range(0, 15):  # Loop metalicidades
                    color[j, 0] = table[21 * i + 1491 * j, 0]
                    color[j, i + 1] = table[21 * i + 1491 * j, 3 + band1] - table[21 * i + 1491 * j, 3 + band2]

            cmap = plt.get_cmap('nipy_spectral')
            for i in range(0, 15):
                plt.plot(array_age, color[i, 1:],
                         label="{:.2f}".format(
                             np.log10((color[i, 0] / (1 - color[i, 0] - (0.2485 + (1.78 * color[i, 0])))) / (0.0207))),
                         color=cmap(18 * i))
            plt.ylabel(array_band[band1] + '-' + array_band[band2], fontsize=15)
            plt.xlabel('log(age)', fontsize=15)
            # plt.title('('+array_band[band1]+'-'+array_band[band2]+') x log(age)')
            plt.tick_params(labelsize=15)
            plt.legend(title='[Fe/H]', bbox_to_anchor=(1.0, 1.05), loc="upper left", frameon=False)
            plt.savefig('parameters_study/age/' + array_band[band1] + '-' + array_band[band2] + '.jpg', dpi=300,
                        bbox_inches='tight')
            plt.close()


def ebv_study_age_var():
    table = np.loadtxt('tables/mag.dat')
    array_band = (['uJava', 'F378', 'F395', 'F410', 'F430', 'gSDSS', 'F515', 'rSDSS', 'F660', 'iSDSS', 'F861', 'zSDSS'])
    array_ebv = np.arange(0, 0.2, 0.01)

    for band1 in range(0, 11):
        for band2 in range(band1 + 1, 12):
            color = np.zeros((10, 21))
            for j in range(0, 20):  # Loop avermlhamento
                for i in range(0, 70, 7):  # loop idade
                    if j == 0:
                        color[i // 7, 0] = table[21 * i + j, 1]  # idade
                    color[i // 7, j + 1] = table[21 * i + j, 3 + band1] - table[21 * i + j, 3 + band2]
            print(color)
            if band1 == 0:
                return
            cmap = plt.get_cmap('nipy_spectral')
            for i in range(0, 10):
                plt.plot(array_ebv,
                         color[i, 1:], label='{:.2f}'.format(np.log10(color[i, 0])), color=cmap(25 * i), marker='.')
            plt.ylabel(array_band[band1] + '-' + array_band[band2], fontsize=15)
            plt.xlabel('E(B-V)', fontsize=15)
            # plt.title('('+array_band[band1]+'-'+array_band[band2]+') x [Fe/H]')
            plt.legend(title='log(age)', bbox_to_anchor=(1.0, 1.05), loc="upper left", frameon=False)
            plt.tick_params(labelsize=15)
            plt.savefig('parameters_study/E(B-V)/age_var/' + array_band[band1] + '-' + array_band[band2] + '.jpg',
                        dpi=300, bbox_inches='tight')
            plt.close()


def ebv_study_Z_var():
    table = np.loadtxt('tables/mag.dat')
    array_band = (['uJava', 'F378', 'F395', 'F410', 'F430', 'gSDSS', 'F515', 'rSDSS', 'F660', 'iSDSS', 'F861', 'zSDSS'])
    array_ebv = np.arange(0, 0.2, 0.01)

    for band1 in range(0, 11):  # loop idade
        for band2 in range(band1 + 1, 12):
            color = np.zeros((15, 21))
            for i in range(0, 20):
                for j in range(0, 15):  # Loop metalicidades
                    if i == 0:
                        color[j, 0] = table[i + 1491 * j, 0]
                    color[j, i + 1] = table[i + 1491 * j, 3 + band1] - table[i + 1491 * j, 3 + band2]
            print(color)
            if band1 == 0:
                return

            cmap = plt.get_cmap('nipy_spectral')
            for i in range(0, 15):
                plt.plot(array_ebv, color[i, 1:],
                         label="{:.2f}".format(
                             np.log10((color[i, 0] / (1 - color[i, 0] - (0.2485 + (1.78 * color[i, 0])))) / (0.0207))),
                         color=cmap(18 * i), marker='.')
            plt.ylabel(array_band[band1] + '-' + array_band[band2], fontsize=15)
            plt.xlabel('E(B-V)', fontsize=15)
            # plt.title('('+array_band[band1]+'-'+array_band[band2]+') x log(age)')
            plt.tick_params(labelsize=15)
            plt.legend(title='[Fe/H]', bbox_to_anchor=(1.0, 1.05), loc="upper left", frameon=False)
            plt.savefig('parameters_study/E(B-V)/Z_var/' + array_band[band1] + '-' + array_band[band2] + '.jpg',
                        dpi=300, bbox_inches='tight')
            plt.close()


Z_study()
age_study()
ebv_study_age_var()
# print('sep')
ebv_study_Z_var()
