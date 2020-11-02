# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
from astropy.stats import biweight_scale, biweight_location
from tabulate import tabulate


def plot_errors(param, x, y, ymin=-3, ymax=3, label_x='label_x', label_y='label_y', file='test', ticks=0, est='mmm'):
    fig, ax = plt.subplots(ncols=1)
    ax.set_xlabel(label_x)
    ax.set_ylabel(label_y, fontsize=15)
    points = plt.scatter(x, y, marker='o', c=ebv_results, cmap='inferno_r')
    plt.colorbar(points, label='E(B-V)')

    if (param == 'Z'):
        location = biweight_location(y)
        scale = biweight_scale(y)
        if (ticks == 1):
            ax.hlines(location + scale, 0, len(results_floats[:, 1]), label='biweight scale=' + str(scale))
            ax.hlines(location - scale, 0, len(results_floats[:, 1]))
            ax.hlines(location, 0, len(results_floats[:, 1]), color='blue',
                      label='biweight location=' + str(location))
            ax.set_xticks(x)
            ax.set_xticklabels((results_strings[:, 0]), rotation='vertical', fontsize=8)
        else:
            ax.hlines(location + scale, np.min(results_floats[:, 0]), np.max(results_floats[:, 0]),
                      label='biweight scale=' + str(scale))
            ax.hlines(location - scale, np.min(results_floats[:, 0]), np.max(results_floats[:, 0]))
            ax.hlines(location, np.min(results_floats[:, 0]), np.max(results_floats[:, 0]), color='blue',
                      label='biweight location=' + str(location))

    if (param == 'age'):
        location = biweight_location(y)
        scale = biweight_scale(y)
        if (ticks == 1):
            ax.hlines(location + scale, 0, len(results_floats[:, 1]), label='biweight scale=' + str(scale))
            ax.hlines(location - scale, 0, len(results_floats[:, 1]))
            ax.hlines(location, 0, len(results_floats[:, 1]), color='blue', label='biweight location=' + str(location))
            ax.set_xticks(x)
            ax.set_xticklabels((results_strings[:, 0]), rotation='vertical', fontsize=8)
        else:
            ax.hlines(location + scale, np.min(results_floats[:, 5]), np.max(results_floats[:, 5]),
                      label='biweight scale=' + str(scale))
            ax.hlines(location - scale, np.min(results_floats[:, 5]), np.max(results_floats[:, 5]))
            ax.hlines(location, np.min(results_floats[:, 5]), np.max(results_floats[:, 5]), color='blue',
                      label='biweight location=' + str(location))

    plt.ylim(ymin, ymax)
    plt.legend()
    plt.tight_layout()
    plt.savefig('results/' + pasta + '/' + file + '.png', dpi=300, bbox_inches='tight')
    plt.close()


def plot_hist_errors(y_Z, y_age, label_Z='label_Z', label_age='label_age', file='test'):
    fig = plt.figure(constrained_layout=False)
    gs = fig.add_gridspec(nrows=1, ncols=2)  # , width_ratios=[2, 1], height_ratios=[1, 2])
    a = fig.add_subplot(gs[0, 0])
    n, bins_Z, patches = plt.hist(y_Z, orientation='vertical', color='khaki')
    plt.xlabel(label_Z)
    a = fig.add_subplot(gs[0, 1])
    n, bins_age, patches = plt.hist(y_age, orientation='vertical', color='khaki')
    plt.xlabel(label_age)
    '''
    if (est == 'mmm'):
        location = [biweight_location(results_floats[:, 0] - mmm_Z),
                    biweight_location(results_floats[:, 5] - mmm_age)]
        scale = [biweight_scale(results_floats[:, 0] - mmm_Z),
                 biweight_scale(results_floats[:, 5] - mmm_age)]
    if (est == 'mean'):
        location = [biweight_location(results_floats[:, 0] - results_floats[:, 2]),
                    biweight_location(results_floats[:, 5] - results_floats[:, 7])]
        scale = [biweight_scale(results_floats[:, 0] - results_floats[:, 2]),
                 biweight_scale(results_floats[:, 5] - results_floats[:, 7])]
    if (est == 'mode'):
        location = [biweight_location(results_floats[:, 0] - results_floats[:, 1]),
                    biweight_location(results_floats[:, 5] - results_floats[:, 6])]
        scale = [biweight_scale(results_floats[:, 0] - results_floats[:, 1]),
                 biweight_scale(results_floats[:, 5] - results_floats[:, 6])]
    '''
    # plt.legend()
    plt.tight_layout()
    plt.savefig('results/' + pasta + '/' + file + '.png', dpi=300, bbox_inches='tight')
    plt.close()


# pasta='all'
pasta = 'all-noreg'
# pasta = 'all-fullreg'
# pasta = 'all_reg9-10'
results_strings = np.loadtxt('results/' + pasta + '/clusters_output.txt', usecols=[0, 1], dtype=str)
results_floats = np.loadtxt('results/' + pasta + '/clusters_output.txt', usecols=[2, 3, 4, 5, 6, 7, 8, 9, 10, 11],
                            dtype=float)
param_table_name = np.loadtxt('results/parameter-cat.txt', usecols=[0], dtype=str)
param_table_ebv = np.loadtxt('results/parameter-cat.txt', usecols=[7], dtype=float)
results = np.hstack((results_strings, results_floats))

mmm_Z = np.zeros((len(results_floats)))
mmm_age = np.zeros((len(results_floats)))
ebv_results = np.zeros((len(results_floats)))
for i in range(len(results_floats)):
    cluster_locator = 0  # locating the cluster in the literature table
    while (results_strings[i, 0] != param_table_name[cluster_locator]):
        cluster_locator = cluster_locator + 1
    ebv_results[i] = param_table_ebv[cluster_locator]
    mmm_Z[i] = np.mean((results_floats[i, 1], results_floats[i, 2]))
    mmm_age[i] = np.mean((results_floats[i, 6], results_floats[i, 7]))

results_latex = tabulate(results, tablefmt="latex",
                         headers=["Object", "field", "Lit_[Fe/H]", "Mode_[Fe/H]", "Mean_[Fe/H]" "-_[Fe/H]",
                                  "+_[Fe/H]", "Lit_log(Age)", "Mode_log(Age)", "Mean_log(Age)", "-_log(Age)",
                                  "+_log(Age)"])

with open('results/' + pasta + '/clusters_output_latex.txt', 'w') as f:
    f.write("{}".format(results_latex))

# Estimator: MMM
plot_errors(x=np.arange(0, len(results_strings)), y=(results_floats[:, 0] - mmm_Z[:]),
            label_x='Clusters', label_y='ref_[Fe/H]-mmm_[Fe/H]', param='Z', ticks=1, file='clusters_Z_mmm')
plot_errors(x=np.arange(0, len(results_strings)), y=(results_floats[:, 5] - mmm_age[:]),
            label_x='Clusters', label_y='ref_log(age)-mmm_log(age)', param='age', ticks=1, file='clusters_age_mmm')
plot_errors(x=results_floats[:, 0], y=(results_floats[:, 0] - mmm_Z[:]),
            label_x='Clusters', label_y='ref_[Fe/H]-mmm_[Fe/H]', param='Z', ticks=0, file='clusters_Z-error_mmm')
plot_errors(x=results_floats[:, 5], y=(results_floats[:, 5] - mmm_age[:]),
            label_x='Clusters', label_y='ref_log(age)-mmm_log(age)', param='age', ticks=0,
            file='clusters_age-error_mmm')
plot_hist_errors(y_Z=(results_floats[:, 0] - mmm_Z[:]), y_age=(results_floats[:, 5] - mmm_age[:]),
                 label_Z='ref_[Fe/H]-mmm_[Fe/H]', label_age='ref_log(age)-mmm_log(age)', file='hist_mmm')

# Estimator: Mode
plot_errors(x=np.arange(0, len(results_strings)), y=(results_floats[:, 0] - results_floats[:, 1]),
            label_x='Clusters', label_y='ref_[Fe/H]-mmm_[Fe/H]', param='Z', ticks=1, file='clusters_Z_mode')
plot_errors(x=np.arange(0, len(results_strings)), y=(results_floats[:, 5] - results_floats[:, 6]),
            label_x='Clusters', label_y='ref_log(age)-mode_log(age)', param='age', ticks=1, file='clusters_age_mode')
plot_errors(x=results_floats[:, 0], y=(results_floats[:, 0] - results_floats[:, 1]),
            label_x='ref_[Fe/H]', label_y='ref_[Fe/H]-mode_[Fe/H]', param='Z', ticks=0, file='clusters_Z-error_mode')
plot_errors(x=results_floats[:, 5], y=(results_floats[:, 5] - results_floats[:, 6]),
            label_x='Clusters', label_y='ref_log(age)-mode_log(age)', param='age', ticks=0,
            file='clusters_age-error_mode')
plot_hist_errors(y_Z=(results_floats[:, 0] - results_floats[:, 1]), y_age=(results_floats[:, 5] - results_floats[:, 6]),
                 label_Z='ref_[Fe/H]-mode_[Fe/H]', label_age='ref_log(age)-mode_log(age)', file='hist_mode')

# Estimator: Median
plot_errors(x=np.arange(0, len(results_strings)), y=(results_floats[:, 0] - results_floats[:, 2]),
            label_x='Clusters', label_y='ref_[Fe/H]-median_[Fe/H]', param='Z', ticks=1, file='clusters_Z_median')
plot_errors(x=np.arange(0, len(results_strings)), y=(results_floats[:, 5] - results_floats[:, 7]),
            label_x='Clusters', label_y='ref_log(age)-median_log(age)', param='age', ticks=1,
            file='clusters_age_median')
plot_errors(x=results_floats[:, 0], y=(results_floats[:, 0] - results_floats[:, 2]),
            label_x='Clusters', label_y='ref_[Fe/H]-median_[Fe/H]', param='Z', ticks=0, file='clusters_Z-error_median')
plot_errors(x=results_floats[:, 5], y=(results_floats[:, 5] - results_floats[:, 7]),
            label_x='ref_log(age)', label_y='ref_log(age)-median_log(age)', param='age', ticks=0,
            file='clusters_age-error_median')
plot_hist_errors(y_Z=(results_floats[:, 0] - results_floats[:, 2]), y_age=(results_floats[:, 5] - results_floats[:, 7]),
                 label_Z='ref_[Fe/H]-median_[Fe/H]', label_age='ref_log(age)-median_log(age)', file='hist_median')
