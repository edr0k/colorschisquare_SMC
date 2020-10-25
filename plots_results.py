# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
from astropy.stats import biweight_scale, biweight_location
from tabulate import tabulate

results_strings = np.loadtxt('results/Z-age/clusters_output.txt', usecols=[0, 1], dtype=str)
results_floats = np.loadtxt('results/Z-age/clusters_output.txt', usecols=[2, 3, 4, 5, 6, 7, 8, 9, 10, 11], dtype=float)
results = np.hstack((results_strings, results_floats))
# erros=[results_floats[i,5] - results_floats[i,1], results_floats[i,5] - results_floats[i,6]]
location = [biweight_location(results_floats[:, 0] - results_floats[:, 1]),
            biweight_location(results_floats[:, 5] - results_floats[:, 6])]
scale = [biweight_scale(results_floats[:, 0] - results_floats[:, 1]),
         biweight_scale(results_floats[:, 5] - results_floats[:, 6])]
# bicor=[biweight_midcorrelation(results_floats[:,1]),biweight_midcorrelation(results_floats[:,6])]
print(location, scale)
results_latex = tabulate(results, tablefmt="latex",
                         headers=["Object", "field", "Lit_[Fe/H]", "Mode_[Fe/H]", "Mean_[Fe/H]" "-_[Fe/H]",
                                  "+_[Fe/H]", "Lit_log(Age)", "Mode_log(Age)", "Mean_log(Age)", "-_log(Age)",
                                  "+_log(Age)"])

x = np.arange(0, len(results_strings))
ax0.set_xlabel('Clusters')
ax0.set_ylabel('ref_[Fe/H]-calc_[Fe/H]', fontsize=15)
ax0.set_xticks(x)
ax0.set_xticklabels((results_strings[:, 0]), rotation='vertical', fontsize=8)
# plt.tick_params(labelsize=15)

for i in range(len(results_strings)):
    # print(i,np.abs(results_floats[i, 0] - results_floats[i, 1]))

    ax0.scatter(results_strings[i, 0], (results_floats[i, 0] - results_floats[i, 1]), marker='o')
    # if(results_floats[i,2]!=results_floats[i,3]):
    #    ax0.vlines(results_strings[i,0],results_floats[i,2],results_floats[i,3])

ax0.hlines(location[0] + scale[0], 0, len(results_floats[:, 1]), label='biweight scale')
ax0.hlines(location[0] - scale[0], 0, len(results_floats[:, 1]))
ax0.hlines(location[0], 0, len(results_floats[:, 1]), color='blue', label='biweight location')
plt.legend()
plt.tight_layout()
plt.savefig('results/Z-age/clusters_Z.png', dpi=300, bbox_inches='tight')

fig1, ax1 = plt.subplots(ncols=1)
# plt.subplots(121)
x = np.arange(0, len(results_strings))
ax1.set_xlabel('Clusters')
ax1.set_ylabel('ref_log(age)-calc_log(age)', fontsize=15)
ax1.set_xticks(x)
ax1.set_xticklabels(results_strings[:, 0], rotation='vertical', fontsize=8)  # adicionar campo no tick também?

for i in range(len(results_strings)):
    ax1.scatter(i, (results_floats[i, 5] - results_floats[i, 6]), marker='o')
    # plt.xticks(results_strings[i,1],rotation=55)
    # if(results_floats[i,2]!=results_floats[i,3]):
    #    ax1.vlines(results_strings[i,0],results_floats[i,2],results_floats[i,3])

ax1.hlines(location[1] + scale[1], 0, len(results_floats[:, 6]), label='biweight scale')
ax1.hlines(location[1] - scale[1], 0, len(results_floats[:, 6]))
ax1.hlines(location[1], 0, len(results_floats[:, 6]), color='blue', label='biweight location')

plt.legend()
plt.tight_layout()
plt.savefig('results/Z-age/clusters_age.png', dpi=300, bbox_inches='tight')
# plt.show()

with open('results/Z-age/clusters_output_latex.txt', 'w') as f:
    f.write("{}".format(results_latex))

fig0, ax0 = plt.subplots(ncols=1)
# x=np.arange(0,len(results_strings))
ax0.set_xlabel('ref_[Fe/H]')
ax0.set_ylabel('ref_[Fe/H]-calc_[Fe/H]', fontsize=15)
# ax0.set_xticks(x)
# ax0.set_xticklabels((results_strings[:,0]),rotation='vertical',fontsize=8)
# plt.tick_params(labelsize=15)

for i in range(len(results_strings)):
    # print(i,np.abs(results_floats[i, 0] - results_floats[i, 1]))

    ax0.scatter(results_floats[i, 0], (results_floats[i, 0] - results_floats[i, 1]), marker='o')
    # if(results_floats[i,2]!=results_floats[i,3]):
    #    ax0.vlines(results_strings[i,0],results_floats[i,2],results_floats[i,3])

ax0.hlines(location[0] + scale[0], np.min(results_floats[:, 0]), np.max(results_floats[:, 0]), label='biweight scale')
ax0.hlines(location[0] - scale[0], np.min(results_floats[:, 0]), np.max(results_floats[:, 0]))
ax0.hlines(location[0], np.min(results_floats[:, 0]), np.max(results_floats[:, 0]), color='blue',
           label='biweight location')
plt.legend()
plt.tight_layout()
plt.savefig('results/Z-age/clusters_Z-error.png', dpi=300, bbox_inches='tight')

fig1, ax1 = plt.subplots(ncols=1)
# plt.subplots(121)
# x=np.arange(0,len(results_strings))
ax1.set_xlabel('ref_log(age)')
ax1.set_ylabel('ref_log(age)-calc_log(age)', fontsize=15)
# ax1.set_xticks(x)
# ax1.set_xticklabels(results_strings[:,0],rotation='vertical',fontsize=8) #adicionar campo no tick também?


for i in range(len(results_strings)):
    ax1.scatter(results_floats[i, 5], (results_floats[i, 5] - results_floats[i, 7]), marker='o')
    # plt.xticks(results_strings[i,1],rotation=55)
    # if(results_floats[i,2]!=results_floats[i,3]):
    #    ax1.vlines(results_strings[i,0],results_floats[i,2],results_floats[i,3])

ax1.hlines(location[1] + scale[1], np.min(results_floats[:, 5]), np.max(results_floats[:, 5]), label='biweight scale')
ax1.hlines(location[1] - scale[1], np.min(results_floats[:, 5]), np.max(results_floats[:, 5]))
ax1.hlines(location[1], np.min(results_floats[:, 5]), np.max(results_floats[:, 5]), color='blue',
           label='biweight location')
plt.xlim(9, 10)
plt.ylim(-1, 1)
plt.legend()
plt.tight_layout()
plt.savefig('results/Z-age/clusters_age-error.png', dpi=300, bbox_inches='tight')
# plt.show()
