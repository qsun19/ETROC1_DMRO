#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import time
import numpy as np
from scipy.stats import norm
from collections import Counter
from scipy.fftpack import fft
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from matplotlib import rcParams
from matplotlib import gridspec
'''
TDC Converted Data Plot
@author: Wei Zhang
@date: January 3, 2020
@address: SMU
'''
#===================================================================================#
## plot parameters
lw_grid = 0.5                   # grid linewidth
fig_dpi = 800                   # save figure's resolution
lw = 0.8

#===================================================================================#
def Code_Distribution():
    file_list = []
    for root, dirs, files in os.walk(".", topdown=False):
        for name in files:
            print(os.path.join(root, name))
            file_list += [os.path.join(root, name)]
    print(file_list)

    PhaseAdj_vs_TOA = [[] for x in range(3)]
    for i in range(len(file_list)-3):
        PhaseAdj = int(file_list[i].split("=")[8].split('_')[0])
        TOA_Code = []
        TOT_Code = []
        Cal_Code = []
        HitFlag = []

        with open(file_list[i], 'r') as outfile:
            j = 0
            for line in outfile.readlines():
                if int(line.split()[3]) == 1:
                    TOA_Code += [int(line.split()[0])]
                    TOT_Code += [int(line.split()[1])]
                    Cal_Code += [int(line.split()[2])]
                    HitFlag += [int(line.split()[3])]
                j += 1
        if len(TOA_Code) != 0:
            (mu, sigma) = norm.fit(TOA_Code)
            PhaseAdj_vs_TOA[0] += [PhaseAdj]
            PhaseAdj_vs_TOA[1] += [mu]
            PhaseAdj_vs_TOA[2] += [sigma]

    coef = np.polyfit(PhaseAdj_vs_TOA[0], PhaseAdj_vs_TOA[1], 1)
    print(coef)
    poly1d_fn = np.poly1d(coef)

    print(PhaseAdj_vs_TOA)
    plt.figure(figsize=(7,5))
    plt.errorbar(PhaseAdj_vs_TOA[0], PhaseAdj_vs_TOA[1], PhaseAdj_vs_TOA[2], color='g', elinewidth=0.8, linewidth=0.5, fmt='--o', ecolor='b', capthick=0.5, capsize=2, markersize=0.2, marker='.', label='TOA_Code Mean')
    # plt.plot(PhaseAdj_vs_TOA[0], poly1d_fn(PhaseAdj_vs_TOA[0], '--k')
    plt.title("PhaseAdj vs TOA", family="Times New Roman", fontsize=12)
    plt.xlabel("PhaseAdj Code", family="Times New Roman", fontsize=10)
    plt.ylabel("TOA_Code Mean", family="Times New Roman", fontsize=10)

    # plt.xscale('log')
    plt.xticks(family="Times New Roman", fontsize=8)
    plt.yticks(family="Times New Roman", fontsize=8)
    plt.grid(linestyle='-.', linewidth=lw_grid)
    plt.legend(fontsize=8, edgecolor='green')
    plt.savefig("TOA_vs_PhaseAdj.png", dpi=fig_dpi)         # save figure
    plt.clf()
#===================================================================================#
def main():
    print("OK!")
    # TDC_Converted_Data_Plot()
    Code_Distribution()
#===================================================================================#
if __name__ == '__main__':
    main()
