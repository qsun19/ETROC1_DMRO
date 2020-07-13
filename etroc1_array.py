import os, sys
import h5py
import numpy as np
import pandas as pd
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy import optimize
import matplotlib.pyplot as plt
from scipy import integrate
from scipy.optimize import curve_fit
from scipy.stats import norm
import matplotlib.font_manager as font_manager
from tqdm.autonotebook import tqdm
from matplotlib.backends.backend_pdf import PdfPages
import time
from scipy.optimize import curve_fit
from scipy.stats import norm
import scipy.stats

def plot_distribution_toa(list_in, num_bins= 20, range_default = None, xaxis = 'Time Resolution(ns)',
                        ylable = 'Occurrence', title = 'r$\delta$', pic = True, pdf = False):
    entries = len(list_in)
    mu, std = norm.fit(list_in)
    fig, ax4= plt.subplots(dpi=200)
    n,bins,patches=ax4.hist(list_in, bins=num_bins, range=range_default, density=False, 
                            label = 'entries: %d\nstd:%f\nmean:%f'%(entries, std, mu))
    xmin = np.min(list_in)
    xmax = np.max(list_in)
    x = np.linspace(xmin, xmax, 100)
    p = norm.pdf(x, mu, std)
    ax4.plot(x, np.max(n)*p/np.max(p), 'g', linewidth=2)
    ax4.grid()
    ax4.legend()
    ax4.set(xlabel=xaxis, ylabel=ylable,
               title=title)
    if pdf==True:
        pp.savefig(fig)
    if pic==True:
        plt.show()
    plt.close(fig)