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

def gaus(x,a,x0,sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

def gaus_fit(to_fit, num_bins=50):
    mean = np.mean(to_fit)
    sigma = np.std(to_fit, ddof=1)
    bins, edges = np.histogram(to_fit, num_bins, density=False)
    centers = 0.5*(edges[1:] + edges[:-1])
    popt, pcov = curve_fit(gaus,centers,bins,p0=[1,mean,sigma])
    mean = popt[1]
    sigma = popt[2]
    return mean,sigma, popt, pcov

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
 

def plot_distribution_time(list_in, file_item, num_bins= 20, range_default = None, xaxis = 'Time Resolution(ns)',
                        ylable = 'Occurrence', title = 'r$\delta$', pic = False, pdf = True):
    entries = len(list_in)
    mu, std = norm.fit(list_in)
    fig, ax4= plt.subplots(dpi=200)
    n,bins,patches=ax4.hist(list_in, bins=num_bins, range=range_default, density=False, 
                            label = 'entries: %d\nstd:%f ns\nmean:%f ns'%(entries, std, mu))
    xmin = np.min(list_in)
    xmax = np.max(list_in)
    x = np.linspace(xmin, xmax, 100)
    p = norm.pdf(x, mu, std)
    ax4.plot(x, np.max(n)*p/np.max(p), 'g', linewidth=2)
    ax4.grid()
    ax4.legend()
    ax4.set(xlabel=xaxis, ylabel='Occurrence',
               title=title)
    if pdf==True:
        pp.savefig(fig)
    if pic==True:
        plt.show()
    plt.close(fig)
    
    
def plot_distribution_time_nofit(list_in, file_item, num_bins= 20, range_default = None, xaxis = 'Time Resolution(ns)',
                        ylable = 'Occurrence', title = 'r$\delta$', pic = False, pdf = True):
    entries = len(list_in)
    mu = np.mean(list_in)
    std = np.std(list_in)
    fig, ax4= plt.subplots(dpi=200)
    n,bins,patches=ax4.hist(list_in, bins=num_bins, range=range_default, density=False, 
                            label = 'entries: %d\nstd:%f ps\nmean:%f ps'%(entries, std, mu))
    xmin = np.min(list_in)
    xmax = np.max(list_in)
    x = np.linspace(xmin, xmax, 100)
    p = norm.pdf(x, mu, std)
    #ax4.plot(x, np.max(n)*p/np.max(p), 'g', linewidth=2)
    ax4.grid()
    ax4.legend()
    ax4.set(xlabel=xaxis, ylabel='Occurrence',
               title=title)
    if pdf==True:
        pp.savefig(fig)
    if pic==True:
        plt.show()
    plt.close(fig)


def plot_distribution_charge_gaus(list_in, file_item, num_bins= 20, range_default = None, xaxis = 'Time Resolution(ns)',
                        ylable = 'Occurrence', title = 'r$\delta$', pic = True, pdf = False):
    entries = len(list_in)
    mu, std = norm.fit(list_in)
    fig, ax4= plt.subplots(dpi=200)
    n,bins,patches=ax4.hist(list_in, bins=num_bins, range=range_default, density=False, 
                            label = 'entries: %d\nstd:%f fC\nmean:%f fC'%(entries, std, mu))
    xmin = np.min(list_in)
    xmax = np.max(list_in)
    x = np.linspace(xmin, xmax, 100)
    p = norm.pdf(x, mu, std)
    ax4.plot(x, np.max(n)*p/np.max(p), 'g', linewidth=2)
    ax4.grid()
    ax4.legend()
    ax4.set(xlabel=xaxis, ylabel='Occurrence',
               title='Charge Distribution,'+title)
    if pdf==True:
        pp.savefig(fig)
    if pic==True:
        plt.show()
    plt.close(fig)


def plot_distribution_time_Lindsey(list_in, file_item, num_bins= 20, range_default = None, xaxis = 'Time Resolution(ns)',
                          ylable = 'Occurrence', title = 'r$\delta$', pic = False, pdf = True):
    
    to_fit = list_in
    mean = np.mean(to_fit)
    sigma = np.std(to_fit, ddof=1)
    bins, edges = np.histogram(to_fit, num_bins, density=False)
    centers = 0.5*(edges[1:] + edges[:-1])
    popt, pcov = curve_fit(gaus,centers,bins,p0=[1,mean,sigma])
    mean = popt[1]
    sigma = popt[2]
    
    entries=len(to_fit)
    fig, ax4= plt.subplots(dpi=200)
    ax4.hist(to_fit, bins=num_bins, range=range_default, density=False, 
                            label = 'entries: %d\nstd: %.4f\nmean: %.2f'%(entries, abs(sigma), mean))
    ax4.plot(centers, gaus(centers,popt[0], popt[1], popt[2]), 'g', linewidth=2)
    ax4.grid()
    ax4.legend()
    ax4.set(xlabel=xaxis, ylabel='Occurrence',
               title=title)
    if pdf==True:
        pp.savefig(fig)
    if pic==True:
        plt.show()
    plt.close(fig)


def func1(x, a, b, c, d):
    return a + b*(x**1) + c*(x**2) + d*(x**3)

def func0(x, a, b):
    return a + b*(x**1)    

import matplotlib
def heatmap(data, row_labels, col_labels, ax=None,
            cbar_kw={}, cbarlabel="", **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (N, M).
    row_labels
        A list or array of length N with the labels for the rows.
    col_labels
        A list or array of length M with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    """

    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # We want to show all ticks...
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    # ... and label them with the respective list entries.
    ax.set_xticklabels(col_labels)
    ax.set_yticklabels(row_labels)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=True, bottom=False,
                   labeltop=True, labelbottom=False)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=-30, ha="right",
             rotation_mode="anchor")

    # Turn spines off and create white grid.
    for edge, spine in ax.spines.items():
        spine.set_visible(False)

    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar


def annotate_heatmap(im, data=None, valfmt="{x:.2f}",
                     textcolors=["black", "white"],
                     threshold=None, **textkw):
    """
    A function to annotate a heatmap.

    Parameters
    ----------
    im
        The AxesImage to be labeled.
    data
        Data used to annotate.  If None, the image's data is used.  Optional.
    valfmt
        The format of the annotations inside the heatmap.  This should either
        use the string format method, e.g. "$ {x:.2f}", or be a
        `matplotlib.ticker.Formatter`.  Optional.
    textcolors
        A list or array of two color specifications.  The first is used for
        values below a threshold, the second for those above.  Optional.
    threshold
        Value in data units according to which the colors from textcolors are
        applied.  If None (the default) uses the middle of the colormap as
        separation.  Optional.
    **kwargs
        All other arguments are forwarded to each call to `text` used to create
        the text labels.
    """

    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(data.max())/2.

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center",
              verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
            text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
            texts.append(text)

    return texts




def tw_correction (toa, tot, poly3rd=True, LindseyFit=True):
    if poly3rd == True:
        popt_toa, pcov_toa = curve_fit(func1, tot, toa)
    else:
        popt_toa, pcov_toa = curve_fit(func0, tot, toa)
    
    #################################
    #time walk correction and fit
    #################################



    # Generate data 
    if poly3rd == True:
        toa_fitted=func1(tot, *popt_toa)
    else:
        toa_fitted=func0(tot, *popt_toa)
    
    toa_corrected = toa - toa_fitted
    
    
    if LindseyFit == True:
        mu_corrected, sigma_corrected, popt_corrected, pcov_corrected = gaus_fit(to_fit=toa_corrected, num_bins=50)
    else:
        mu_corrected, sigma_corrected = norm.fit(toa_corrected)

    return mu_corrected, abs(sigma_corrected), toa_corrected