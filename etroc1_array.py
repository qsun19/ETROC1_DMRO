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