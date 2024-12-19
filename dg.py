####Analysis of free energy calculations performed with GROMACS####
####script takes GROMACS .xvg files and computes and plots the cumulative dG estimate at every ns 
####Plot can then be used to check convergence
####dG is estimated using MBAR (alchemlyb https://github.com/alchemistry/alchemlyb/tree/master/src/alchemlyb and pymbar https://github.com/choderalab/pymbar/tree/master/pymbar)
#
# Use the following versions in your conda environment
# pymbar                    4.0.3                h7900ff3_1    conda-forge
# alchemlyb                 2.3.1              pyhd8ed1ab_0    conda-forge

import pickle
from alchemlyb.parsing.gmx import extract_dHdl, extract_u_nk
import os.path
from os.path import join
import pandas as pd
from alchemlyb.parsing import gmx
from alchemlyb.preprocessing import (slicing, statistical_inefficiency,
                                     equilibrium_detection)
from alchemlyb.estimators import TI
from alchemlyb.estimators import MBAR, BAR
import argparse
import pymbar
from pymbar.timeseries import subsample_correlated_data
from pymbar.timeseries import statistical_inefficiency as PTstatistical_inefficiency
import numpy as np
import matplotlib.pyplot as plt
#from alchemlyb.visualisation import plot_mbar_overlap_matrix

def plot_mbar_overlap_matrix(matrix, skip_lambda_index=[], ax=None):
    '''Plot the MBAR overlap matrix.
    Parameters
    ----------
    matrix : numpy.matrix
        DataFrame of the overlap matrix obtained from
        :attr:`~alchemlyb.estimators.MBAR.overlap_matrix`
    skip_lambda_index : List
        list of lambda indices to be omitted from plotting process.
        Default: ``[]``.
    ax : matplotlib.axes.Axes
        Matplotlib axes object where the plot will be drawn on. If ``ax=None``,
        a new axes will be generated.
    Returns
    -------
    matplotlib.axes.Axes
        An axes with the overlap matrix drawn.
    Note
    ----
    The code is taken and modified from
    `Alchemical Analysis <https://github.com/MobleyLab/alchemical-analysis>`_.
    .. versionadded:: 0.4.0
    '''
    # Compute the size of the figure, if ax is not given.
    max_prob = matrix.max()
    size = len(matrix)
    if ax is None:
        fig, ax = plt.subplots(figsize=(size / 2, size / 2))
    ax.set_xticks([])
    ax.set_yticks([])
    ax.axis('off')
    for i in range(size):
        if i != 0:
            ax.axvline(x=i, ls='-', lw=0.5, color='k', alpha=0.25)
            ax.axhline(y=i, ls='-', lw=0.5, color='k', alpha=0.25)
        for j in range(size):
            if matrix[j, i] < 0.005:
                ii = ''
            elif matrix[j, i] > 0.995:
                ii = '1.00'
            else:
                ii = ("{:.2f}".format(matrix[j, i])[1:])
            alf = matrix[j, i] / max_prob
            ax.fill_between([i, i + 1], [size - j, size - j],
                            [size - (j + 1), size - (j + 1)], color='k',
                            alpha=alf)
            ax.annotate(ii, xy=(i, j), xytext=(i + 0.5, size - (j + 0.5)),
                        size=8, textcoords='data', va='center',
                        ha='center', color=('k' if alf < 0.5 else 'w'))

    if skip_lambda_index:
        ks = [int(l) for l in skip_lambda_index]
        ks = np.delete(np.arange(size + len(ks)), ks)
    else:
        ks = range(size)
    for i in range(size):
        ax.annotate(ks[i], xy=(i + 0.5, 1), xytext=(i + 0.5, size + 0.5),
                    size=10, textcoords=('data', 'data'),
                    va='center', ha='center', color='k')
        ax.annotate(ks[i], xy=(-0.5, size - (size - 0.5)),
                    xytext=(-0.5, size - (i + 0.5)),
                    size=10, textcoords=('data', 'data'),
                    va='center', ha='center', color='k')
    ax.annotate(r'$\lambda$', xy=(-0.5, size - (size - 0.5)),
                xytext=(-0.5, size + 0.5),
                size=10, textcoords=('data', 'data'),
                va='center', ha='center', color='k')
    ax.plot([0, size], [0, 0], 'k-', lw=4.0, solid_capstyle='butt')
    ax.plot([size, size], [0, size], 'k-', lw=4.0, solid_capstyle='butt')
    ax.plot([0, 0], [0, size], 'k-', lw=2.0, solid_capstyle='butt')
    ax.plot([0, size], [size, size], 'k-', lw=2.0, solid_capstyle='butt')

    cx = np.repeat(range(size + 1), 2)
    cy = sorted(np.repeat(range(size + 1), 2), reverse=True)
    ax.plot(cx[2:-1], cy[1:-2], 'k-', lw=2.0)
    ax.plot(np.array(cx[2:-3]) + 1, cy[1:-4], 'k-', lw=2.0)
    ax.plot(cx[1:-2], np.array(cy[:-3]) - 1, 'k-', lw=2.0)
    ax.plot(cx[1:-4], np.array(cy[:-5]) - 2, 'k-', lw=2.0)

    ax.set_xlim(-1, size)
    ax.set_ylim(0, size + 1)
    return ax

def test(path, lambdanr):
    counts = []
    for i in range(0, lambdanr):
        count = 0
        start = False
        xvg = 'lambda_%i/Production_MD/prod.xvg' %(i)
        # xvg = '%i/lambda_%i/Production_MD/prod.xvg' %(repeat,i)
        u_nk = open(xvg,'r')
        for x in u_nk:
            if x.startswith('0.0000'):
                start = True
            else:
                # divide by 5 in order to only get it every nanosecond
                if start==True:
                    count += 1
                else:
                    continue
        counts.append(count)
    return min(counts)

#specify how many ps to discard for equilibration
discard = 1000

start = discard+1000
#evaluate the dG every 1000ps for cumulative dG analysis (assess convergence)
protocols = ['step1']
#repeats = [1]
#for r in repeats:
for p in protocols:
    if p == 'step1':
        lambdanr = 20
    try:
        path = ""
        upper = test(path,lambdanr)
        #upper = test(path,r,lambdanr)
        upper = int(upper)
        end = upper+1
        print("start", start, "end", end)
        timeslices = list(range(start,end,1000))
        print("timeslices:", timeslices)
        u_nks_all = []
        #print(r)
        for i in range(0, lambdanr):
            folder = 'lambda_%i/Production_MD'%(i)
            xvg = '%s/prod.xvg' %(folder)
            u_nk = extract_u_nk(xvg, T=303.15)
            print(len(u_nk))
            u_nks_all.append(u_nk)
        dGs = []
        dG_error = []
        for time in timeslices:
            print("TIME:", time)
            u_nks = []
            for i in range(0, lambdanr):
                #Slice data up to that time frame
                u_nk = slicing(u_nks_all[i], upper=time, lower=discard)
                #choose the potential energy of the lambda state of the trajectory for autocorreltation analysis
                #could be changed to always using the potential energy of the same lambda state
                series = u_nk[u_nk.columns[i]]
                #run autocorrelation analysis
                statinef  = PTstatistical_inefficiency(u_nk[u_nk.columns[i]], fast=False)
                indices = subsample_correlated_data(series, g=statinef,
                                                  conservative=True)
                #If there are less than 100 uncorrelated samples: use 100 samples as default
                if len(indices) >= 100:
                    u_nk = u_nk.iloc[indices]
                else:
                    #A low number of uncorrelated samples can indicate a slow motion
                    print('time %i: Low number of uncorrelated samples (%i) in lambda %i. Check for slow motion.'%(time, len(indices), i)) 
                    #number of data points in this lambda window
                    number = int(len(u_nk.index))
                    #extract 100 data points
                    skip = int(number/100)
                    #Slice data to obtain 100 samples
                    u_nk = slicing(u_nk, step=skip)
                u_nks.append(u_nk)
            #Estimate dG using alchemlyb/pymbar
            u_nk = pd.concat(u_nks)
            index_2 = u_nk.index.values
            mbar = MBAR()
            print(u_nk.head())
            fit = mbar.fit(u_nk)
            groups = u_nk.groupby(level=u_nk.index.names[1:])
            N_k = [(len(groups.get_group(i)) if i in groups.groups else 0) for i in u_nk.columns]       
            lam =  mbar.delta_f_.columns.values.tolist()
            endpoint_diff = 0.596 * mbar.delta_f_.loc[[lam[0]], [lam[-1]]]
            dg_diff = round(float(endpoint_diff.iloc[0][lam[-1]]),3)
            dGs.append(dg_diff)
            uncertainty = 0.596 * mbar.d_delta_f_.loc[[lam[0]], [lam[-1]]]
            uncertainty = round(float(uncertainty.iloc[0][lam[-1]]),3)
            dG_error.append(uncertainty)
            if time == timeslices[-1]:
                ####Overlap Matrix
                u_kn = u_nk.transpose()
                mbar_pymbar = pymbar.mbar.MBAR(u_kn, N_k)
                overlap_matrix = mbar_pymbar.compute_overlap()
                overlap_matrix = overlap_matrix['matrix']
                ax = plot_mbar_overlap_matrix(overlap_matrix)
                ax.figure.savefig('O_MBAR.pdf', bbox_inches='tight', pad_inches=0.0)    
        print(dGs)
        print(dG_error)
        file = open('dg.pickle', 'wb')
        pickle.dump(dGs, file)
        file.close()
        file = open('dg_error.pickle', 'wb')
        pickle.dump(dG_error, file)
        file.close()
    except (ValueError, ) as e:
        print(f"!!!Error: {e}")
        continue

#Plot the cumulative dG estimate as a function of simulation time in ns
x_e = np.arange(start/1000,len(dGs)+start/1000,1)
plt.errorbar(x_e, dGs, yerr=dG_error, alpha=0.5, fmt = '-o', mfc='w', mew=2.5, lw=5, ms=10)
plt.xlabel('simulation time in ns',fontsize=20)
plt.ylabel('ΔG in kcal/mol',fontsize=20)
plt.grid(True)    
plt.savefig('dG.png', bbox_inches = "tight")
plt.show()

