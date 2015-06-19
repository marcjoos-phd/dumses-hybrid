#===============================================================================
## \file debug.py
# \brief
# \b DUMSES-Hybrid:
# These are debug and comparison functions for DUMSES outputs and DUMSES history
# \author
# Marc Joos <marc.joos@cea.fr>, Sebastien Fromang
# \copyright
# Copyrights 2013-2015, CEA.
# This file is distributed under the CeCILL-A & GNU/GPL licenses, see
# <http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html> and
# <http://www.gnu.org/licenses/>
# \date
# \b created:          05-04-2015
# \b last \b modified: 05-12-2015

#===============================================================================
import sys
import numpy as np
import pylab as pl
import matplotlib as mp
from dumpy import data as dp

def comp(ndump, var="rho", dir1="old/", dir2="new/", xmin=None, xmax=None \
        , ymin=None, ymax=None, zmin=0, zmax=1, ratio=False, save=False \
        , *args, **kwargs):
    """Compare two outputs from two different DUMSES runs.

    Usage:
       comp(ndump, [var=var, dir1=dir1, dir2=dir2, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, zmin=zmin, zmax=zmax, ratio=ratio])

    With:
       ndump      index of the output to read
       var        name (or list of names) of variable to plot. Use the keyword 'all' to plot all the fields. Default: 'rho' 
       dir1       directory of the first DUMSES output to read. Default: 'old/'
       dir2       directory of the second DUMSES output to read. Default: 'new/'
       xmin, xmax boundaries of the domain to plot in the x-direction. Default: None
       ymin, ymax boundaries of the domain to plot in the y-direction. Default: None
       zmin, zmax boundaries of the domain to plot in the x-direction. Default: 0, 1
       ratio      set to 'True' to plot the ratio of each field. Default: False
       save       set to 'True' to save the plot in PNG, PDF and EPS formats
       
    Examples:
       comp(42)
       comp(42, 'all', dir1="ot/oacc/", dir2="ot/ompi/")
       comp(42, ['rho', 'Bx', 'E'], xmin=0, xmax=1, zmin=0, zmax=16)
       """
    do = dp.DumsesData()
    do.load(ndump, dir1)
    dn = dp.DumsesData()
    dn.load(ndump, dir2)
    dod = dp.dictionarize(do)
    dnd = dp.dictionarize(dn)

    if var == "all":
        var = dod.keys()
        if "Brx" in var:
            var.remove("Brx")
            var.remove("Bry")
            var.remove("Brz")
    if not isinstance(var, list): var = [var]

    if not xmin: xmin = 0
    if not xmax: xmax = dod[dod.keys()[0]].shape[0]
    if not ymin: ymin = 0
    if not ymax: ymax = dod[dod.keys()[0]].shape[1]
    if not zmin: zmin = 0
    if not zmax: zmax = dod[dod.keys()[0]].shape[2]
    ldim = np.array([xmax - xmin, ymax - ymin, zmax - zmin])
    ldim = ldim[ldim > 1]
    try:
        assert(ldim.size == 2)
    except AssertionError:
        print('Error: Wrong dimensions! Check your {x,y,z} boundaries.')
        sys.exit(42)
    xdim, ydim = ldim

    if zmax - zmin == 1: 
        xlab = "x"
        ylab = "y"
    if ymax - ymin == 1: 
        xlab = "x"
        ylab = "z"
    if xmax - xmin == 1: 
        xlab = "y"
        ylab = "z"

    if ratio:
        fig, ax = pl.subplots(3, len(var), figsize=(3.5*len(var),6))
    else:
        fig, ax = pl.subplots(2, len(var), figsize=(4.*len(var),6))
    pl.subplots_adjust(hspace=0.5)
    if len(ax.shape) == 1: ax = ax[:, np.newaxis]

    for i, cvar in enumerate(var):
        if cvar == "rho":  
            title = r"$\rho$"
        elif "rhou" in cvar: 
            title = r"$\rho v_" + cvar[-1] + "$"
        elif "B" in cvar:    
            title = "B$_" + cvar[-1] + "$"
        else:
            title = cvar

        dop = dod[cvar][xmin:xmax,ymin:ymax,zmin:zmax].reshape((xdim, ydim))
        ax[0,i].set_title(title)
        if i == 0: ax[0,i].set_ylabel(ylab)
        map1 = ax[0,i].imshow(dop, origin='lower', *args, **kwargs)
        ax[0,i].xaxis.set_major_locator(mp.ticker.MaxNLocator(4))
        ax[0,i].yaxis.set_major_locator(mp.ticker.MaxNLocator(4))
        fig.colorbar(map1, ax=ax[0,i])

        dnp = dnd[cvar][xmin:xmax,ymin:ymax,zmin:zmax].reshape((xdim, ydim))
        if i == 0: ax[1,i].set_ylabel(ylab)
        map2 = ax[1,i].imshow(dnp, origin='lower', *args, **kwargs)
        ax[1,i].xaxis.set_major_locator(mp.ticker.MaxNLocator(4))
        ax[1,i].yaxis.set_major_locator(mp.ticker.MaxNLocator(4))
        fig.colorbar(map2, ax=ax[1,i])

        if ratio:
            if i == 0: ax[2,i].set_ylabel(ylab)
            ax[2,i].set_xlabel(xlab)
            map3 = ax[2,i].imshow(dnp/dop, origin='lower', *args, **kwargs)
            ax[2,i].xaxis.set_major_locator(mp.ticker.MaxNLocator(4))
            ax[2,i].yaxis.set_major_locator(mp.ticker.MaxNLocator(4))
            fig.colorbar(map3, ax=ax[2,i])
        else:
            ax[1,i].set_xlabel(xlab)

    if save:
        for format_ in ['png', 'pdf', 'eps']:
            pl.savefig("comp_" + "".join(name[:3].title() for name in var) \
                           + "." + format_)
    else:
        pl.show()

def compHistory(var="maxwell", dir1="old/", dir2="new/", tmin=None, tmax=None \
                    , log=False, save=False, *args, **kwargs):
    """Compare the history.txt of two different DUMSES runs.

    Usage:
        compHistory([var=var, dir1=dir1, dir2=dir2, tmin=tmin, tmax=tmax])

    With:
        var        name (or list of names) of variable to plot. Use the keyword 'all' to plot all the fields. Default: 'rho' 
       dir1        directory of the first DUMSES output to read. Default: 'old/'
       dir2        directory of the second DUMSES output to read. Default: 'new/'
       tmin, tmax  boundaries of the time domain to plot. Default: 'None'
       log         plot in logscale. Default: 'False'
       save        set to 'True' to save the plot in PNG, PDF and EPS formats
    """
    mp.rcParams['legend.fontsize'] = 'small'

    do = dp.DumsesHistory()
    do.load(dir1)
    dn = dp.DumsesHistory()
    dn.load(dir2)

    label1 = dir1
    label2 = "reference" if "reference" in dir2 else dir2

    if var == "all":
        var = dn.dict.keys()
    if not isinstance(var, list): var = [var]

    fig, ax = pl.subplots(len(var),1, figsize=(3.5,2.5*len(var)))
    if not isinstance(ax, np.ndarray): ax = np.array([ax])
    for i, cvar in enumerate(var):
        if "time" in dn.dict:
            to0 = do.dict['time'] > tmin
            to1 = do.dict['time'] < tmax
            tn0 = dn.dict['time'] > tmin
            tn1 = dn.dict['time'] < tmax
            if not to1.any(): to1 = True
            if not tn1.any(): tn1 = True
            ymin = min(do.dict[cvar][to0 & to1].min() \
                     , dn.dict[cvar][tn0 & tn1].min())
            ymax = max(do.dict[cvar][to0 & to1].max() \
                     , dn.dict[cvar][tn0 & tn1].max())
            vmin = ymin - abs(ymax - ymin)/10.
            vmax = ymax + abs(ymax - ymin)/10.
            if log:
                ax[i].semilogy(do.dict["time"], do.dict[cvar], label='from: ' \
                                   + label1)
                ax[i].semilogy(dn.dict["time"], dn.dict[cvar], label='from: ' \
                                   + label2)
            else:
                ax[i].plot(do.dict["time"], do.dict[cvar], label='from: ' \
                               + label1)
                ax[i].plot(dn.dict["time"], dn.dict[cvar], label='from: ' \
                               + label2)
            ax[i].legend(*args, **kwargs)
            ax[i].set_xlim(xmin=tmin, xmax=tmax)
            ax[i].set_ylim(ymin=vmin, ymax=vmax)
            if (i == len(var) - 1): ax[i].set_xlabel("time")
        else:
            if log:
                ax[i].semilogy(do.dict[cvar], label='from: ' + label1)
                ax[i].semilogy(dn.dict[cvar], label='from: ' + label2)
            else:
                ax[i].plot(do.dict[cvar], label='from: ' + label1)
                ax[i].plot(dn.dict[cvar], label='from: ' + label2)
            ax[i].legend(*args, **kwargs)
            if (i == len(var) - 1): ax[i].set_xlabel(r"n$_{\sf hist}$")
        ax[i].set_ylabel(cvar)
        ax[i].xaxis.set_major_locator(mp.ticker.MaxNLocator(4))
    fig.tight_layout()
    if save:
        for format_ in ['png', 'pdf', 'eps']:
            pl.savefig("compHist_" + "".join(name[:3].title() for name in var) \
                           + "." + format_)
    else:
        pl.show()
