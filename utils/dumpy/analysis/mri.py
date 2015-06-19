#===============================================================================
## \file mri.py
# \brief
# \b DUMSES-Hybrid:
# These are functions to compute alpha parameter in MRI simulations
# \author
# Marc Joos <marc.joos@cea.fr>, Sebastien Fromang
# \copyright
# Copyrights 2013-2015, CEA.
# This file is distributed under the CeCILL-A & GNU/GPL licenses, see
# <http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html> and
# <http://www.gnu.org/licenses/>
# \date
# \b created:          05-13-2015 
# \b last \b modified: 05-13-2015

#===============================================================================
import sys
import numpy as np
import pylab as pl
from dumpy import data as dp

bold    = "\033[1m"
reset   = "\033[0;0m"

## Formatter for alpha print in compAlpha
def formatter(x, nbd=2):
    sciexp   = '%e' %x
    realpart = sciexp.split('e')[0]
    exppart  = sciexp.split('e')[1]
    if len(exppart.split('+')) == 1:
        sign    = '-'
        exppart = exppart.split('-')[1]
    else:
        sign    = '+'
        exppart = exppart.split('-')[1]

    if nbd == 1: realpart = '%2.1f' %eval(realpart)
    if nbd == 2: realpart = '%3.2f' %eval(realpart)
    if nbd >= 3: realpart = '%4.3f' %eval(realpart)

    expr = realpart + r'$ \times $10$^{\sf %s%s}$' %(sign, '%d' %eval(exppart))
    return(expr)

## compAlpha function, to compute Alpha out of \c 'history.txt'
# @param data a DumsesHistory object
# @param nb_orbit number of orbits on which to perform the computation
# @param t_orb elapsed time during one orbit (in code unit)
# @param ciso sound speed
# @param plot set to \c True to plot the results
# @param save set to \c True to save the plot in PNG, PDF and EPS formats
## If \c time is not defined, plots are done with respect \c n_hist.
## if \c maxwell, \c reynolds and \c max+rey variables are not defined,
## the computation will fail.
def compAlpha(data, nb_orb=1., t_orb=6283., ciso=1e-3, plot=False \
                  , save=False):
    """compAlpha method, to compute Alpha out of \c 'history.txt'.

    Usage:
      compAlpha(data, [nb_orb, t_orb, ciso, plot])
    With:
      data:     a DumsesHistory object
      nb_orbit: number of orbits on which to perform the computation
      t_orb:    elapsed time during one orbit (in code unit)
      ciso:     sound speed
      plot:     set to True to plot the results
      save:     set to True to save the plot in PNG, PDF and EPS formats
    """
    try:
        assert(isinstance(data, dp.DumsesHistory))
    except AssertionError:
        print("Error: 'data' has to be a DumsesHistory object!")
        sys.exit(42)

    try:
        assert("maxwell" in data.dict)
    except AssertionError:
        print(bold + "Error: " + reset + "'maxwell' is not in 'history.txt'")
    try:
        assert("reynolds" in data.dict)
    except AssertionError:
        print(bold + "Error: " + reset + "'reynolds' is not in 'history.txt'")
    try:
        assert("max+rey" in data.dict)
    except AssertionError:
        print(bold + "Error: " + reset + "'max+rey' is not in 'history.txt'")
    maxwell = data.dict["maxwell"]
    reynolds = data.dict["reynolds"]
    maxnolds = data.dict["max+rey"]

    if "time" in data.dict:
        time = data.dict["time"]
        tmax = time[-1]
        tlim = nb_orb*t_orb
        tn   = pl.find(time >= (tmax-tlim))
        ntn  = tn.size
        xlabel = "time"
    else:
        tmax = data.dict.values()[0].size-1
        time = np.linspace(0, tmax, tmax+1)
        tlim = tmax - nb_orb
        tn   = map(int, np.linspace(tlim, tmax, nb_orb))
        xlabel = r"n$_{\sf hist}$"

    alpha     = np.sum(maxnolds[tn])/ntn/ciso**2
    alpha_max = np.sum(maxwell[tn])/ntn/ciso**2
    alpha_rey = np.sum(reynolds[tn])/ntn/ciso**2
    sigma_alpha = np.sum((maxnolds[tn]/ciso**2-alpha)**2)/ntn
    sigma_alpha_max = np.sum((maxwell[tn]/ciso**2-alpha_max)**2)/ntn
    sigma_alpha_rey = np.sum((reynolds[tn]/ciso**2-alpha_rey)**2)/ntn
    sigma_alpha = np.sqrt(sigma_alpha)
    sigma_alpha_max = np.sqrt(sigma_alpha_max)
    sigma_alpha_rey = np.sqrt(sigma_alpha_rey)
    print("For %d orbits, " %nb_orb + "Alpha = %5.4f" %alpha \
              + ' +/- %3.2e' %sigma_alpha)
    print("           " + "Alpha_Max = %5.4f" %alpha_max \
              + ' +/- %3.2e' %sigma_alpha_max)
    print("           " + "Alpha_Rey = %5.4f" %alpha_rey \
              + ' +/- %3.2e' %sigma_alpha_rey)
    
    if plot:
        fig = pl.figure(figsize=(8,12))
        # Alpha total
        ax = pl.subplot(311)
        pl.subplots_adjust(bottom=0.15)
        pl.ylabel(r'$\alpha$')
        pl.title(r'$\alpha_{\sf mean}$ = %5.4f' %alpha \
              + r' $\pm$ %s' %formatter(sigma_alpha))
        pl.plot(time/t_orb, maxnolds/ciso**2, color='w')
        x0,x1,y0,y1 = pl.axis()
        vert = [((tmax-tlim)/t_orb,y0)] + [((tmax-tlim)/t_orb,y1)] \
             + [(x1,y1)] + [(x1,y0)]
        poly = pl.Polygon(vert, fill=False, edgecolor='grey', hatch='//')
        ax.add_patch(poly)
        pl.plot(time/t_orb, maxnolds/ciso**2)
        pl.plot(time/t_orb \
              , pl.linspace(alpha,alpha,time.size) \
              , color='k', lw=1.)
        vert = [(min(time), alpha-sigma_alpha)] \
             + [(min(time), alpha+sigma_alpha)] \
             + [(x1, alpha+sigma_alpha)] \
             + [(x1, alpha-sigma_alpha)]
        poly = pl.Polygon(vert, facecolor='lightgrey' \
                              , edgecolor='lightgrey')
        ax.add_patch(poly)
        pl.axis([x0, x1, y0, y1])

        # Alpha Maxwell
        ax = pl.subplot(312)
        pl.subplots_adjust(bottom=0.15)
        pl.ylabel(r'$\alpha_{\sf Maxwell}$')
        pl.title(r'$\alpha_{\sf Maxwell\ mean}$ = %5.4f' %alpha_max \
              + r' $\pm$ %s' %formatter(sigma_alpha_max))
        pl.plot(time/t_orb, maxwell/ciso**2, color='w')
        x0,x1,y0,y1 = pl.axis()
        vert = [((tmax-tlim)/t_orb,y0)] + [((tmax-tlim)/t_orb,y1)] \
             + [(x1,y1)] + [(x1,y0)]
        poly = pl.Polygon(vert, fill=False, edgecolor='grey', hatch='//')
        ax.add_patch(poly)
        pl.plot(time/t_orb, maxwell/ciso**2)
        pl.plot(time/t_orb \
              , pl.linspace(alpha_max,alpha_max,time.size) \
              , color='k', lw=1.)
        vert = [(min(time), alpha_max-sigma_alpha_max)] \
             + [(min(time), alpha_max+sigma_alpha_max)] \
             + [(x1, alpha_max+sigma_alpha_max)] \
             + [(x1, alpha_max-sigma_alpha_max)]
        poly = pl.Polygon(vert, facecolor='lightgrey' \
                              , edgecolor='lightgrey')
        ax.add_patch(poly)
        pl.axis([x0, x1, y0, y1])

        # Alpha Reynolds
        ax = pl.subplot(313)
        pl.subplots_adjust(bottom=0.15)
        pl.xlabel(xlabel)
        pl.ylabel(r'$\alpha_{\sf Reynolds}$')
        pl.title(r'$\alpha_{\sf Reynolds\ mean}$ = %5.4f' %alpha_rey \
              + r' $\pm$ %s' %formatter(sigma_alpha_rey))
        pl.plot(time/t_orb, reynolds/ciso**2, color='w')
        x0,x1,y0,y1 = pl.axis()
        vert = [((tmax-tlim)/t_orb,y0)] + [((tmax-tlim)/t_orb,y1)] \
             + [(x1,y1)] + [(x1,y0)]
        poly = pl.Polygon(vert, fill=False, edgecolor='grey', hatch='//')
        ax.add_patch(poly)
        pl.plot(time/t_orb, reynolds/ciso**2)
        pl.plot(time/t_orb \
              , pl.linspace(alpha_rey,alpha_rey,time.size) \
              , color='k', lw=1.)
        vert = [(min(time), alpha_rey-sigma_alpha_rey)] \
             + [(min(time), alpha_rey+sigma_alpha_rey)] \
             + [(x1, alpha_rey+sigma_alpha_rey)] \
             + [(x1, alpha_rey-sigma_alpha_rey)]
        poly = pl.Polygon(vert, facecolor='lightgrey' \
                              , edgecolor='lightgrey')
        ax.add_patch(poly)
        pl.axis([x0, x1, y0, y1])

        if save:
            for form in ['png', 'pdf', 'eps']:
                namefig = 'alpha.' + form
                pl.savefig(namefig)
            pl.close()
        else:
            pl.show()
