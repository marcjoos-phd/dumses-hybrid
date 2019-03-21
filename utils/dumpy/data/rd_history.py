#===============================================================================
## \file rd_history.py
# \brief
# \b DUMSES-Hybrid:
# These are data class for DUMSES history file
# \author
# Marc Joos <marc.joos@cea.fr>, Sebastien Fromang
# \copyright
# Copyrights 2013-2015, CEA.
# This file is distributed under the CeCILL-A & GNU/GPL licenses, see
# <http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html> and
# <http://www.gnu.org/licenses/>
# \date
# \b created:          06-17-2014 
# \b last \b modified: 05-13-2015

#===============================================================================
import sys
import numpy as np
import pylab as pl

bold    = "\033[1m"
reset   = "\033[0;0m"

## DUMSES history class
class DumsesHistory:
    """DumsesHistory class:
    Use this class to read your DUMSES history.

    Usage:
      History = DumsesHistory()
      History.load([filedir, listKeys])

    With:
      filedir:  path of the directory containing the history.txt file
      listKeys: use if you don't use the default variables of DUMSES (time', 'dt', 'mass', 'maxwell', 'reynolds', 'max+rey', 'magp', 'meanBx', 'meanBy', 'meanBz', 'divB' for 3D simulations), to specify the variables to read.

    It is then possible to plot variables with the plotting methods:
      History.plot(variable [, tmin, tmax])
      History.plotVs(variable1, variable2 [, tmin, tmax])

    With:
      variable[1,2]: the name of the variable(s) to plot
      tmin, tmax:    lower and upper boundary in time (only if 'time' is defined)
    """
    ## Initialization method
    def __init__(self):
        ## directory of the \c 'history.txt' file
        self.dir = ""
        ## dictionary to store data from \c 'history.txt'
        self.dict = {}

    ## load method to read \c 'history.txt'
    # @param filedir directory of the history file
    # @param listKeys list of variable names (for backward compatibility only)
    def load(self, filedir="./", listKeys=None):
        """DumsesHistory load method.
        
        see DumsesHistory docstring for more informations."""
        fname = filedir + "history.txt"
        f = open(fname)
        keys = f.readline()
        if keys[0] != "#":
            print("Old history format detected")
            if not listKeys:
                try:
                    assert(len(keys.split()) in [4, 11])
                except AssertionError:
                    print("Error: unknown number of variables!")
                    print("Please specify the name of your variables with the 'listKeys' argument.")
                    sys.exit(2)
                if len(keys.split()) == 4:
                    keys = ['time', 'rho', 'rhoux', 'p']
                elif len(keys.split()) == 11:
                    keys = ['time', 'dt', 'mass', 'maxwell', 'reynolds' \
                          , 'max+rey', 'magp', 'meanBx', 'meanBy', 'meanBz' \
                          , 'divB']
                print("Assume variable names to be:")
                print(", ".join(keys))
            else:
                try:
                    assert(isinstance(listKeys, list))
                except AssertionError:
                    print("'listKeys' as to be a list!")
                    sys.exit(2)
                keys = listKeys
        else:
            keys = keys.split()[1:]
        f.close()
        for i in range(len(keys)):
            self.dict[keys[i]] = np.loadtxt(fname, usecols=(i,))

    ## listvar method to list the variables written in \c 'history.txt'
    def listvar(self):
        """DumsesHistory listvar method.

        Returns the list of variables read from 'history.txt'."""
        sys.stdout.write(bold + "'history.txt' contains the following variables: \n" + reset)
        for i, ivar in enumerate(self.dict.iterkeys()):
            sys.stdout.write(ivar)
            if i < self.dict.__len__() - 1: sys.stdout.write(', ')
            
    ## plot method
    # @param var variable to plot; check variable names with DumsesHistory.listvar()
    # @param tmin set lower boundary in time range
    # @param tmax set upper boundary in time range
    def plot(self, var="maxwell", tmin=None, tmax=None, save=False):
        """DumsesHistory plot method.

        Usage:
          History.plot([var, tmin, tmax])
        With:
          var:  name of the variable to plot. Default: 'maxwell'
          tmin: if time is defined, lower time boundary. Default: 'None'
          tmax: if time is defined, upper time boundary. Default: 'None'
          save: if 'True', save the plot in PNG, PDF and EPS formats
        """
        if var in self.dict:
            pl.figure(figsize=(6,4))
            pl.subplots_adjust(bottom=0.15)
            if "time" in self.dict:
                t0 = self.dict['time'] > tmin
                t1 = self.dict['time'] < tmax
                if not t1.any(): t1 = True
                vmin = self.dict[var][t0 & t1].min()
                vmax = self.dict[var][t0 & t1].max()
                pl.plot(self.dict["time"], self.dict[var])
                pl.xlim(xmin=tmin, xmax=tmax)
                pl.ylim(ymin=vmin, ymax=vmax)
                pl.xlabel("time")
            else:
                pl.plot(self.dict[var])
                pl.xlabel(r"n$_{\sf hist}$")
            pl.ylabel(var)
            if save:
                for format_ in ['png', 'pdf', 'eps']:
                    pl.savefig(var + "." + format_)
            else:
                pl.show()
        else:
            print(bold + "Error: " + reset \
                       + "the variable %s does not exist" %var)

    ## plotVs method: plot two variables simultaneously
    # @param var1 first variable to plot
    # @param var2 second variable to plot
    # @param tmin set lower boundary in time range
    # @param tmax set upper boundary in time range
    def plotVs(self, var1="maxwell", var2="reynolds", tmin=None, tmax=None \
                   , save=False):
        """DumsesHistory plotVs method: plot two variables from 'history.txt'

        Usage:
          History.plot([var1, var2, tmin, tmax])
        With:
          var1: name of the variable to plot. Default: 'maxwell'
          var2: name of the variable to plot. Default: 'reynolds'
          tmin: if time is defined, lower time boundary. Default: 'None'
          tmax: if time is defined, upper time boundary. Default: 'None'
          save: if 'True', save the plot in PNG, PDF and EPS formats
        """
        if (var1 in self.dict) & (var2 in self.dict):
            fig, ax1 = pl.subplots(figsize=(6,4))
            pl.subplots_adjust(bottom=0.15)
            if "time" in self.dict:
                t0 = self.dict['time'] > tmin
                t1 = self.dict['time'] < tmax
                if not t1.any(): t1 = True
                vmin = self.dict[var1][t0 & t1].min()
                vmax = self.dict[var1][t0 & t1].max()
                ax1.plot(self.dict["time"], self.dict[var1], "b")
                ax1.set_xlim(xmin=tmin, xmax=tmax)
                ax1.set_ylim(ymin=vmin, ymax=vmax)
                ax1.set_xlabel("time")
            else:
                ax1.plot(self.dict[var1], "b")
                ax1.set_xlabel(r"n$_{\sf hist}$")
            ax1.set_ylabel(var1)
            ax2 = ax1.twinx()
            if "time" in self.dict:
                t0 = self.dict['time'] > tmin
                t1 = self.dict['time'] < tmax
                if not t1.any(): t1 = True
                vmin = self.dict[var2][t0 & t1].min()
                vmax = self.dict[var2][t0 & t1].max()
                ax2.plot(self.dict["time"], self.dict[var2], "g")
                ax2.set_xlim(xmin=tmin, xmax=tmax)
                ax2.set_ylim(ymin=vmin, ymax=vmax)
                ax2.set_xlabel("time")
            else:
                ax2.plot(self.dict[var2], "g")
                ax2.set_xlabel(r"n$_{\sf hist}$")
            ax2.set_ylabel(var2)
            if save:
                for format_ in ['png', 'pdf', 'eps']:
                    pl.savefig(var1 + "_" + var2 + "." + format_)
            else:
                pl.show()
        else:
            print(bold + "Error: " + reset \
                       + "the variable %s does not exist" %var)

