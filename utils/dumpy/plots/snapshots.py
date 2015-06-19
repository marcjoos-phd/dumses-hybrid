#===============================================================================
## \file snapshots.py
# \brief
# \b DUMSES-Hybrid:
# These are functions to plot snapshots of DUMSES-Hybrid simulations
# \author
# Marc Joos <marc.joos@cea.fr>, Sebastien Fromang
# \copyright
# Copyrights 2013-2015, CEA.
# This file is distributed under the CeCILL-A & GNU/GPL licenses, see
# <http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html> and
# <http://www.gnu.org/licenses/>
# \date
# \b created:          05-13-2015 
# \b last \b modified: 06-16-2015

#===============================================================================
import re
import sys
import parser
import numpy as np
import pylab as pl
import matplotlib as mp
from dumpy import data as dp
try:
    import sympy
except ImportError:
    print("Warning: package 'sympy' not found, won't be able to use pretty formatting")

## plot1d function to plot a slice in a DUMSES-Hybrid simulation
# @param data a DumsesData object
# @param direction direction in which to plot, in {'x', 'y', 'z'}
# @param xpos position of the slice in the x-direction.
# @param ypos position of the slice in the y-direction.
# @param zpos position of the slice in the z-direction.
# @param xmin lower x-boundary of the plot
# @param xmax upper x-boundary of the plot
# @param ymin lower y-boundary of the plot
# @param ymax upper y-boundary of the plot
# @param zmin lower z-boundary of the plot
# @param zmax upper z-boundary of the plot
# @param physUnit if \c True, use physical units for positions and boundaries
# @param physTick if \c True, use physical units for plotting
# @param log if \c True, plot in logscale
# @param save if \c True, save the plot in PNG, PDF and EPS
def plot1d(data, var='rho', direction='x', xpos=None, ypos=None, zpos=None \
         , xmin=None, xmax=None, ymin=None, ymax=None, zmin=None, zmax=None \
         , physUnit=True, physTick=True, log=False, save=False \
         , *args, **kwargs):
    """plot1d function to plot variables in a 1D DUMSES-Hybrid simulation

    Usage:
      plot1d(data [, var, direction, xpos, ypos, zpos, xmin, xmax, xmin, xmax, zmin, zmax, physUnit, physTick, log, save, *args, **kwargs])
    With:
      data       a DumsesData object
      var        a variable or list of variables to plot
      direction  direction in which to plot, in {'x', 'y', 'z'}
      xpos       position of the slice in the x-direction. Default: None
      ypos       position of the slice in the y-direction. Default: None
      zpos       position of the slice in the z-direction. Default: None
      xmin lower x-boundary of the plot
      xmax upper x-boundary of the plot
      ymin lower y-boundary of the plot
      ymax upper y-boundary of the plot
      zmin lower z-boundary of the plot
      zmax upper z-boundary of the plot
      physUnit   if True, use physical units for positions and boundaries
      physTick   if True, use physical units for plotting
      log        if True, plot in logscale
      save       if True, save the plot in PNG, PDF and EPS
    """
    # Ensure that data is a DumsesData object
    try:
        assert(isinstance(data, dp.DumsesData))
    except AssertionError:
        print("Error: 'data' must be a DumsesData object!")
        sys.exit(42)

    dictDat = dp.dictionarize(data)
    
    if var == "all":
        var = dictDat.keys()
        if "Brx" in var:
            var.remove("Brx")
            var.remove("Bry")
            var.remove("Brz")
    if not isinstance(var, list): var = [var]
    
    try:
        assert(direction in ['x', 'y', 'z'])
    except AssertionError:
        print("Error: direction must be in {'x', 'y', 'z'}")
        sys.exit(42)

    # Determine the position, the boundaries and the labels of the plot
    if direction == "x":
        xdir = True; ydir = False; zdir = False
    elif direction == "y":
        xdir = False; ydir = True; zdir = False
    else:
        xdir = False; ydir = False; zdir = True

    if xdir:
        try:
            assert(data.nx > 1)
        except AssertionError:
            print('Error: you should specify another direction for your plot (xdim = 1)')
            sys.exit(42)
        try:
            assert(xpos == None)
        except AssertionError:
            print("Error: you cannot specify 'xpos' in a plot in the x-direction")
            sys.exit(42)
        if physUnit:
            if xmin: xmin = findPosition(xmin, data.x)
            if xmax: xmax = findPosition(xmax, data.x)
        if not xmin: xmin = 0
        if not xmax: xmax = dictDat[dictDat.keys()[0]].shape[0]
        if physTick:
            xlab  = "x/H"
            xlim  = (data.x[xmin], data.x[xmax-1])
            xplot = data.x[xmin:xmax]
        else:
            xlab = "x"
            xlim = (xmin, xmax-1)
    else:
        xmin = 0
        xmax = 1
    if ydir:
        try:
            assert(data.ny > 1)
        except AssertionError:
            print('Error: you should specify another direction for your plot (ydim = 1)')
            sys.exit(42)
        try:
            assert(ypos == None)
        except AssertionError:
            print("Error: you cannot specify 'ypos' in a plot in the y-direction")
            sys.exit(42)
        if physUnit:
            if ymin: ymin = findPosition(ymin, data.y)
            if ymax: ymax = findPosition(ymax, data.y)
        if not ymin: ymin = 0
        if not ymax: ymax = dictDat[dictDat.keys()[0]].shape[1]
        if physTick:
            xlab  = "y/H"
            xlim  = (data.y[ymin], data.y[ymax-1])
            xplot = data.y[ymin:ymax]
        else:
            xlab = "y"
            xlim = (ymin, ymax-1)
    else:
        ymin = 0
        ymax = 1
    if zdir:
        try:
            assert(data.nz > 1)
        except AssertionError:
            print('Error: you should specify another direction for your plot (zdim = 1)')
            sys.exit(42)
        try:
            assert(zpos == None)
        except AssertionError:
            print("Error: you cannot specify 'zpos' in a plot in the z-direction")
            sys.exit(42)
        if physUnit:
            if zmin: zmin = findPosition(zmin, data.z)
            if zmax: zmax = findPosition(zmax, data.z)
        if not zmin: zmin = 0
        if not zmax: zmax = dictDat[dictDat.keys()[0]].shape[2]
        if physTick:
            xlab  = "z/H"
            xlim  = (data.z[zmin], data.z[zmin-1])
            xplot = data.z[zmin:zmax]
        else:
            xlab = "z"
            xlim = (zmin, zmax-1)
    else:
        zmin = 0
        zmax = 1

    ldim = np.array([xmax - xmin, ymax - ymin, zmax - zmin])
    ldim = ldim[ldim > 1]
    try:
        assert(ldim.size == 1)
    except AssertionError:
        print('Error: Wrong dimensions! Check your {x,y,z} boundaries.')
        sys.exit(42)

    if xpos:
        if physUnit:
            xmin = findPosition(xpos, data.x)
        else:
            xmin = xpos
        xmax = xmin + 1
    if ypos:
        if physUnit:
            ymin = findPosition(ypos, data.y)
        else:
            ymin = ypos
        ymax = ymin + 1
    if zpos:
        if physUnit:
            zmin = findPosition(zpos, data.z)
        else:
            zmin = zpos
        zmax = zmin + 1

    try:
        assert(xmax <= data.x.shape[0])
    except AssertionError:
        print('Error: the position provided does not fit the data')
        print('xmax > dimx(data)')
        sys.exit(42)
    try:
        assert(ymax <= data.y.shape[0])
    except AssertionError:
        print('Error: the position provided does not fit the data')
        print('ymax > dimy(data)')
        sys.exit(42)
    try:
        assert(zmax <= data.z.shape[0])
    except AssertionError:
        print('Error: the position provided does not fit the data')
        print('zmax > dimz(data)')
        sys.exit(42)

    # Define array of figure
    fig, ax = pl.subplots(len(var),1, figsize=(6,3.5*len(var)))
    if not isinstance(ax, np.ndarray): ax = np.array([ax])
    for i, cvar in enumerate(var):
        if cvar == "rho":  
            title = r"$\rho$"
        elif "rhou" in cvar: 
            title = r"$\rho v_" + cvar[-1] + "$"
        elif "B" in cvar:    
            title = "B$_" + cvar[-1] + "$"
        else:
            title = cvar
        if log: title = "log(" + title + ")"
        ax[i].set_title(title)
        if (i == len(var) - 1): ax[i].set_xlabel(xlab)
        
        datPlot = dictDat[cvar][xmin:xmax,ymin:ymax,zmin:zmax].reshape(ldim)

        if physTick:
            if log:
                ax[i].semilogy(xplot, datPlot, *args, **kwargs)
            else:
                ax[i].plot(xplot, datPlot, *args, **kwargs)
        else:
            if log:
                ax[i].semilogy(datPlot, *args, **kwargs)
            else:
                ax[i].plot(datPlot, *args, **kwargs)
        ax[i].set_xlim(xlim)
        ax[i].xaxis.set_major_locator(mp.ticker.MaxNLocator(5))
        ax[i].yaxis.set_major_locator(mp.ticker.MaxNLocator(4))
    fig.tight_layout()

    if save:
        basename = "snapshot_"
        for format_ in ['png', 'pdf', 'eps']:
            pl.savefig(basename + "".join(name[:3].title() for name in var) \
                           + "." + format_)
    else:
        pl.show()

## plot2d function to plot a slice in a DUMSES-Hybrid simulation
# @param data a DumsesData object
# @param var a variable or a list of variables to plot. Variables can be a combination of ['rho', 'rhou{x,y,z}', 'E', 'B{x,y,z}'] and any simple mathematical operations.
# @param xpos position of the slice in the x-direction
# @param ypos position of the slice in the y-direction
# @param zpos position of the slice in the z-direction
# @param xmin lower x-boundary of the slice
# @param xmax upper x-boundary of the slice
# @param ymin lower y-boundary of the slice
# @param ymax upper y-boundary of the slice
# @param zmin lower z-boundary of the slice
# @param zmax upper z-boundary of the slice
# @param physUnit if \c True, use physical units for positions and boundaries
# @param physTick if \c True, use physical units for plotting
# @param polar if \c True, plot in polar coordinates
# @param log if \c True, plot in logscale
# @param save if \c True, save the plot in PNG, PDF and EPS
# @param cm colormap
# @param showTitle if \c True, display a title for each plot
def plot2d(data, var='rho', xpos=None, ypos=None, zpos=0 \
         , xmin=None, xmax=None, ymin=None, ymax=None, zmin=None, zmax=None \
         , physUnit=True, physTick=True, polar=False, log=False, save=False \
         , cm=None, showTitle=True, *args, **kwargs):
    """plot2d function to plot a slice in a DUMSES-Hybrid simulation

    Usage:
      plot2d(data [, var, xpos, ypos, zpos, xmin, xmax, xmin, xmax, zmin, zmax, physUnit, physTick, polar, log, save, cm, *args, **kwargs])
    With:
      data      a DumsesData object
      var       a variable or list of variables to plot. Variables can be a combination of ['rho', 'rhou{x,y,z}', 'E', 'B{x,y,z}'] and any simple mathematical operations. Default: 'rho'
      xpos      position of the slice in the x-direction. Default: None
      ypos      position of the slice in the y-direction. Default: None
      zpos      position of the slice in the z-direction. Default: 0
      xmin      lower x-boundary of the slice. Default: None
      xmax      upper x-boundary of the slice. Default: None
      ymin      lower y-boundary of the slice. Default: None
      ymax      upper y-boundary of the slice. Default: None
      zmin      lower z-boundary of the slice. Default: None
      zmax      upper z-boundary of the slice. Default: None
      physUnit  if True, use physical units for positions and boundaries. Default: True
      physTick  if True, use physical units for plotting. Default: True
      polar     if True, plot in polar coordinates. Default: False
      log       if True, plot in logscale. Default: False
      save      if True, save the plot in PNG, PDF and EPS. Default: False
      cm        colormap. Default: None
      showTitle if True, display title for each plot. Default: True
    """
    # Ensure that data is a DumsesData object
    try:
        assert(isinstance(data, dp.DumsesData))
    except AssertionError:
        print("Error: 'data' must be a DumsesData object!")
        sys.exit(42)

    dictDat = dp.dictionarize(data)
    
    if var == "all":
        var = dictDat.keys()
        if "Brx" in var:
            var.remove("Brx")
            var.remove("Bry")
            var.remove("Brz")
    if not isinstance(var, list): var = [var]

    # Parse the variables
    vname = re.compile(r'(?!\d)([\w_]+)')
    code  = []
    for cvar in var:
        st  = parser.expr(vname.sub(r"dictDat['\1']", cvar))
        code.append(st.compile())

    # Determine the position and the boundaries of the slice to plot
    if xpos != None: zpos = None
    if ypos != None: zpos = None
    try:
        assert([xpos, ypos, zpos].count(None) == 2)
    except AssertionError:
        print('Error: you have to specify one (and only one) direction and its position where to cut for the snapshot.')
        sys.exit(42)

    if polar and (xpos != None or ypos != None):
        print("Warning: you probably want to make your polar plot in a z-plane")

    if xpos == None:
        if physUnit:
            if xmin: xmin = findPosition(xmin, data.x)
            if xmax: xmax = findPosition(xmax, data.x)
        if not xmin: xmin = 0
        if not xmax: xmax = dictDat[dictDat.keys()[0]].shape[0]
    else:
        if physUnit:
            xmin = findPosition(xpos, data.x)
            xmax = xmin + 1
        else:
            xmin = xpos
            xmax = xpos + 1
    if ypos == None:
        if physUnit:
            if ymin: ymin = findPosition(ymin, data.y)
            if ymax: ymax = findPosition(ymax, data.y)
        if not ymin: ymin = 0
        if not ymax: ymax = dictDat[dictDat.keys()[0]].shape[1]
    else:
        if physUnit:
            ymin = findPosition(ypos, data.y)
            ymax = ymin + 1
        else:
            ymin = ypos
            ymax = ypos + 1
    if zpos == None:
        if physUnit:
            if zmin: zmin = findPosition(zmin, data.z)
            if zmax: zmax = findPosition(zmax, data.z)
        if not zmin: zmin = 0
        if not zmax: zmax = dictDat[dictDat.keys()[0]].shape[2]
    else:
        if physUnit:
            zmin = findPosition(zpos, data.z)
            zmax = zmin + 1
        else:
            zmin = zpos
            zmax = zpos + 1
    ldim = np.array([xmax - xmin, ymax - ymin, zmax - zmin])
    ldim = ldim[ldim > 1]
    try:
        assert(ldim.size == 2)
    except AssertionError:
        print('Error: Wrong dimensions! Check your {x,y,z} boundaries.')
        sys.exit(42)
    xdim, ydim = ldim

    # Define labels and extent of the plot
    if zmax - zmin == 1: 
        xlab = "x"
        ylab = "y"
        if physTick:
            extent = [data.x[xmin], data.x[xmax-1], data.y[ymin] \
                          , data.y[ymax-1]]
        else:
            extent = [xmin, xmax-1, ymin, ymax-1]
    if ymax - ymin == 1: 
        xlab = "x"
        ylab = "z"
        if physTick:
            extent = [data.x[xmin], data.x[xmax-1], data.z[zmin] \
                          , data.z[zmax-1]]
        else:
            extent = [xmin, xmax-1, zmin, zmax-1]
    if xmax - xmin == 1: 
        xlab = "y"
        ylab = "z"
        if physTick:
            extent = [data.y[ymin], data.y[ymax-1], data.z[zmin] \
                          , data.z[zmax-1]]
        else:
            extent = [ymin, ymax-1, zmin, zmax-1]

    if polar:
        xlab = "R"
        ylab = ""
        if not physTick: print("Warning: you probably want to make your polar plot in physical units")
    
    # Define array of figures
    fig, ax = pl.subplots(1, len(var), figsize=(6.*len(var), 4.))
    pl.subplots_adjust(bottom=0.15, left=0.15, wspace=0.5)
    if not isinstance(ax, np.ndarray): ax = np.array([ax])

    for i, cvar in enumerate(var):
        # Use SymPy to simplify and TeXify the title
        if "sympy" in sys.modules:
            title = sympy.simplify(cvar.replace('rhou', 'rho*v').replace('S', 'S_'))
            title = sympy.latex(title)
            title = '$' + title.replace('x', r'_x').replace('y', r'_y').replace('z', r'_z').replace('S_', 'S') + '$'
        else:
            if cvar == "rho":  
                title = r"$\rho$"
            elif "rhou" in cvar and len(cvar) == 5: 
                title = r"$\rho v_" + cvar[-1] + "$"
            elif "B" in cvar and len(cvar) == 2:    
                title = "B$_" + cvar[-1] + "$"
            else:
                title = cvar
        if log: title = "log(" + title + ")"

        # Compute the actual variable to plot
        res = eval(code[i])
        datPlot = res[xmin:xmax,ymin:ymax,zmin:zmax].reshape((xdim,ydim))
        if not polar:
            datPlot = datPlot.transpose()
        if showTitle:
            ax[i].set_title(title)
        if i == 0: ax[i].set_ylabel(ylab)
        ax[i].set_xlabel(xlab)
        # if polar, defines plot coordinates
        if polar:
            ax[i].yaxis.set_ticks_position('none')
            ax[i].set_yticks(())
            rtile = np.ones(ydim)
            ttile = np.ones(xdim)
            r     = np.linspace(extent[2], extent[3], xdim)
            theta = np.linspace(extent[0], extent[1], ydim)
            r     = r[:, np.newaxis]*rtile
            theta = theta*ttile[:, np.newaxis]
            xp    = r*np.cos(theta)
            yp    = r*np.sin(theta)
        # if in logscale, hatched negative regions
        if log:
            if polar:
                map1 = ax[i].pcolormesh(xp, yp, np.log10(abs(datPlot)) \
                                    , cmap=cm, *args, **kwargs)
                cs = ax[i].contourf(xp, yp, np.where(datPlot > 0, 1, -1), 2 \
                                    , hatches=['//', ''], alpha=0.)
            else:
                map1 = ax[i].imshow(np.log10(abs(datPlot)), extent=extent \
                                    , origin='lower', cmap=cm, *args, **kwargs)
                cs = ax[i].contourf(np.where(datPlot > 0, 1, -1), 2 \
                                , extent=extent, hatches=['//', ''], alpha=0.)
        else:
            if polar:
                map1 = ax[i].pcolormesh(xp, yp, datPlot, cmap=cm \
                                            , *args, **kwargs)
            else:
                map1 = ax[i].imshow(datPlot, extent=extent, origin='lower' \
                                        , cmap=cm, *args, **kwargs)
        ax[i].xaxis.set_major_locator(mp.ticker.MaxNLocator(4))
        if not polar:
            ax[i].yaxis.set_major_locator(mp.ticker.MaxNLocator(4))
        fig.colorbar(map1, ax=ax[i])
        
    if log:
        pl.suptitle('hatched where $<$ 0', fontsize=10, x=0.2)

    if save:
        basename = "snapshot_" + ("polar_" if polar else "")
        for format_ in ['png', 'pdf', 'eps']:
            pl.savefig(basename + "_".join(name[:3].title() \
               if name in dictDat.keys() else name.translate(None, ' +-*/()') \
               for name in var) + "." + format_)
        pl.clf()
    else:
        pl.show()

## findPostion function, to find the index of a given position in an array of positions
def findPosition(xphy, xpos):
    low = xpos[xpos < xphy]
    if low.size > 0:
        lbound = low[-1]
    else:
        lbound = xpos[0]
    high = xpos[xpos > xphy]
    
    if high.size > 0:
        hbound = high[0]
    else:
        hbound = xpos[-1]
    distLow  = abs(lbound - xphy)
    distHigh = abs(hbound - xphy)
    if distLow == distHigh or distLow < distHigh:
        returnPos = np.where(xpos == lbound)[0][0]
    else:
        returnPos = np.where(xpos == hbound)[0][0]
    return returnPos
