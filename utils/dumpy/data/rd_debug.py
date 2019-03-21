#===============================================================================
## \file rd_debug.py
# \brief
# \b DUMSES-Hybrid:
# These are data classes for DUMSES debug outputs
# \author
# Marc Joos <marc.joos@cea.fr>, Sebastien Fromang
# \copyright
# Copyrights 2013-2015, CEA.
# This file is distributed under the CeCILL-A & GNU/GPL licenses, see
# <http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html> and
# <http://www.gnu.org/licenses/>
# \date
# \b created:          06-03-2014 
# \b last \b modified: 05-25-2015

#===============================================================================
import sys
import tables
import numpy as np
import pylab as pl
import netCDF4 as nc

## DUMSES debugging data class
class DumsesDebug:
    ## Initialization method
    def __init__(self):
        ## directory of the outputs
        self.dir = ""
        ## list of output files
        self.filename = []
        ## number of cells of the sub-domain in the x-direction
        self.nx = 0
        ## number of cells of the sub-domain in the y-direction
        self.ny = 0
        ## number of cells of the sub-domain in the z-direction
        self.nz = 0
        ## number of MPI processes in the x-direction
        self.nxslice = 0
        ## number of MPI processes in the y-direction
        self.nyslice = 0
        ## number of MPI processes in the z-direction
        self.nzslice = 0
        ## dictionary containing debug variables
        self.dict = {}
        ## number of debug variables
        self.nvar = 0

    ## load method to read DUMSES debug data
    # @param idump number of the output to read
    # @param filedebug name of the debug file to read
    # @param filedir directory of the outputs
    # @param nghost number of ghost cells
    # @param io_type type of I/O in (\c binary, \c hdf5, \c phdf5, \c pnetcdf)
    def load(self, idump = 1, filedebug='fgodunov', filedir='./', nghost=3 \
                 , io_type='phdf5'):
        self.dir = filedir + "output_%06d/" %idump
        if io_type == 'phdf5':
            self.filename = self.dir + filedebug + ".000000.h5"

            f = tables.openFile(self.filename)

            meta_int = f.root.meta_int[:]
            nx = meta_int[0]; nxslice = meta_int[3]
            ny = meta_int[1]; nyslice = meta_int[4]
            nz = meta_int[2]; nzslice = meta_int[5]
            if nghost:
                nx -= 2*nghost
                if ny != 1:
                    ny -= 2*nghost
                if nz != 1:
                    nz -= 2*nghost
            nproc = nxslice*nyslice*nzslice
            self.nx = nx; self.nxslice = nxslice
            self.ny = ny; self.nyslice = nyslice
            self.nz = nz; self.nzslice = nzslice

            tempdict = {}
            for dset in f.root.__iter__():
                if dset.name != "meta_int":
                    tempdict[dset.name] = dset[:]
                    self.dict[dset.name] \
                        = np.zeros((nx*nxslice,ny*nyslice,nz*nzslice))
            f.close()
        elif io_type == 'pnetcdf':
            self.filename = self.dir + filedebug + ".000000.nc"
        
            f = nc.Dataset(self.filename)
        
            meta_int = f.variables['meta_int'][:]
            nx = meta_int[0]; nxslice = meta_int[3]
            ny = meta_int[1]; nyslice = meta_int[4]
            nz = meta_int[2]; nzslice = meta_int[5]
            if nghost:
                nx -= 2*nghost
                if ny != 1:
                    ny -= 2*nghost
                if nz != 1:
                    nz -= 2*nghost
            nproc = nxslice*nyslice*nzslice
            self.nx = nx; self.nxslice = nxslice
            self.ny = ny; self.nyslice = nyslice
            self.nz = nz; self.nzslice = nzslice
        
            tempdict = {}
            for dset in f.variables.iterkeys():
                if dset.name != "meta_int":
                    tempdict[dset] = f.variables[dset][:]
                    self.dict[dset] \
                        = np.zeros((nx*nxslice,ny*nyslice,nz*nzslice))
            f.close()
        else:
            raise NotImplementedError
                    
        selx = [i + nghost for i in range(nx)]
        sely = [i + nghost for i in range(ny)] if ny != 1 else [0]
        selz = [i + nghost for i in range(nz)] if nz != 1 else [0]

        nxs = nx + 2*nghost
        nys = ny + 2*nghost if ny != 1 else ny
        nzs = nz + 2*nghost if nz != 1 else nz

        for dset in tempdict.iterkeys():
            tempdict[dset] = np.reshape(tempdict[dset] \
                , (nproc, nzs, nys, nxs))
            tempdict[dset] = np.transpose(tempdict[dset], (0, 3, 2, 1))

        for n in range(nproc):
            xpos = n%nxslice
            ypos = n/(nxslice*nzslice)%nyslice
            zpos = (n/nxslice)%nzslice

            i0 = xpos*nx; i1 = (xpos + 1)*nx
            j0 = ypos*ny; j1 = (ypos + 1)*ny
            k0 = zpos*nz; k1 = (zpos + 1)*nz

            for dset in tempdict.iterkeys():
                self.dict[dset][i0:i1,j0:j1,k0:k1] \
                    = tempdict[dset][n][selx,:,:][:,sely,:][:,:,selz]

        self.nvar = self.dict.__len__()

    ## plot method for DUMSES debug data
    # @param var variable(s) to plot; can be a single string or a list of string. Default: \c None (plot all variables)
    # @param direction direction of the slicing. Default: \c "z"
    # @param prefix prefix for variable labels
    # @param suffix suffix for variable labels
    def plot(self, var=None, direction="z", prefix="", suffix="" \
                 , *args, **kwargs):
        if direction == "x":
            srange = [0,1,0,self.ny*self.nyslice,0,self.nz*self.nzslice]
        elif direction == "y":
            srange = [0,self.nx*self.nxslice,0,1,0,self.nz*self.nzslice]
        else:
            srange = [0,self.nx*self.nxslice,0,self.ny*self.nyslice,0,1]
        if not var:
            pl.figure(figsize=(self.nvar*4,4))
            pl.subplots_adjust(bottom=0.15)
            for i in range(self.nvar):
                ivar = self.dict.keys()[i]
                if len(ivar) == 2:
                    title = ivar[0] + "$_" + ivar[1] + "$"
                elif len(ivar) == 3:
                    if ivar == 'rho':
                        title = "$\\" + ivar + "$"
                    else:
                        title = ivar[0] + '$_{' + ivar[1:] + "}$"
                else:
                    title = ivar
                pl.subplot(1, self.nvar, i+1)
                pl.ylabel(prefix + title + suffix)
                if direction == "x":
                    pl.imshow(self.dict.get(ivar)[0,:,:], origin='lower' \
                                  , *args, **kwargs)
                elif direction == "y":
                    pl.imshow(self.dict.get(ivar)[:,0,:], origin='lower' \
                                  , *args, **kwargs)
                else:
                    pl.imshow(self.dict.get(ivar)[:,:,0], origin='lower' \
                                  , *args, **kwargs)
                pl.colorbar()
        elif isinstance(var, list):
            pl.figure(figsize=(len(var)*4,4))
            pl.subplots_adjust(bottom=0.15)
            i = 0
            for ivar in var:
                if len(ivar) == 2:
                    title = ivar[0] + "$_" + ivar[1] + "$"
                elif len(ivar) == 3:
                    if ivar == 'rho':
                        title = "$\\" + ivar + "$"
                    else:
                        title = ivar[0] + '$_{' + ivar[1:] + "}$"
                else:
                    title = ivar
                pl.subplot(1, len(var), i+1)
                pl.ylabel(prefix + title + suffix)
                if direction == "x":
                    pl.imshow(self.dict.get(ivar)[0,:,:], origin='lower' \
                                  , *args, **kwargs)
                elif direction == "y":
                    pl.imshow(self.dict.get(ivar)[:,0,:], origin='lower' \
                                  , *args, **kwargs)
                else:
                    pl.imshow(self.dict.get(ivar)[:,:,0], origin='lower' \
                                  , *args, **kwargs)
                pl.colorbar()
                i += 1
        else:
            pl.figure(figsize=(4,4))
            pl.subplots_adjust(bottom=0.15)
            if len(var) == 2:
                title = var[0] + "$_" + var[1] + "$"
            elif len(var) == 3:
                if var == 'rho':
                    title = "$\\" + var + "$"
                else:
                    title = var[0] + '$_{' + var[1:] + "}$"
            else:
                title = var
            pl.ylabel(prefix + title + suffix)
            if direction == "x":
                pl.imshow(self.dict.get(var)[0,:,:], origin='lower' \
                                  , *args, **kwargs)
            elif direction == "y":
                pl.imshow(self.dict.get(var)[:,0,:], origin='lower' \
                                  , *args, **kwargs)
            else:
                pl.imshow(self.dict.get(var)[:,:,0], origin='lower' \
                                  , *args, **kwargs)
            pl.colorbar()
        pl.show()
