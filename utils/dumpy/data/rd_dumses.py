#===============================================================================
## \file rd_dumses.py
# \brief
# \b DUMSES-Hybrid:
# These are data classes for DUMSES outputs
# \author
# Marc Joos <marc.joos@cea.fr>, Sebastien Fromang
# \copyright
# Copyrights 2013-2015, CEA.
# This file is distributed under the CeCILL-A & GNU/GPL licenses, see
# <http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html> and
# <http://www.gnu.org/licenses/>
# \date
# \b created:          04-15-2013 
# \b last \b modified: 06-25-2015

#===============================================================================
import sys
import glob
import numpy as np
import h5py  as h5
import netCDF4 as nc
from .formatChecker import formatChecker

## DUMSES data class factory
class DumsesDataFactory(object):
    metadataMakers = {}
    dataMakers = {}

    ## Wrapper for metadata class
    @classmethod
    def registerMetadata(cls, typ):
        def wrapper(maker):
            cls.metadataMakers[typ] = maker
            return cls
        return wrapper

    ## Wrapper for data class
    @classmethod
    def registerData(cls, typ):
        def wrapper(maker):
            cls.dataMakers[typ] = maker
            return cls
        return wrapper

    ## Maker for metadata class
    @classmethod
    def makeMetadata(cls, typ, *args, **kwargs):
        try:
            return cls.metadataMakers[typ](*args, **kwargs)
        except KeyError:
            print('FiletypeError\n' + typ + ': This filetype does not exist!')

    ## Maker for data class
    @classmethod
    def makeData(cls, typ, *args, **kwargs):
        try:
            return cls.dataMakers[typ](*args, **kwargs)
        except KeyError:
            print('FiletypeError\n' + typ + ': This filetype does not exist!')

## DUMSES data general class
class DumsesData(object):
    """DumsesData class:
    Use this class to read your DUMSES outputs.

    Usage:
      data = DumsesData()
      data.load(idump [, dir_, nghost, ndim])
    With:
      idump:              index of the output to read
      dir_ (optional):    directory containing the outputs
      nghost (optional):  the number of ghost cells
      ndim (optional):    the actual number of dimension of the problem, as set in the Makefile. Use it only if 'NDIM==3' for a 1D or 2D run.

    dumses_mpi outputs should be read with no problem using DumsesData; in this case, the contiguity of the data is automatically checked."""
    ## Initialization method
    def __init__(self):
        ## directory of the current output
        self.dir = ""
        ## list of metadata attributes
        self._metaData = ['nx', 'ny', 'nz', 'nxslice', 'nyslice', 'nzslice' \
            , 'nxglob', 'nyglob', 'nzglob', 'dim', 'dim_glob', 'nproc' \
            , 'time', 'dt', 'dx', 'dy', 'dz', 'xmin', 'xmax', 'ymin', 'ymax' \
            , 'zmin', 'zmax', 'gamma', 'ciso']
        ## list of data attributes
        self._data = ['rho', 'rhou', 'E', 'B', 'Br']

    ## load method to read DUMSES outputs
    # @param idump number of the output to read
    # @param filedir directory of the outputs
    # @param nghost number of ghost cells
    # @param ndim actual dimension of the problem if the code was compiled with
    # \c NDIM==3 (useful for tests in 1D or 2D in other directions than x)
    def load(self, idump=1, dir_='./', nghost=3, ndim=None):
        """DumsesData load method.
        For more information, see DumsesData docstring."""
        ## Check output format
        old, ext = formatChecker(dir_, idump, verbose=False)
        filename = ('data.*' if not(old) else 'slices.*')
        self.dir = dir_ + 'output_%06d/' %idump
        lfiles   = glob.glob(self.dir + filename)
        nfiles   = len(lfiles)
        self._getMetadata(idump, dir_, mype=0 \
                        , ext=('oldhdf5' if ext == 'oldphdf5' else ext))
        if nfiles > 1 and ext == 'h5': ext = 'h5seq'
        self._getData(idump, dir_, ext=ext, nghost=nghost, ndim=ndim)

    ## _getMetadata method to load metadata from DUMSES output
    def _getMetadata(self, idump=1, dir_='./', mype=0, ext=''):
        metadataReader = DumsesDataFactory.makeMetadata(ext)
        metadataReader.load(idump, dir_, mype)
        for attr in self._metaData:
            self.__setattr__(attr, metadataReader.__getattribute__(attr))
        # If it is a dumses_hybrid file format or a paralle dumses_mpi file
        # format, reconstruct {x,y,z} from metadata
        if not(ext == 'oldbinary' or ext == 'oldhdf5'):
            self.x = np.linspace(self.xmin + self.dx/2, self.xmax - self.dx/2 \
                                     , self.nxglob)
            self.y = np.linspace(self.ymin + self.dy/2, self.ymax - self.dy/2 \
                                     , self.nyglob)
            self.z = np.linspace(self.zmin + self.dz/2, self.zmax - self.dz/2 \
                                     , self.nzglob)

    ## _getData method to load data from DUMSES output
    def _getData(self, idump=1, dir_='./', ext='', nghost=3, ndim=None):
        dataReader = DumsesDataFactory.makeData(ext, self.nx, self.ny, self.nz \
                     , self.nxslice, self.nyslice, self.nzslice \
                     , self.nxglob, self.nyglob, self.nzglob, self.nproc)
        dataReader.load(idump=idump, dir_=dir_, nghost=nghost, ndim=ndim)
        for attr in self._data:
            self.__setattr__(attr, dataReader.__getattribute__(attr))
        # If it is a dumses_mpi sequential file format, read {x,y,z}
        if ext == 'oldbinary' or ext == 'oldhdf5':
            self.x = dataReader.x; self.y = dataReader.y; self.z = dataReader.z
            self.xmin = dataReader.xmin - self.dx/2
            self.xmax = dataReader.xmax + self.dx/2
            self.ymin = dataReader.ymin - self.dy/2
            self.ymax = dataReader.ymax + self.dy/2
            self.zmin = dataReader.zmin - self.dz/2
            self.zmax = dataReader.zmax + self.dz/2

    ## getPrimitive method to compute primitive variables
    # @param iso set to 'True' if isothermal computation
    def getPrimitive(self, iso=False):
        """getPrimitive method, to compute primitive variables.

        Usage:
          data.getPrimitive([iso])
        With:
          iso: set to True if isothermal computation
        """
        nxglob, nyglob, nzglob = self.nxglob, self.nyglob, self.nzglob
        self.q_rho = self.rho
        self.q_u   = self.rhou/self.rho[:,:,:,None]
        q_Bx  = (self.B[1:,1:,1:,0] \
                   + self.B[:nxglob-1,:nyglob-1,:nzglob-1,0])/2.
        q_By  = (self.B[1:,1:,1:,1] \
                   + self.B[:nxglob-1,:nyglob-1,:nzglob-1,1])/2.
        q_Bz  = (self.B[1:,1:,1:,2] \
                   + self.B[:nxglob-1,:nyglob-1,:nzglob-1,2])/2.
        self.q_B = np.array([q_Bx, q_By, q_Bz]).transpose(1,2,3,0)
        if iso:
            self.q_E = self.q_rho*self.ciso**2
        else:
            Eken = (self.q_u[:,:,:,0]**2 + self.q_u[:,:,:,1]**2 \
                              + self.q_u[:,:,:,2]**2)/2.
            Emag = (self.q_B[:,:,:,0]**2 + self.q_B[:,:,:,1]**2 \
                              + self.q_B[:,:,:,2]**2)/2.
            Eint = (self.E[:nxglob-1,:nyglob-1,:nzglob-1] \
                  - Emag[:nxglob-1,:nyglob-1,:nzglob-1]) \
                  /self.rho[:nxglob-1,:nyglob-1,:nzglob-1] \
                  - Eken[:nxglob-1,:nyglob-1,:nzglob-1]
            self.q_E = (self.gamma - 1.) \
                *self.q_rho[:nxglob-1,:nyglob-1,:nzglob-1]*Eint

    ## get_1d method, to extract a 1D array out of DumsesData data
    # @param var:       name of the variable to extract, in ['rho', 'vx', 'vy', 'vz', 'Bx', 'By', 'Bz', 'E']
    # @param direction: direction in which the cut is done
    # @param pos:       position in the two other directions of the cut
    def get_1d(self, var='rho', direction='x', pos=[0,0]):
        """get_1d method, to extract a 1D array out of DumsesData data.

        Usage:
          data.get_1d([var, direction, pos])
        With:
          var:       name of the variable to extract, in ['rho', 'vx', 'vy', 'vz', 'Bx', 'By', 'Bz', 'E']. Default: 'rho'
          direction: direction in which the cut is done. Default: 'x'
          pos:       position in the two other directions of the cut. Default: [0,0]
        """
        ypos, zpos = pos

        try:
            assert(var in ['rho', 'vx', 'vy', 'vz', 'Bx', 'By', 'Bz', 'E'])
        except AssertionError:
            print("'var' must be in ['rho', 'vx', 'vy', 'vz', 'Bx', 'By', 'Bz', 'E']")
            sys.exit(42)

        if direction == 'x':
            try:
                assert(0 <= ypos < self.nyglob)
            except AssertionError:
                print('pos[0] (along y) is out of bounds!')
                sys.exit(42)
            try:
                assert(0 <= zpos < self.nzglob)
            except AssertionError:
                print('pos[1] (along z) is out of bounds!')
                sys.exit(42)
            xmin = 0
            xmax = self.nxglob
            ymin = ypos
            ymax = ypos + 1
            zmin = zpos
            zmax = zpos + 1
        elif direction == 'y':
            try:
                assert(0 <= ypos < self.nxglob)
            except AssertionError:
                print('pos[0] (along x) is out of bounds!')
                sys.exit(42)
            try:
                assert(0 <= zpos < self.nzglob)
            except AssertionError:
                print('pos[1] (along z) is out of bounds!')
                sys.exit(42)
            xmin = ypos    
            xmax = ypos + 1
            ymin = 0     
            ymax = self.nyglob
            zmin = zpos
            zmax = zpos + 1
        else:
            try:
                assert(0 <= ypos < self.nxglob)
            except AssertionError:
                print('pos[0] (along x) is out of bounds!')
                sys.exit(42)
            try:
                assert(0 <= zpos < self.nyglob)
            except AssertionError:
                print('pos[1] (along y) is out of bounds!')
                sys.exit(42)
            xmin = ypos    
            xmax = ypos + 1
            ymin = zpos    
            ymax = zpos + 1
            zmin = 0     
            zmax = self.nzglob
        ldim = np.array([xmax - xmin, ymax - ymin, zmax - zmin])
        ldim = ldim[ldim > 1]
        
        if var[0] == 'v': var = 'rhou' + var[-1]
        if var[-1] in ['x', 'y', 'z']:
            if var[-1] == 'x':
                extract = self.__getattribute__(var[:-1])[xmin:xmax,ymin:ymax,zmin:zmax,0]
            if var[-1] == 'y':
                extract = self.__getattribute__(var[:-1])[xmin:xmax,ymin:ymax,zmin:zmax,1]
            if var[-1] == 'z':
                extract = self.__getattribute__(var[:-1])[xmin:xmax,ymin:ymax,zmin:zmax,2]
            if var[:-1] == 'rhou':
                extract = extract/self.rho[xmin:xmax,ymin:ymax,zmin:zmax]
        else:
            extract = self.__getattribute__(var)[xmin:xmax,ymin:ymax,zmin:zmax]

        return extract.reshape(ldim)

    ## get_2d method, to extract a 2D array out of DumsesData data
    # @param var:       name of the variable to extract, in ['rho', 'vx', 'vy', 'vz', 'Bx', 'By', 'Bz', 'E']
    # @param direction: direction in which the cut is done
    # @param pos:       position in the given position
    def get_2d(self, var='rho', direction='z', pos=0):
        """get_2d method, to extract a 2D array out of DumsesData data.

        Usage:
          data.get_2d([var, direction, pos])
        With:
          var:       name of the variable to extract, in ['rho', 'vx', 'vy', 'vz', 'Bx', 'By', 'Bz', 'E']. Default: 'rho'
          direction: direction in which the cut is done. Default: 'z'
          pos:       position in the given direction. Default: 0
        """
        xmin = 0; xmax = self.nxglob
        ymin = 0; ymax = self.nyglob
        zmin = 0; zmax = self.nzglob

        try:
            assert(var in ['rho', 'vx', 'vy', 'vz', 'Bx', 'By', 'Bz', 'E'])
        except AssertionError:
            print("'var' must be in ['rho', 'vx', 'vy', 'vz', 'Bx', 'By', 'Bz', 'E']")

        if direction == 'x':
            try:
                assert(0 <= pos < self.nxglob)
            except AssertionError:
                print('position along x is out of bounds!')
                sys.exit(42)
            xmin = pos
            xmax = pos + 1
        elif direction == 'y':
            try:
                assert(0 <= pos < self.nyglob)
            except AssertionError:
                print('position along y is out of bounds!')
                sys.exit(42)
            ymin = pos
            ymax = pos + 1
        else:
            try:
                assert(0 <= pos < self.nzglob)
            except AssertionError:
                print('position along z is out of bounds!')
                sys.exit(42)
            zmin = pos
            zmax = pos + 1
        ldim = np.array([xmax - xmin, ymax - ymin, zmax - zmin])
        ldim = ldim[ldim > 1]

        if var[0] == 'v': var = 'rhou' + var[-1]
        if var[-1] in ['x', 'y', 'z']:
            if var[-1] == 'x':
                extract = self.__getattribute__(var[:-1])[xmin:xmax,ymin:ymax,zmin:zmax,0]
            if var[-1] == 'y':
                extract = self.__getattribute__(var[:-1])[xmin:xmax,ymin:ymax,zmin:zmax,1]
            if var[-1] == 'z':
                extract = self.__getattribute__(var[:-1])[xmin:xmax,ymin:ymax,zmin:zmax,2]
            if var[:-1] == 'rhou':
                extract = extract/self.rho[xmin:xmax,ymin:ymax,zmin:zmax]
        else:
            extract = self.__getattribute__(var)[xmin:xmax,ymin:ymax,zmin:zmax]
        return extract.reshape(ldim)

## Metadata general class
class MetadataReader(object):
    def __init__(self):
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
        ## total number of cells in the x-direction
        self.nxglob = 0
        ## total number of cells in the y-direction
        self.nyglob = 0
        ## total number of cells in the z-direction
        self.nzglob = 0
        ## dimension of a sub-domain (nx, ny, nz)
        self.dim = 0
        ## global dimension (nxglob, nyglob, nzglob)
        self.dim_glob = 0
        ## number of MPI processes
        self.nproc = 0
        ## time of the output
        self.time = 0.
        ## timestep of the output
        self.dt = 0.
        ## resolution in the x-direction
        self.dx = 0.
        ## resolution in the y-direction
        self.dy = 0.
        ## resolution in the z-direction
        self.dz = 0.
        ## \f$x_{min}\f$
        self.xmin = 0.
        ## \f$x_{max}\f$
        self.xmax = 0.
        ## \f$y_{min}\f$
        self.ymin = 0.
        ## \f$y_{max}\f$
        self.ymax = 0.
        ## \f$z_{min}\f$
        self.zmin = 0.
        ## \f$z_{max}\f$
        self.zmax = 0.
        ## \f$\gamma\f$
        self.gamma = None
        ## \f$c_{iso}\f$
        self.ciso = None

    ## load method; returns filename of the current outputfile
    def load(self, idump=1, dir_='./', mype=0, ext=''):
        dir_ = dir_ + 'output_%06d/' %idump
        fname = ('data.%06d' + '.' + ext if ext != '' else 'slices.%06d') %mype
        return dir_ + fname

    ## read method to set all metadata attributes
    def read(self, meta_int, meta_real):
        self.nx = meta_int[0]; self.nxslice = meta_int[3]
        self.ny = meta_int[1]; self.nyslice = meta_int[4]
        self.nz = meta_int[2]; self.nzslice = meta_int[5]
        self.dim = meta_int[0:3]; self.dim_glob = meta_int[0:3]*meta_int[3:6]
        self.nxglob, self.nyglob, self.nzglob = self.dim_glob[:]
        self.nproc = self.nxslice*self.nyslice*self.nzslice

        self.time = meta_real[0]; self.dt = meta_real[1]
        self.dx = meta_real[2]; self.dy = meta_real[3]; self.dz = meta_real[4]
        self.xmin = meta_real[5]; self.xmax = meta_real[6]
        self.ymin = meta_real[7]; self.ymax = meta_real[8]
        self.zmin = meta_real[9]; self.zmax = meta_real[10]
        self.gamma = meta_real[11]
        self.ciso  = meta_real[12]

## Metadata binary class; reads metadata from Fortran binary outputs
@DumsesDataFactory.registerMetadata('bin')
class MetadataBinaryReader(MetadataReader):
    def __init__(self):
        MetadataReader.__init__(self)

    ## load method; actually reads metadata in file
    def load(self, idump=1, dir_='./', mype=0, ext='bin'):
        filename = MetadataReader.load(self, idump, dir_, mype, ext)

        f = open(filename, 'rb')
        meta_int  = getArray(f, 6, 'i4')
        meta_real = np.zeros(13)
        meta_real[:5]   = getArray(f, 5, 'f8')
        meta_real[5:11] = getArray(f, 6, 'f8')
        meta_real[11:]  = getArray(f, 2, 'f8')
        f.close()
        
        self.read(meta_int, meta_real)

## Metadata HDF5 class; reads metadata from HDF5 outputs
@DumsesDataFactory.registerMetadata('h5')
class MetadataHDF5Reader(MetadataReader):
    def __init__(self):
        MetadataReader.__init__(self)

    ## load method; actually reads metadata in file
    def load(self, idump=1, dir_='./', mype=0, ext='h5'):
        filename = MetadataReader.load(self, idump, dir_, mype, ext)

        f = h5.File(filename, 'r')
        meta_int  = f['meta_int'][:]
        meta_real = f['meta_real'][:]
        f.close()

        self.read(meta_int, meta_real)

## Metadata NetCDF class; reads metadata from NetCDF outputs
@DumsesDataFactory.registerMetadata('nc')
class MetadataNetCDFReader(MetadataReader):
    def __init__(self):
        MetadataReader.__init__(self)

    ## load method; actually reads metadata in file
    def load(self, idump=1, dir_='./', mype=0, ext='nc'):
        filename = MetadataReader.load(self, idump, dir_, mype, ext)

        f = nc.Dataset(filename, 'r')
        meta_int  = f.variables['meta_int'][:]
        meta_real = f.variables['meta_real'][:]
        f.close()

        self.read(meta_int, meta_real)

## Metadata old binary class: read metadata from dumses_mpi outputs
@DumsesDataFactory.registerMetadata('oldbinary')
class MetadataOldHDF5Reader(MetadataReader):
    def __init__(self):
        MetadataReader.__init__(self)

    ## load method
    def load(self, idump=1, dir_='./', mype=0, ext=''):
        filename = MetadataReader.load(self, idump, dir_, mype, ext)

        f = open(filename, 'rb')
        meta_real = getArray(f, 5, 'f8')
        tmp = getArray(f, 3, 'i4')
        meta_int = getArray(f, 3, 'i4')
        meta_int = np.append(meta_int, getArray(f, 3, 'i4'))
        meta_real = np.append(meta_real, [None for i in xrange(8)])
        f.close()

        self.read(meta_int, meta_real)

## Metadata old HDF5 class: read metadata from dumses_mpi outputs
@DumsesDataFactory.registerMetadata('oldhdf5')
class MetadataOldHDF5Reader(MetadataReader):
    def __init__(self):
        MetadataReader.__init__(self)

    ## load method
    def load(self, idump=1, dir_='./', mype=0, ext=''):
        filename = MetadataReader.load(self, idump, dir_, mype, ext)

        f = h5.File(filename, 'r')
        meta_int  = f['para_int'][:]
        meta_real = f['para_real'][:]
        if 'boxSize' in f.keys(): box_size  = f['boxSize'][:]
        meta_int[:6] = meta_int[3:9]
        if 'boxSize' in f.keys(): 
            meta_real = np.append(meta_real, box_size)
        else:
            meta_real = np.append(meta_real, [None]*6)
        meta_real = np.append(meta_real, [None]*2)
        f.close()

        self.read(meta_int, meta_real)

## Metadata old NetCDF class: read metadata from dumses_mpi outputs
@DumsesDataFactory.registerMetadata('oldpnetcdf')
class MetadataOldNetCDFReader(MetadataReader):
    def __init__(self):
        MetadataReader.__init__(self)

    ## load method
    def load(self, idump=1, dir_='./', mype=0, ext=''):
        filename = MetadataReader.load(self, idump, dir_, mype, ext)

        f = nc.Dataset(filename, 'r')
        meta_int  = f.variables['para_int'][:]
        meta_real = f.variables['para_real'][:]
        box_size  = f.variables['boxSize'][:]
        meta_int[:6] = meta_int[3:9]
        meta_real = np.append(meta_real, box_size)
        meta_real = np.append(meta_real, [None]*2)
        f.close()

        self.read(meta_int, meta_real)

## Data general class
class DataReader(object):
    def __init__(self, nx=1, ny=1, nz=1, nxslice=1, nyslice=1, nzslice=1 \
                     , nxglob=1, nyglob=1, nzglob=1, nproc=1):
        ## density variable of size (nxglob, nyglob, nzglob)
        self.rho  = np.ndarray(shape=(nxglob,nyglob,nzglob))
        ## velocity variable of size (nxglob, nyglob, nzglob, 3)
        self.rhou = np.ndarray(shape=(nxglob,nyglob,nzglob,3))
        ## energy variable of size (nxglob, nyglob, nzglob)
        self.E    = np.ndarray(shape=(nxglob,nyglob,nzglob))
        ## magnetic field variable of size (nxglob, nyglob, nzglob, 3)
        self.B    = np.ndarray(shape=(nxglob,nyglob,nzglob,3))
        ## magnetic field variable of size (nxglob, nyglob, nzglob, 3)
        self.Br   = np.ndarray(shape=(nxglob,nyglob,nzglob,3))
        ## number of cells of the sub-domain in the x-direction
        self.nx = nx
        ## number of cells of the sub-domain in the y-direction
        self.ny = ny
        ## number of cells of the sub-domain in the z-direction
        self.nz = nz
        ## number of MPI processes in the x-direction
        self.nxslice = nxslice
        ## number of MPI processes in the y-direction
        self.nyslice = nyslice
        ## number of MPI processes in the z-direction
        self.nzslice = nzslice
        ## number of MPI processes
        self.nproc = nproc
        ## total number of cells in the x-direction
        self.nxglob = nxglob
        ## total number of cells in the y-direction
        self.nyglob = nyglob
        ## total number of cells in the z-direction
        self.nzglob = nzglob

    ## _getFile method; returns filename of the current outputfile
    def _getFile(self, idump=1, dir_='./', mype=0, ext=''):
        dir_ = dir_ + 'output_%06d/' %idump
        fname = ('data.%06d' + '.' + ext if ext != '' else 'slices.%06d') %mype
        return dir_ + fname

## Sequential data general class
class DataSequentialReader(DataReader):
    def __init__(self, nx, ny, nz, nxslice, nyslice, nzslice, nxglob, nyglob \
                     , nzglob, nproc):
        DataReader.__init__(self, nx, ny, nz, nxslice, nyslice, nzslice \
                                , nxglob, nyglob, nzglob, nproc)
    
    def _load(self, *args, **kwargs):
        pass

    ## load method; loops over all files and retrieve data
    def load(self, idump=1, dir_='./', nghost=3, ndim=None):
        x = np.ndarray(shape=self.nxglob)
        y = np.ndarray(shape=self.nyglob)
        z = np.ndarray(shape=self.nzglob)
        old = False
        for mype in range(self.nproc):
            data = self._load(idump, dir_, nghost, mype, ndim)

            xpos = mype%self.nxslice
            ypos = int(mype/(self.nxslice*self.nzslice)%self.nyslice)
            zpos = int(mype/(self.nxslice)%self.nzslice)

            i0 = xpos*self.nx; i1 = (xpos + 1)*self.nx
            j0 = ypos*self.ny; j1 = (ypos + 1)*self.ny
            k0 = zpos*self.nz; k1 = (zpos + 1)*self.nz

            a0, a1, b0, b1, c0, c1 = 0, 1, 0, 1, 0, 1
            if(self.nx != 1 or ndim > 1): a0 = nghost; a1 = nghost + self.nx
            if(self.ny != 1 or ndim > 1): b0 = nghost; b1 = nghost + self.ny
            if(self.nz != 1 or ndim > 1): c0 = nghost; c1 = nghost + self.nz
            
            self.rho[i0:i1,j0:j1,k0:k1]    = data['rho'][a0:a1,b0:b1,c0:c1]
            self.E[i0:i1,j0:j1,k0:k1]      = data['E'][a0:a1,b0:b1,c0:c1]
            self.rhou[i0:i1,j0:j1,k0:k1,0] = data['rhou'][a0:a1,b0:b1,c0:c1,0]
            self.rhou[i0:i1,j0:j1,k0:k1,1] = data['rhou'][a0:a1,b0:b1,c0:c1,1]
            self.rhou[i0:i1,j0:j1,k0:k1,2] = data['rhou'][a0:a1,b0:b1,c0:c1,2]
            self.B[i0:i1,j0:j1,k0:k1,0]    = data['B'][a0:a1,b0:b1,c0:c1,0]
            self.B[i0:i1,j0:j1,k0:k1,1]    = data['B'][a0:a1,b0:b1,c0:c1,1]
            self.B[i0:i1,j0:j1,k0:k1,2]    = data['B'][a0:a1,b0:b1,c0:c1,2]
            self.Br[i0:i1,j0:j1,k0:k1,0]   = data['Br'][a0:a1,b0:b1,c0:c1,0]
            self.Br[i0:i1,j0:j1,k0:k1,1]   = data['Br'][a0:a1,b0:b1,c0:c1,1]
            self.Br[i0:i1,j0:j1,k0:k1,2]   = data['Br'][a0:a1,b0:b1,c0:c1,2]

            if 'x' in data:
                old = True
                x[i0:i1] = data['x'][a0:a1]
                y[j0:j1] = data['y'][b0:b1]
                z[k0:k1] = data['z'][c0:c1]

        if old:
            self.x = x; self.y = y; self.z = z
            self.xmin = min(x); self.xmax = max(x)
            self.ymin = min(y); self.ymax = max(y)
            self.zmin = min(z); self.zmax = max(z)

## Sequential Fortran binary data class
@DumsesDataFactory.registerData('bin')
class DataBinaryReader(DataSequentialReader):
    def __init__(self, nx, ny, nz, nxslice, nyslice, nzslice, nxglob, nyglob \
                     , nzglob, nproc):
        DataSequentialReader.__init__(self, nx, ny, nz, nxslice, nyslice \
                                    , nzslice, nxglob, nyglob, nzglob, nproc)

    ## _load method; reads one Fortran binary output file
    def _load(self, idump=1, dir_='./', nghost=3, mype=0, ndim=None, ext='bin'):
        filename = self._getFile(idump, dir_, mype, ext)

        nvar = 11
        nxtot, nytot, nztot = self.nx, self.ny, self.nz
        if(self.nx != 1 or ndim > 1): nxtot = self.nx + 2*nghost
        if(self.ny != 1 or ndim > 1): nytot = self.ny + 2*nghost
        if(self.nz != 1 or ndim > 1): nztot = self.nz + 2*nghost
        shape = (nvar, nztot, nytot, nxtot)
        ntot  = nxtot*nytot*nztot*nvar
        
        f = open(filename, 'rb')
        # Offset is #meta_int*4bytes + #meta_real*8bytes + 2*#record*4bytes
        offset = 6*4 + 13*8 + 2*4*4
        f.seek(offset)
        uin = getArray(f, ntot, 'f8')
        f.close()
        uin = uin.reshape(shape)
        rho = uin[0,:,:,:].transpose(); E   = uin[4,:,:,:].transpose()
        vx  = uin[1,:,:,:].transpose(); Bx  = uin[5,:,:,:].transpose() 
        vy  = uin[2,:,:,:].transpose(); By  = uin[6,:,:,:].transpose() 
        vz  = uin[3,:,:,:].transpose(); Bz  = uin[7,:,:,:].transpose()
        Bxr = uin[8,:,:,:].transpose(); Byr = uin[9,:,:,:].transpose()
        Bzr = uin[10,:,:,:].transpose()
        return dict({'rho': rho, 'E': E \
                  , 'rhou': np.array([vx,vy,vz]).transpose(1,2,3,0) \
                  , 'B': np.array([Bx,By,Bz]).transpose(1,2,3,0) \
                  , 'Br': np.array([Bxr,Byr,Bzr]).transpose(1,2,3,0)})

## Sequential HDF5 data class
@DumsesDataFactory.registerData('h5seq')
class DataHDF5Reader(DataSequentialReader):
    def __init__(self, nx, ny, nz, nxslice, nyslice, nzslice, nxglob, nyglob \
                     , nzglob, nproc):
        DataSequentialReader.__init__(self, nx, ny, nz, nxslice, nyslice \
                                    , nzslice, nxglob, nyglob, nzglob, nproc)

    ## _load method; reads one HDF5 output file
    def _load(self, idump=1, dir_='./', nghost=3, mype=0, ndim=None, ext='h5'):
        filename = self._getFile(idump, dir_, mype, ext)

        f = h5.File(filename, 'r')
        rho = f['rho'][:].transpose(); E  = f['E'][:].transpose()
        vx  = f['vx'][:].transpose();  Bx = f['Bx'][:].transpose()
        vy  = f['vy'][:].transpose();  By = f['By'][:].transpose()
        vz  = f['vz'][:].transpose();  Bz = f['Bz'][:].transpose()
        Bxr = f['Bxr'][:].transpose(); Byr = f['Byr'][:].transpose()
        Bzr = f['Bzr'][:].transpose()
        f.close()
        return dict({'rho': rho, 'E': E \
                  , 'rhou': np.array([vx,vy,vz]).transpose(1,2,3,0) \
                  , 'B': np.array([Bx,By,Bz]).transpose(1,2,3,0) \
                  , 'Br': np.array([Bxr,Byr,Bzr]).transpose(1,2,3,0)})

## Old binary data class
@DumsesDataFactory.registerData('oldbinary')
class DataOldBinaryReader(DataSequentialReader):
    def __init__(self, nx, ny, nz, nxslice, nyslice, nzslice, nxglob, nyglob \
                     , nzglob, nproc):
        DataSequentialReader.__init__(self, nx, ny, nz, nxslice, nyslice \
                                  , nzslice, nxglob, nyglob, nzglob, nproc)
    
    ## _load method; reads one Fortran binary output file
    def _load(self, idump=1, dir_='./', nghost=3, mype=0, ndim=None, ext=''):
        filename = self._getFile(idump, dir_, mype, ext)

        nvar = 11
        nxtot, nytot, nztot = self.nx, self.ny, self.nz
        if(self.nx != 1 or ndim > 1): nxtot = self.nx + 2*nghost
        if(self.ny != 1 or ndim > 1): nytot = self.ny + 2*nghost
        if(self.nz != 1 or ndim > 1): nztot = self.nz + 2*nghost
        shape = (nvar, nztot, nytot, nxtot)
        ntot  = nxtot*nytot*nztot*nvar
        xtot  = nxtot + nytot + nztot
        
        f = open(filename, 'rb')
        # Offset is #meta_int*4bytes + #meta_real*8bytes + 2*#record*4bytes
        offset = 9*4 + 5*8 + 2*4*4
        f.seek(offset)
        # x, y, z have to be retrieve from the file
        positions = getArray(f, xtot, 'f8')
        x = positions[0:nxtot]
        y = positions[nxtot:nxtot+nytot]
        z = positions[nxtot+nytot:nxtot+nytot+nztot]
        uin = getArray(f, ntot, 'f8')
        f.close()
        uin = uin.reshape(shape)
        rho = uin[0,:,:,:].transpose(); E   = uin[4,:,:,:].transpose()
        vx  = uin[1,:,:,:].transpose(); Bx  = uin[5,:,:,:].transpose() 
        vy  = uin[2,:,:,:].transpose(); By  = uin[6,:,:,:].transpose() 
        vz  = uin[3,:,:,:].transpose(); Bz  = uin[7,:,:,:].transpose()
        Bxr = uin[8,:,:,:].transpose(); Byr = uin[9,:,:,:].transpose()
        Bzr = uin[10,:,:,:].transpose()
        return dict({'x': x, 'y': y, 'z': z, 'rho': rho, 'E': E \
                  , 'rhou': np.array([vx,vy,vz]).transpose(1,2,3,0) \
                  , 'B': np.array([Bx,By,Bz]).transpose(1,2,3,0) \
                  , 'Br': np.array([Bxr,Byr,Bzr]).transpose(1,2,3,0)})

## Old HDF5 data class
@DumsesDataFactory.registerData('oldhdf5')
class DataOldHDF5Reader(DataSequentialReader):
    def __init__(self, nx, ny, nz, nxslice, nyslice, nzslice, nxglob, nyglob \
                     , nzglob, nproc):
        DataSequentialReader.__init__(self, nx, ny, nz, nxslice, nyslice \
                                  , nzslice, nxglob, nyglob, nzglob, nproc)
    
    ## _load method; reads one HDF5 output file
    def _load(self, idump=1, dir_='./', nghost=3, mype=0, ndim=None, ext=''):
        filename = self._getFile(idump, dir_, mype, ext)

        nvar = 11
        nxtot, nytot, nztot = self.nx, self.ny, self.nz
        if(self.nx != 1 or ndim > 1): nxtot = self.nx + 2*nghost
        if(self.ny != 1 or ndim > 1): nytot = self.ny + 2*nghost
        if(self.nz != 1 or ndim > 1): nztot = self.nz + 2*nghost
        shape = (nvar, nztot, nytot, nxtot)

        f = h5.File(filename, 'r')
        # x, y, z have to be retrieve from the file
        x = f['x'][:]; y = f['y'][:]; z = f['z'][:]
        uin = f['uin'][:]
        f.close()
        uin = uin.reshape(shape)
        rho = uin[0,:,:,:].transpose(); E   = uin[4,:,:,:].transpose()
        vx  = uin[1,:,:,:].transpose(); Bx  = uin[5,:,:,:].transpose() 
        vy  = uin[2,:,:,:].transpose(); By  = uin[6,:,:,:].transpose() 
        vz  = uin[3,:,:,:].transpose(); Bz  = uin[7,:,:,:].transpose()
        Bxr = uin[8,:,:,:].transpose(); Byr = uin[9,:,:,:].transpose()
        Bzr = uin[10,:,:,:].transpose()
        return dict({'x': x, 'y': y, 'z': z, 'rho': rho, 'E': E \
                  , 'rhou': np.array([vx,vy,vz]).transpose(1,2,3,0) \
                  , 'B': np.array([Bx,By,Bz]).transpose(1,2,3,0) \
                  , 'Br': np.array([Bxr,Byr,Bzr]).transpose(1,2,3,0)})

## Parallel data general class
class DataParallelReader(DataReader):
    def __init__(self, nx, ny, nz, nxslice, nyslice, nzslice, nxglob, nyglob \
                     , nzglob, nproc):
        DataReader.__init__(self, nx, ny, nz, nxslice, nyslice, nzslice \
                                , nxglob, nyglob, nzglob, nproc)

    ## _load method to read chunk of data in a parallel file
    def _load(self, *args, **kwargs):
        pass

    ## _inline method to check contiguity of data. For backward compatibility.
    def _inline(self, *args, **kwargs):
        return False

    ## load method; loops over all the MPI processes to retrieve data
    def load(self, idump=1, dir_='./', nghost=3, ndim=None):
        notInline = self._inline(idump, dir_, 0)
        for mype in range(self.nproc):
            xpos = mype%self.nxslice
            ypos = mype/(self.nxslice*self.nzslice)%self.nyslice
            zpos = mype/(self.nxslice)%self.nzslice

            i0 = xpos*self.nx; i1 = (xpos + 1)*self.nx
            j0 = ypos*self.ny; j1 = (ypos + 1)*self.ny
            k0 = zpos*self.nz; k1 = (zpos + 1)*self.nz

            if notInline:
                a0, a1, b0, b1, c0, c1 = 0, 1, 0, 1, 0, 1
                if(self.nx != 1 or ndim > 1): 
                    a0 = xpos*(self.nx + 2*nghost) + nghost
                    a1 = (xpos + 1)*(self.nx + 2*nghost) - nghost
                if(self.ny != 1 or ndim > 1): 
                    b0 = ypos*(self.ny + 2*nghost) + nghost
                    b1 = (ypos + 1)*(self.ny + 2*nghost) - nghost
                if(self.nz != 1 or ndim > 1): 
                    c0 = zpos*(self.nz + 2*nghost) + nghost 
                    c1 = (zpos + 1)*(self.nz + 2*nghost) - nghost
            else:
                a0, a1, b0, b1, c0, c1 = 0, 1, 0, 1, 0, 1
                if(self.nx != 1 or ndim > 1): a0 = nghost; a1 = nghost + self.nx
                if(self.ny != 1 or ndim > 1): b0 = nghost; b1 = nghost + self.ny
                if(self.nz != 1 or ndim > 1): c0 = nghost; c1 = nghost + self.nz
                if(self.nz != 1): 
                    offset = mype*(self.nz + 2*nghost)
                else:
                    offset = mype*self.nz
                c0 += offset; c1 += offset

            data = self._load(a0, a1, b0, b1, c0, c1, idump, dir_, nghost, 0)

            self.rho[i0:i1,j0:j1,k0:k1]    = data['rho']
            self.E[i0:i1,j0:j1,k0:k1]      = data['E']
            self.rhou[i0:i1,j0:j1,k0:k1,0] = data['rhou'][:,:,:,0]
            self.rhou[i0:i1,j0:j1,k0:k1,1] = data['rhou'][:,:,:,1]
            self.rhou[i0:i1,j0:j1,k0:k1,2] = data['rhou'][:,:,:,2]
            self.B[i0:i1,j0:j1,k0:k1,0]    = data['B'][:,:,:,0]
            self.B[i0:i1,j0:j1,k0:k1,1]    = data['B'][:,:,:,1]
            self.B[i0:i1,j0:j1,k0:k1,2]    = data['B'][:,:,:,2]
            self.Br[i0:i1,j0:j1,k0:k1,0]   = data['Br'][:,:,:,0]
            self.Br[i0:i1,j0:j1,k0:k1,1]   = data['Br'][:,:,:,1]
            self.Br[i0:i1,j0:j1,k0:k1,2]   = data['Br'][:,:,:,2]

## Parallel HDF5 data class
@DumsesDataFactory.registerData('h5')
class DataPHDF5Reader(DataParallelReader):
    def __init__(self, nx, ny, nz, nxslice, nyslice, nzslice, nxglob, nyglob \
                     , nzglob, nproc):
        DataParallelReader.__init__(self, nx, ny, nz, nxslice, nyslice, nzslice\
                                , nxglob, nyglob, nzglob, nproc)
    
    ## _load method; reads a chunck of data from the parallel output file
    def _load(self, i0, i1, j0, j1, k0, k1, idump=1, dir_='./', nghost=3 \
                  , mype=0, ext='h5'):
        filename = self._getFile(idump, dir_, mype, ext)
        
        f = h5.File(filename, 'r')
        rho = f['rho'][k0:k1,j0:j1,i0:i1].transpose()
        vx  = f['vx'][k0:k1,j0:j1,i0:i1].transpose()
        vy  = f['vy'][k0:k1,j0:j1,i0:i1].transpose()
        vz  = f['vz'][k0:k1,j0:j1,i0:i1].transpose()
        E   = f['E'][k0:k1,j0:j1,i0:i1].transpose()
        Bx  = f['Bx'][k0:k1,j0:j1,i0:i1].transpose()
        By  = f['By'][k0:k1,j0:j1,i0:i1].transpose()
        Bz  = f['Bz'][k0:k1,j0:j1,i0:i1].transpose()
        Bxr = f['Bxr'][k0:k1,j0:j1,i0:i1].transpose()
        Byr = f['Byr'][k0:k1,j0:j1,i0:i1].transpose()
        Bzr = f['Bzr'][k0:k1,j0:j1,i0:i1].transpose()
        f.close()
        return dict({'rho': rho, 'E': E \
                  , 'rhou': np.array([vx,vy,vz]).transpose(1,2,3,0) \
                  , 'B': np.array([Bx,By,Bz]).transpose(1,2,3,0) \
                  , 'Br': np.array([Bxr,Byr,Bzr]).transpose(1,2,3,0)})

## Parallel NetCDF data class
@DumsesDataFactory.registerData('nc')
class DataNetCDFReader(DataParallelReader):
    def __init__(self, nx, ny, nz, nxslice, nyslice, nzslice, nxglob, nyglob \
                     , nzglob, nproc):
        DataParallelReader.__init__(self, nx, ny, nz, nxslice, nyslice, nzslice\
                                , nxglob, nyglob, nzglob, nproc)
    
    ## _load method; reads a chunck of data from the parallel output file
    def _load(self, i0, i1, j0, j1, k0, k1, idump=1, dir_='./', nghost=3 \
                  , mype=0, ext='nc'):
        filename = self._getFile(idump, dir_, mype, ext)
        
        f = nc.Dataset(filename, 'r')
        rho = f.variables['rho'][k0:k1,j0:j1,i0:i1].transpose()
        vx  = f.variables['vx'][k0:k1,j0:j1,i0:i1].transpose()
        vy  = f.variables['vy'][k0:k1,j0:j1,i0:i1].transpose()
        vz  = f.variables['vz'][k0:k1,j0:j1,i0:i1].transpose()
        E   = f.variables['E'][k0:k1,j0:j1,i0:i1].transpose()
        Bx  = f.variables['Bx'][k0:k1,j0:j1,i0:i1].transpose()
        By  = f.variables['By'][k0:k1,j0:j1,i0:i1].transpose()
        Bz  = f.variables['Bz'][k0:k1,j0:j1,i0:i1].transpose()
        Bxr = f.variables['Bxr'][k0:k1,j0:j1,i0:i1].transpose()
        Byr = f.variables['Byr'][k0:k1,j0:j1,i0:i1].transpose()
        Bzr = f.variables['Bzr'][k0:k1,j0:j1,i0:i1].transpose()
        f.close()
        return dict({'rho': rho, 'E': E \
                  , 'rhou': np.array([vx,vy,vz]).transpose(1,2,3,0) \
                  , 'B': np.array([Bx,By,Bz]).transpose(1,2,3,0) \
                  , 'Br': np.array([Bxr,Byr,Bzr]).transpose(1,2,3,0)})

## Old parallel HDF5 data class
@DumsesDataFactory.registerData('oldphdf5')
class DataOldPHDF5Reader(DataParallelReader):
    def __init__(self, nx, ny, nz, nxslice, nyslice, nzslice, nxglob, nyglob \
                     , nzglob, nproc):
        DataParallelReader.__init__(self, nx, ny, nz, nxslice, nyslice \
                                  , nzslice, nxglob, nyglob, nzglob, nproc)
    
    ## _inline method to check the contiguity of data
    def _inline(self, idump=1, dir_='./', mype=0, nghost=3):
        filename = self._getFile(idump, dir_, mype, ext='')
        f    = h5.File(filename, 'r')
        rho  = f['rho']
        zdim = rho.shape[0]
        f.close()
        zdimAct = ((self.nz + 2*nghost)*self.nzslice if self.nz > 1 \
               else self.nzslice)
        return zdim == zdimAct

    ## _load method; reads a chunck of data from the parallel output file
    def _load(self, i0, i1, j0, j1, k0, k1, idump=1, dir_='./', nghost=3 \
                  , mype=0, ext=''):
        filename = self._getFile(idump, dir_, mype, ext)
        
        f = h5.File(filename, 'r')
        rho = f['rho'][k0:k1,j0:j1,i0:i1].transpose()
        vx  = f['rho_vx'][k0:k1,j0:j1,i0:i1].transpose()
        vy  = f['rho_vy'][k0:k1,j0:j1,i0:i1].transpose()
        vz  = f['rho_vz'][k0:k1,j0:j1,i0:i1].transpose()
        E   = f['E'][k0:k1,j0:j1,i0:i1].transpose()
        Bx  = f['Bx'][k0:k1,j0:j1,i0:i1].transpose()
        By  = f['By'][k0:k1,j0:j1,i0:i1].transpose()
        Bz  = f['Bz'][k0:k1,j0:j1,i0:i1].transpose()
        Bxr = f['Bxr'][k0:k1,j0:j1,i0:i1].transpose()
        Byr = f['Byr'][k0:k1,j0:j1,i0:i1].transpose()
        Bzr = f['Bzr'][k0:k1,j0:j1,i0:i1].transpose()
        f.close()
        return dict({'rho': rho, 'E': E \
                  , 'rhou': np.array([vx,vy,vz]).transpose(1,2,3,0) \
                  , 'B': np.array([Bx,By,Bz]).transpose(1,2,3,0) \
                  , 'Br': np.array([Bxr,Byr,Bzr]).transpose(1,2,3,0)})

## Old parallel NetCDF data class
@DumsesDataFactory.registerData('oldpnetcdf')
class DataOldNetCDFReader(DataParallelReader):
    def __init__(self, nx, ny, nz, nxslice, nyslice, nzslice, nxglob, nyglob \
                     , nzglob, nproc):
        DataParallelReader.__init__(self, nx, ny, nz, nxslice, nyslice \
                                  , nzslice, nxglob, nyglob, nzglob, nproc)
    
    ## _inline method to check the contiguity of data
    def _inline(self, idump=1, dir_='./', mype=0, nghost=3):
        filename = self._getFile(idump, dir_, mype, ext='')
        f    = nc.Dataset(filename, 'r')
        rho  = f.variables['rho']
        zdim = rho.shape[0]
        f.close()
        zdimAct = ((self.nz + 2*nghost)*self.nzslice if self.nz > 1 \
               else self.nzslice)
        return zdim == zdimAct

    ## _load method; reads a chunck of data from the parallel output file
    def _load(self, i0, i1, j0, j1, k0, k1, idump=1, dir_='./', nghost=3 \
                  , mype=0, ext=''):
        filename = self._getFile(idump, dir_, mype, ext)
        
        f = nc.Dataset(filename, 'r')
        rho = f.variables['rho'][k0:k1,j0:j1,i0:i1].transpose()
        vx  = f.variables['rho_vx'][k0:k1,j0:j1,i0:i1].transpose()
        vy  = f.variables['rho_vy'][k0:k1,j0:j1,i0:i1].transpose()
        vz  = f.variables['rho_vz'][k0:k1,j0:j1,i0:i1].transpose()
        E   = f.variables['E'][k0:k1,j0:j1,i0:i1].transpose()
        Bx  = f.variables['Bx'][k0:k1,j0:j1,i0:i1].transpose()
        By  = f.variables['By'][k0:k1,j0:j1,i0:i1].transpose()
        Bz  = f.variables['Bz'][k0:k1,j0:j1,i0:i1].transpose()
        Bxr = f.variables['Bxr'][k0:k1,j0:j1,i0:i1].transpose()
        Byr = f.variables['Byr'][k0:k1,j0:j1,i0:i1].transpose()
        Bzr = f.variables['Bzr'][k0:k1,j0:j1,i0:i1].transpose()
        f.close()
        return dict({'rho': rho, 'E': E \
                  , 'rhou': np.array([vx,vy,vz]).transpose(1,2,3,0) \
                  , 'B': np.array([Bx,By,Bz]).transpose(1,2,3,0) \
                  , 'Br': np.array([Bxr,Byr,Bzr]).transpose(1,2,3,0)})

## Function to read a binary array from a Fortran binary file
def getArray(fid, nbCount, dtype):
    bitsType = 'i4'
    padBegin = np.fromfile(fid, count=1, dtype=bitsType)
    array    = np.fromfile(fid, count=nbCount, dtype=dtype)
    padEnd   = np.fromfile(fid, count=1, dtype=bitsType)
    try:
        assert padBegin == padEnd
    except AssertionError:
        print('Something went wrong reading binary data!')
    return array

## Function to return a dictionary out of a DumsesData object
def dictionarize(data):
    try:
        assert(isinstance(data, DumsesData))
    except AssertionError:
        print("Error: 'data' must be a DumsesData object")
        sys.exit(2)
    dict_ = {}
    for key in data._data:
        if len(data.__getattribute__(key).shape) == 4:
            dict_[key + 'x'] = data.__getattribute__(key)[:,:,:,0]
            dict_[key + 'y'] = data.__getattribute__(key)[:,:,:,1]
            dict_[key + 'z'] = data.__getattribute__(key)[:,:,:,2]
        else:
            dict_[key] = data.__getattribute__(key)
    return dict_
