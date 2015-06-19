#===============================================================================
## \file formatChecker
# \brief
# \b DUMSES-Hybrid:
# These is a format checker for DUMSES outputs
# \author
# Marc Joos <marc.joos@cea.fr>, Sebastien Fromang
# \copyright
# Copyrights 2013-2015, CEA.
# This file is distributed under the CeCILL-A & GNU/GPL licenses, see
# <http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html> and
# <http://www.gnu.org/licenses/>
# \date
# \b created:          03-27-2015
# \b last \b modified: 06-15-2015

#===============================================================================
import re
import sys
import glob
import numpy as np
import h5py  as h5
import netCDF4 as nc

## Function to check the output format
def formatChecker(dir_='./', ndump=None, verbose=False):
    """Check the DUMSES output format.

    Usage:
      old, ext = formatChecker([dir_, ndump, verbose])
    With:
      dir_ (optional):    directory containing the outputs. Default: ./
      ndump (optional):   index of the output to read. Default: first output found
      verbose (optional): verbose mode
    Returns:
      old: True if the simulation was made with dumses_mpi
      ext: file format
    """
    
    ext = ''
    old = False
    format_ = ''

    loutput = np.sort(glob.glob(dir_ + '/output_*'))
    try:
        assert loutput
    except AssertionError:
        print('The directory provided contains no output.')
        sys.exit(0)
    except ValueError:
        pass
    if ndump:
        output = dir_ + '/output_%06d' %ndump
    else:
        output = loutput[0]

    tmpLfiles = np.sort(glob.glob(output + '/*'))
    outputFile = re.compile("(.*)(data|slices)\.(\d{6})(\.\w+)?")
    lfiles = []
    for file_ in tmpLfiles:
        matched = outputFile.match(file_)
        if matched:
            lfiles.append(file_)
            ismatched = matched
    nfiles = len(lfiles)

    if ismatched:
        fname = ismatched.group(2)

        if fname == "data":
            extension = ismatched.group(4)
            if extension == ".bin":
                format_ = "sequential binary"
                ext = "bin"
            elif extension == ".nc":
                format_ = "parallel NetCDF"
                ext = "nc"
            elif extension == ".h5":
                if nfiles > 1:
                    format_ = "sequential HDF5"
                else:
                    format_ = "(parallel) HDF5"
                ext = "h5"
            else:
                raise IOError('File format not recognized')

        elif fname == "slices":
            old = True
            filename = output + '/' + lfiles[0].split('/')[-1]
            try:
                with h5.File(filename, 'r') as f:
                    rho = f['uin']
            except IOError:
                pass
            except KeyError:
                pass
            else:
                format_ = 'sequential HDF5'
                ext = 'hdf5' 
            
            try:
                with h5.File(filename, 'r') as f:
                    rho = f['rho']
            except IOError:
                pass
            except KeyError:
                pass
            else:
                format_ = 'parallel HDF5'
                ext = 'phdf5'
            
            try:
                major = int(nc.__version__.split('.')[0])
                if major >= 1:
                    with nc.Dataset(filename) as f:
                        rho = f.variables['rho']
                else:
                    f   = nc.Dataset(filename, 'r')
                    rho = f.variables['rho']
                    f.close()
            except IOError:
                pass
            except KeyError:
                pass
            except RuntimeError:
                pass
            except AttributeError:
                pass
            else:
                if not format_: 
                    format_ = 'parallel NetCDF'
                    ext = 'pnetcdf'
            
            if not format_: 
                format_ = 'binary'
                ext = 'binary'
            ext = 'old' + ext

        else:
            raise IOError('File format not recognized')

    if verbose: 
        print("This is a " + format_ \
           + (" dumses_mpi" if old else " dumses_hybrid") + " file.")
    return old, ext
