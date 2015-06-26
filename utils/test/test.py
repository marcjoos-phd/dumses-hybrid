#!/usr/bin/python
#===============================================================================
## \file test.py
# \brief
# \b DUMSES-Hybrid:
# This is test suite for DUMSES-hybrid
# \author
# Marc Joos <marc.joos@cea.fr>
# \copyright
# Copyrights 2013-2015, CEA.
# This file is distributed under the CeCILL-A & GNU/GPL licenses, see
# <http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html> and
# <http://www.gnu.org/licenses/>
# \date
# \b created:          12-20-2013 
# \b last \b modified: 06-26-2015

#===============================================================================
import os
import sys
import glob
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import pylab as pl
if sys.version_info.major <= 2:
    from commands import getstatusoutput as cmd
else:
    from subprocess import getstatusoutput as cmd

bold    = "\033[1m"
reset   = "\033[0;0m"

## Compare current configuration with the last execution of configure script
def compConf(fname='conf.log', mpi='', hdf='/', cdf='/', fcc='gfortran' \
           , oacc=0, pro='shock_tube', ndim=1, iso=0):
    try:
        with open(fname) as logfile:
            st, out = cmd("ls Makefile")
            if st != 0:
                return 0
            else:
                logline = logfile.read()
                args = logline.split(" --")
                args = args[1:]
                larg = ['with-mpi', 'with-phdf5', 'with-pnetcdf' \
                , 'with-fortran-compiler', 'with-openacc', 'problem' \
                , 'ndim', 'isothermal']
                linp = ['' for i in xrange(len(larg))]
                linp[larg.index('with-openacc')] = 0
                for i in xrange(len(args)):
                    if args[i] != "with-openacc":
                        targ, tinp = args[i].split('=')
                    else:
                        targ = "with-openacc"; tinp=1
                    linp[larg.index(targ)] = tinp
                linp[-1] = linp[-1].split('\n')[0]
                lnew = np.array([mpi, hdf, cdf, fcc, oacc, pro, ndim, iso])
                if np.array(lnew == linp).all():
                    return 1
                else:
                    return 0
    except IOError:
        print("File 'conf.log' does not exist")
        return 0

## Run configure script
def runConfigure(mpi='/', hdf='/', cdf='/', fcc='gfortran', oacc=0 \
                     , pro='shock_tube', ndim=1, iso=0):
    st, out = cmd('./configure' + ('' if mpi == '' else ' --with-mpi=' + mpi) \
                + ('' if hdf == '' else ' --with-phdf5=' + hdf) \
                + ('' if cdf == '' else ' --with-pnetcdf=' + cdf) \
                + ('' if fcc == '' else ' --with-fortran-compiler=' + fcc) \
                + ('' if oacc == 0 else ' --with-openacc') \
                + ('' if pro == '' else ' --problem=' + pro) \
                + ('' if ndim == '' else ' --ndim=%d' %ndim) \
                + ('' if iso == '' else ' --iso=%d' %iso))
    return st, out

## Create input for the given problem
def createInput(tlim="1.", verbose=".false.", debug=".false." \
              , bdtypex="'zero_gradient'", bdtypey="'zero_gradient'" \
              , bdtypez="'zero_gradient'", riemann="'hlld'", riemann2d="'hlld'"\
              , slope_type="2", courant="0.7", nx="32", ny="32", nz="32" \
              , xmin="-1.d0", xmax="1.d0", ymin="-1.d0", ymax="1.d0" \
              , zmin="0.d0", zmax="1.d0", Omega0="0.d0", ciso="0.d0" \
              , gamma="1.001d0", dtdump=".1", dthist="-1", io_type="'binary'" \
              , nxslice="1", nyslice="1", nzslice="1"):
    st, out = cmd('cp ../input.template input')
    template = open('input', 'rt').read()
    data = {"tlim": tlim,
            "verbose": verbose,
            "debug": debug,
            "bdtypex": bdtypex,
            "bdtypey": bdtypey,
            "bdtypez": bdtypez,
            "riemann": riemann,
            "riemann2d": riemann2d,
            "slope_type": slope_type,
            "courant": courant,
            "nx": nx,
            "ny": ny,
            "nz": nz,
            "xmin": xmin,
            "xmax": xmax,
            "ymin": ymin,
            "ymax": ymax,
            "zmin": zmin,
            "zmax": zmax,
            "Omega0": Omega0,
            "ciso": ciso,
            "gamma": gamma,
            "dtdump": dtdump,
            "dthist": dthist,
            "io_type": io_type,
            "nxslice": nxslice,
            "nyslice": nyslice,
            "nzslice": nzslice
            }

    with open('input', 'wt') as output:
        output.write(template %data)

## Launch the problem
def launchProblem(currentDir, dumsesRoot, mpi='/', hdf='/', cdf='/' \
              , fcc='gfortran', oacc=0, pro='shock_tube', ndim=1, iso=0 \
              , tlim="1.", verbose=".false.", debug=".false." \
              , bdtypex="'zero_gradient'", bdtypey="'zero_gradient'" \
              , bdtypez="'zero_gradient'", riemann="'hlld'", riemann2d="'hlld'"\
              , slope_type="2", courant="0.7", nx="32", ny="32", nz="32" \
              , xmin="-1.d0", xmax="1.d0", ymin="-1.d0", ymax="1.d0" \
              , zmin="0.d0", zmax="1.d0", Omega0="0.d0", ciso="0.d0" \
              , gamma="1.001d0", dtdump=".1", dthist="-1", io_type="'binary'" \
              , nxslice="1", nyslice="1", nzslice="1", direction='x' \
              , beta="400.d0", type_="Radial", amp="0.d0"):
    os.chdir(dumsesRoot)
    same = compConf(fname='conf.log', mpi=mpi, hdf=hdf, cdf=cdf, fcc=fcc \
                        , oacc=oacc, pro=pro, ndim=ndim, iso=iso)
    if not(same):
        st, out = runConfigure(mpi=mpi, hdf=hdf, cdf=cdf, fcc=fcc, oacc=oacc \
                               , pro=pro, ndim=ndim, iso=iso)

        print(out + '\n')
        st, out = cmd('make clean')
        print(out + '\n')
        st, out = cmd('./make.py')
        try:
            assert(st == 0)
        except AssertionError:
            print("Build failed!")
            sys.exit(2)
        print(out + '\n')
        st, out = cmd('cp bin/dumses ' + currentDir + '/tmp/.')
    else:
        st, out = cmd('make')
        try:
            assert(st == 0)
        except AssertionError:
            print("Build failed!")
            sys.exit(2)
        print(out + '\n')
        st, out = cmd('cp bin/dumses ' + currentDir + '/tmp/.')
    
    os.chdir(currentDir + '/tmp')

    createInput(tlim=tlim, verbose=verbose, debug=debug \
              , bdtypex=bdtypex, bdtypey=bdtypey \
              , bdtypez=bdtypez, riemann=riemann, riemann2d=riemann2d \
              , slope_type=slope_type, courant=courant, nx=nx, ny=ny, nz=nz \
              , xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax \
              , zmin=zmin, zmax=zmax, Omega0=Omega0, ciso=ciso \
              , gamma=gamma, dtdump=dtdump, dthist=dthist, io_type=io_type \
              , nxslice=nxslice, nyslice=nyslice, nzslice=nzslice)

    lproDirect = ['shock_tube', 'orszag_tang', 'magnetic_loop', 'wind_tunnel']
    if pro in lproDirect:
        f = open('input', 'a')
        initparam="\n&init_params\n        direction = '%s'\n        /\n" \
            %direction
        f.write(initparam)
        f.close()
    if pro == 'mri':
        f = open('input', 'a')
        initparam="\n&init_params\n        beta = %s\n        type = '%s'\n        amp = %s\n        /\n" \
            %(beta, type_, amp)
        f.write(initparam)
        f.close()

    if mpi == '/' or mpi == '':
        st, out = cmd('./dumses')
    else:
        st, out = cmd(mpi + '/bin/mpirun -np %d ./dumses > log' \
                          %(int(nxslice)*int(nyslice)*int(nzslice)))
    print(out + '\n')
    
    os.chdir(currentDir)

## Plot result of the given test
def plotResult(dumsesRoot, idump=10, ndim=1, fdir='./tmp/', io_type='binary' \
             , pro="", save=False, datdim=None, direction='x', form=['png'] \
             , vmind=None, vmaxd=None, vminE=None, vmaxE=None \
             , vminv=None, vmaxv=None, vminB=None, vmaxB=None, diff=False):
    sys.path.append(dumsesRoot + 'utils/')
    sys.path = sys.path[::-1]
    import dumpy as dp

    if isinstance(form, list):
        listformat = ['.' + f.split('.')[-1] for f in form]
    elif isinstance(form, str):
        if form == 'all':
            listformat = ['.png', '.pdf', '.eps']
        else:
            listformat = ['.' + form.split('.')[-1]]
    
    if diff: 
        dirname = pro.lower()
        dirname = dirname.replace(' ', '_')
        dirname = 'fig/reference/' + dirname + '/' + direction + '/'
        datRef  = dp.DumsesData()
        datRef.load(idump, filedir=dirname, io_type='binary', ndim=datdim)

    dat = dp.DumsesData()
    dat.load(idump, dir_=fdir, ndim=datdim)

    if ndim == 1:
        if diff:
            pl.figure(figsize=(12,8))
        else:
            pl.figure(figsize=(8,4))
    elif ndim == 2:
        if diff:
            pl.figure(figsize=(10,8))
        else:
            pl.figure(figsize=(8,8))
    pl.subplots_adjust(bottom=0.15)
    if diff:
        pl.suptitle(pro + ', $' + direction + '$-direction' \
                        , fontsize=12, x=.5, y=.99)
    else:
        pl.suptitle(pro, fontsize=12, x=.5, y=.99)
    if ndim == 1:
        if diff:
            pl.subplot(331)
            pl.title("Current run")
            pl.ylabel(r'$\rho$')
            if direction == 'x':
                pl.plot(dat.rho[:,0,0])
            elif direction == 'y':
                pl.plot(dat.rho[0,:,0])
            else:
                pl.plot(dat.rho[0,0,:])
            pl.ylim(ymin=vmind, ymax=vmaxd)
            pl.subplot(334)
            pl.ylabel(r'$E$')
            if direction == 'x':
                pl.plot(dat.E[:,0,0])
            elif direction == 'y':
                pl.plot(dat.E[0,:,0])
            else:
                pl.plot(dat.E[0,0,:])
            pl.ylim(ymin=vminE, ymax=vmaxE)
            pl.subplot(337)
            if direction == 'x':
                pl.ylabel(r'$v_x$')
                pl.xlabel('x')
                pl.plot(dat.rhou[:,0,0,0]/dat.rho[:,0,0])
            elif direction == 'y':
                pl.ylabel(r'$v_y$')
                pl.xlabel('y')
                pl.plot(dat.rhou[0,:,0,1]/dat.rho[0,:,0])
            else:
                pl.ylabel(r'$v_z$')
                pl.xlabel('z')
                pl.plot(dat.rhou[0,0,:,2]/dat.rho[0,0,:])
            pl.ylim(ymin=vminv, ymax=vmaxv)

            pl.subplot(332)
            pl.title("Reference run")
            if direction == 'x':
                pl.plot(datRef.rho[:,0,0])
            elif direction == 'y':
                pl.plot(datRef.rho[0,:,0])
            else:
                pl.plot(datRef.rho[0,0,:])
            pl.ylim(ymin=vmind, ymax=vmaxd)
            pl.subplot(335)
            if direction == 'x':
                pl.plot(datRef.E[:,0,0])
            elif direction == 'y':
                pl.plot(datRef.E[0,:,0])
            else:
                pl.plot(datRef.E[0,0,:])
            pl.ylim(ymin=vminE, ymax=vmaxE)
            pl.subplot(338)
            if direction == 'x':
                pl.xlabel('x')
                pl.plot(datRef.rhou[:,0,0,0]/datRef.rho[:,0,0])
            elif direction == 'y':
                pl.xlabel('y')
                pl.plot(datRef.rhou[0,:,0,1]/datRef.rho[0,:,0])
            else:
                pl.xlabel('z')
                pl.plot(datRef.rhou[0,0,:,2]/datRef.rho[0,0,:])
            pl.ylim(ymin=vminv, ymax=vmaxv)

            pl.subplot(333)
            pl.title("Error")
            if direction == 'x':
                pl.plot(abs((dat.rho[:,0,0] - datRef.rho[:,0,0])/datRef.rho[:,0,0]))
            elif direction == 'y':
                pl.plot(abs((dat.rho[0,:,0] - datRef.rho[0,:,0])/datRef.rho[0,:,0]))
            else:
                pl.plot(abs((dat.rho[0,0,:] - datRef.rho[0,0,:])/datRef.rho[0,0,:]))
            pl.ylim(ymin=vmind, ymax=vmaxd)
            pl.subplot(336)
            if direction == 'x':
                pl.plot(abs((dat.E[:,0,0] - datRef.E[:,0,0])/datRef.E[:,0,0]))
            elif direction == 'y':
                pl.plot(abs((dat.E[0,:,0] - datRef.E[0,:,0])/datRef.E[0,:,0]))
            else:
                pl.plot(abs((dat.E[0,0,:] - datRef.E[0,0,:])/datRef.E[0,0,:]))
            pl.ylim(ymin=vminE, ymax=vmaxE)
            pl.subplot(339)
            if direction == 'x':
                pl.xlabel('x')
                pl.plot(abs((dat.rhou[:,0,0,1] - datRef.rhou[:,0,0,0])/datRef.rhou[:,0,0,0]))
            elif direction == 'y':
                pl.xlabel('y')
                pl.plot(abs((dat.rhou[0,:,0,1] - datRef.rhou[0,:,0,1])/datRef.rhou[0,:,0,1]))
            else:
                pl.xlabel('z')
                pl.plot(abs((dat.rhou[0,0,:,2] - datRef.rhou[0,0,:,2])/datRef.rhou[0,0,:,2]))
            pl.ylim(ymin=vminv, ymax=vmaxv)
        else:
            pl.subplot(131)
            pl.title(r'$\rho$')
            if direction == 'x':
                pl.xlabel('x')
                pl.plot(dat.rho[:,0,0])
                pl.ylim(ymin=vmind, ymax=vmaxd)
            elif direction == 'y':
                pl.xlabel('y')
                pl.plot(dat.rho[0,:,0])
                pl.ylim(ymin=vmind, ymax=vmaxd)
            else:
                pl.xlabel('z')
                pl.plot(dat.rho[0,0,:])
                pl.ylim(ymin=vmind, ymax=vmaxd)
            pl.subplot(132)
            pl.title(r'$E$')
            if direction == 'x':
                pl.xlabel('x')
                pl.plot(dat.E[:,0,0])
                pl.ylim(ymin=vminE, ymax=vmaxE)
            elif direction == 'y':
                pl.xlabel('y')
                pl.plot(dat.E[0,:,0])
                pl.ylim(ymin=vminE, ymax=vmaxE)
            else:
                pl.xlabel('z')
                pl.plot(dat.E[0,0,:])
                pl.ylim(ymin=vminE, ymax=vmaxE)
            pl.subplot(133)
            if direction == 'x':
                pl.title(r'$v_x$')
                pl.xlabel('x')
                pl.plot(dat.rhou[:,0,0,0]/dat.rho[:,0,0])
                pl.ylim(ymin=vminv, ymax=vmaxv)
            elif direction == 'y':
                pl.title(r'$v_y$')
                pl.xlabel('y')
                pl.plot(dat.rhou[0,:,0,1]/dat.rho[0,:,0])
                pl.ylim(ymin=vminv, ymax=vmaxv)
            else:
                pl.title(r'$v_z$')
                pl.xlabel('z')
                pl.plot(dat.rhou[0,0,:,2]/dat.rho[0,0,:])
                pl.ylim(ymin=vminv, ymax=vmaxv)
    if ndim == 2:
        if diff:
            pl.subplot(4,3,1)
            pl.title("Current run")
            pl.ylabel(r'$\rho$')
            if direction == 'x':
                pl.imshow(dat.rho[:,:,0], vmin=vmind, vmax=vmaxd)
            elif direction == 'y':
                pl.imshow(dat.rho[0,:,:], vmin=vmind, vmax=vmaxd)
            else:
                pl.imshow(dat.rho[:,0,:], vmin=vmind, vmax=vmaxd)
            pl.colorbar()
            pl.subplot(4,3,4)
            pl.ylabel(r'$E$')
            if direction == 'x':
                pl.imshow(dat.E[:,:,0], vmin=vminE, vmax=vmaxE)
            elif direction == 'y':
                pl.imshow(dat.E[0,:,:], vmin=vminE, vmax=vmaxE)
            else:
                pl.imshow(dat.E[:,0,:], vmin=vminE, vmax=vmaxE)
            pl.colorbar()
            pl.subplot(4,3,7)
            if direction == 'x':
                pl.ylabel(r'$v_x$')
                pl.imshow(dat.rhou[:,:,0,0]/dat.rho[:,:,0], vmin=vminv, vmax=vmaxv)
            elif direction == 'y':
                pl.ylabel(r'$v_y$')
                pl.imshow(dat.rhou[0,:,:,1]/dat.rho[0,:,:], vmin=vminv, vmax=vmaxv)
            else:
                pl.ylabel(r'$v_z$')
                pl.imshow(dat.rhou[:,0,:,2]/dat.rho[:,0,:], vmin=vminv, vmax=vmaxv)
            pl.colorbar()
            pl.subplot(4,3,10)
            if direction == 'x':
                pl.ylabel(r'$B_x$')
                pl.imshow(dat.B[:,:,0,0], vmin=vminB, vmax=vmaxB)
            elif direction == 'y':
                pl.ylabel(r'$B_y$')
                pl.imshow(dat.B[0,:,:,1], vmin=vminB, vmax=vmaxB)
            else:
                pl.ylabel(r'$B_z$')
                pl.imshow(dat.B[:,0,:,2], vmin=vminB, vmax=vmaxB)
            pl.colorbar()

            pl.subplot(4,3,2)
            pl.title("Reference run")
            if direction == 'x':
                pl.imshow(datRef.rho[:,:,0], vmin=vmind, vmax=vmaxd)
            elif direction == 'y':
                pl.imshow(datRef.rho[0,:,:], vmin=vmind, vmax=vmaxd)
            else:
                pl.imshow(datRef.rho[:,0,:], vmin=vmind, vmax=vmaxd)
            pl.colorbar()
            pl.subplot(4,3,5)
            if direction == 'x':
                pl.imshow(datRef.E[:,:,0], vmin=vminE, vmax=vmaxE)
            elif direction == 'y':
                pl.imshow(datRef.E[0,:,:], vmin=vminE, vmax=vmaxE)
            else:
                pl.imshow(datRef.E[:,0,:], vmin=vminE, vmax=vmaxE)
            pl.colorbar()
            pl.subplot(4,3,8)
            if direction == 'x':
                pl.imshow(datRef.rhou[:,:,0,0]/datRef.rho[:,:,0], vmin=vminv, vmax=vmaxv)
            elif direction == 'y':
                pl.imshow(datRef.rhou[0,:,:,1]/datRef.rho[0,:,:], vmin=vminv, vmax=vmaxv)
            else:
                pl.imshow(datRef.rhou[:,0,:,2]/datRef.rho[:,0,:], vmin=vminv, vmax=vmaxv)
            pl.colorbar()
            pl.subplot(4,3,11)
            if direction == 'x':
                pl.imshow(datRef.B[:,:,0,0], vmin=vminB, vmax=vmaxB)
            elif direction == 'y':        
                pl.imshow(datRef.B[0,:,:,1], vmin=vminB, vmax=vmaxB)
            else:                         
                pl.imshow(datRef.B[:,0,:,2], vmin=vminB, vmax=vmaxB)
            pl.colorbar()

            pl.subplot(4,3,3)
            pl.title("Error")
            if direction == 'x':
                pl.imshow(abs((dat.rho[:,:,0] - datRef.rho[:,:,0])/datRef.rho[:,:,0]))
            elif direction == 'y':
                pl.imshow(abs((dat.rho[0,:,:] - datRef.rho[0,:,:])/datRef.rho[0,:,:]))
            else:
                pl.imshow(abs((dat.rho[:,0,:] - datRef.rho[:,0,:])/datRef.rho[:,0,:]))
            pl.colorbar()
            pl.subplot(4,3,6)
            if direction == 'x':
                pl.imshow(abs((dat.E[:,:,0] - datRef.E[:,:,0])/datRef.E[:,:,0]))
            elif direction == 'y':
                pl.imshow(abs((dat.E[0,:,:] - datRef.E[0,:,:])/datRef.E[0,:,:]))
            else:
                pl.imshow(abs((dat.E[:,0,:] - datRef.E[:,0,:])/datRef.E[:,0,:]))
            pl.colorbar()
            pl.subplot(4,3,9)
            if direction == 'x':
                pl.imshow(abs((dat.rhou[:,:,0,0] - datRef.rhou[:,:,0,0])/datRef.rhou[:,:,0,0]))
            elif direction == 'y':                                    
                pl.imshow(abs((dat.rhou[0,:,:,1] - datRef.rhou[0,:,:,1])/datRef.rhou[0,:,:,1]))
            else:                                                     
                pl.imshow(abs((dat.rhou[:,0,:,2] - datRef.rhou[:,0,:,2])/datRef.rhou[:,0,:,2]))
            pl.colorbar()
            pl.subplot(4,3,12)
            if direction == 'x':
                pl.imshow(abs((dat.B[:,:,0,0] - datRef.B[:,:,0,0])/datRef.B[:,:,0,0]))
            elif direction == 'y':                              
                pl.imshow(abs((dat.B[0,:,:,1] - datRef.B[0,:,:,1])/datRef.B[0,:,:,1]))
            else:                                               
                pl.imshow(abs((dat.B[:,0,:,2] - datRef.B[:,0,:,2])/datRef.B[:,0,:,2]))
            pl.colorbar()
        else:

            pl.subplot(221)
            pl.title(r'$\rho$')
            if direction == 'x':
                pl.ylabel('y')
                pl.imshow(dat.rho[:,:,0], vmin=vmind, vmax=vmaxd)
            elif direction == 'y':
                pl.ylabel('z')
                pl.imshow(dat.rho[0,:,:], vmin=vmind, vmax=vmaxd)
            else:
                pl.ylabel('x')
                pl.imshow(dat.rho[:,0,:], vmin=vmind, vmax=vmaxd)
            pl.colorbar()
            pl.subplot(222)
            pl.title(r'$E$')
            if direction == 'x':
                pl.imshow(dat.E[:,:,0], vmin=vminE, vmax=vmaxE)
            elif direction == 'y':
                pl.imshow(dat.E[0,:,:], vmin=vminE, vmax=vmaxE)
            else:
                pl.imshow(dat.E[:,0,:], vmin=vminE, vmax=vmaxE)
            pl.colorbar()
            pl.subplot(223)
            if direction == 'x':
                pl.title(r'$v_x$')
                pl.xlabel('x')
                pl.ylabel('y')
                pl.imshow(dat.rhou[:,:,0,0]/dat.rho[:,:,0], vmin=vminv, vmax=vmaxv)
            elif direction == 'y':
                pl.title(r'$v_y$')
                pl.xlabel('y')
                pl.ylabel('z')
                pl.imshow(dat.rhou[0,:,:,1]/dat.rho[0,:,:], vmin=vminv, vmax=vmaxv)
            else:
                pl.title(r'$v_z$')
                pl.xlabel('z')
                pl.ylabel('x')
                pl.imshow(dat.rhou[:,0,:,2]/dat.rho[:,0,:], vmin=vminv, vmax=vmaxv)
            pl.colorbar()
            pl.subplot(224)
            if direction == 'x':
                pl.title(r'$B_x$')
                pl.xlabel('x')
                pl.imshow(dat.B[:,:,0,0], vmin=vminB, vmax=vmaxB)
            elif direction == 'y':
                pl.title(r'$B_y$')
                pl.xlabel('y')
                pl.imshow(dat.B[0,:,:,1], vmin=vminB, vmax=vmaxB)
            else:
                pl.title(r'$B_z$')
                pl.xlabel('z')
                pl.imshow(dat.B[:,0,:,2], vmin=vminB, vmax=vmaxB)
            pl.colorbar()
    if save:
        lname    = pro.split()
        namefile = lname[0]
        if len(lname) > 1:
            namefile = [namefile + '_' + lname[i] for i in xrange(1, len(lname))][0]
        if ndim == 1 or ndim == 2:
            if diff:
                namefile = namefile + '_Err_' + direction + '_%02d' %idump
            else:
                namefile = namefile + '_' + direction + '_%02d' %idump
        for form in listformat:
            print("Save file: fig/" + namefile + form)
            pl.savefig('fig/' + namefile + form)
        pl.close()
    else:
        pl.show()

## Plot all the results from the given test
def plotAll(dumsesRoot, ndim=1, fdir='./tmp/', io_type='binary', pro="" \
             , save=False, error=False, datdim=None, direction='x', form='png' \
             , vmind=None, vmaxd=None, vminE=None, vmaxE=None \
             , vminv=None, vmaxv=None, vminB=None, vmaxB=None):
    listdir = glob.glob(fdir + "output_*")
    listdir.sort()
    firstnum = int(listdir[0].split('_')[-1])
    lastnum  = int(listdir[-1].split('_')[-1])
    if error:
        plotResult(dumsesRoot, idump=lastnum, ndim=ndim, fdir=fdir \
           , io_type=io_type, pro=pro, save=save, datdim=datdim \
           , direction=direction, form=form, diff=True \
           , vmind=vmind, vmaxd=vmaxd, vminE=vminE, vmaxE=vmaxE \
           , vminv=vminv, vmaxv=vmaxv, vminB=vminB, vmaxB=vmaxB)
    else:
        for i in xrange(lastnum+1):
            plotResult(dumsesRoot, idump=i, ndim=ndim, fdir=fdir \
                , io_type=io_type, pro=pro, save=save, datdim=datdim \
                , direction=direction, form=form \
                , vmind=vmind, vmaxd=vmaxd, vminE=vminE, vmaxE=vmaxE \
                , vminv=vminv, vmaxv=vmaxv, vminB=vminB, vmaxB=vmaxB)
    print('\n')

## Clean tmp directory
def cleanTmp(currentDir):
    os.chdir(currentDir + '/tmp')
    for files in glob.glob("*"):
        out, st = cmd("rm -rf " + files)
    os.chdir(currentDir)

## Run tests
def runTest(currentDir, dumsesRoot, test="all", fcc="gfortran", oacc=0 \
          , io_type="'binary'", hdf='/', cdf='/', mpi='/', error=False \
          , reference=False):
    sys.path.append(dumsesRoot + 'utils/')
    sys.path = sys.path[::-1]
    import dumpy as dp

    if (mpi == '/' or mpi == ''):
        nxslice, nyslice, nzslice = "1", "1", "1"
    if(test == 'all' or test == 'shock_tube'):
        print(bold + "----------| Launching Shock tube test... |----------" + reset)
        if not (mpi == '/' or mpi == ''):
            nxslice, nyslice, nzslice = "2", "1", "1"
        launchProblem(currentDir, dumsesRoot, fcc=fcc, io_type=io_type \
                  , hdf=hdf, cdf=cdf, mpi=mpi, oacc=oacc \
                  , pro="shock_tube", tlim=".245", nx="128", ny="1", nz="1" \
                  , xmin="0.d0", xmax="1.d0", ymin="-0.5d0", ymax="0.5d0" \
                  , zmin="0.d0", zmax="1.d0", Omega0="0.d0", ciso="0.d0" \
                  , gamma="5./3.", dtdump=".0245" \
                  , nxslice = nxslice, nyslice = nyslice, nzslice = nzslice)
        if reference:
            listdir = glob.glob(currentDir + "/tmp/output_*")
            listdir.sort()
            firstnum = int(listdir[0].split('_')[-1])
            lastnum  = int(listdir[-1].split('_')[-1])
            st, out  = cmd("mv " + currentDir + "/tmp/output_%06d " %lastnum \
                         + currentDir + "/fig/reference/shock_tube/x/")
        else:
            plotAll(dumsesRoot, pro="Shock tube", save=True, vminv=-0.05 \
                  , vmaxv=0.9, io_type=eval(io_type), error=error)
        
        if not (mpi == '/' or mpi == ''):
            nxslice, nyslice, nzslice = "1", "2", "1"
        launchProblem(currentDir, dumsesRoot, ndim=3, fcc=fcc, io_type=io_type \
                  , hdf=hdf, cdf=cdf, mpi=mpi, oacc=oacc \
                  , pro="shock_tube", tlim=".245", nx="1", ny="128", nz="1" \
                  , xmin="0.d0", xmax="1.d0", ymin="0.d0", ymax="1.d0" \
                  , zmin="0.d0", zmax="1.d0", Omega0="0.d0", ciso="0.d0" \
                  , gamma="5./3.", dtdump=".0245", direction="y" \
                  , nxslice = nxslice, nyslice = nyslice, nzslice = nzslice)
        if reference:
            listdir = glob.glob(currentDir + "/tmp/output_*")
            listdir.sort()
            firstnum = int(listdir[0].split('_')[-1])
            lastnum  = int(listdir[-1].split('_')[-1])
            st, out  = cmd("mv " + currentDir + "/tmp/output_%06d " %lastnum \
                         + currentDir + "/fig/reference/shock_tube/y/")
        else:
            plotAll(dumsesRoot, pro="Shock tube", save=True, vminv=-0.05 \
                  , vmaxv=0.9, datdim=3, direction='y', io_type=eval(io_type) \
                  , error=error)
        
        if not (mpi == '/' or mpi == ''):
            nxslice, nyslice, nzslice = "1", "1", "2"
        launchProblem(currentDir, dumsesRoot, ndim=3, fcc=fcc, io_type=io_type \
                  , hdf=hdf, cdf=cdf, mpi=mpi, oacc=oacc \
                  , pro="shock_tube", tlim=".245", nx="1", ny="1", nz="128" \
                  , xmin="0.d0", xmax="1.d0", ymin="0.d0", ymax="1.d0" \
                  , zmin="0.d0", zmax="1.d0", Omega0="0.d0", ciso="0.d0" \
                  , gamma="5./3.", dtdump=".0245", direction="z" \
                  , nxslice = nxslice, nyslice = nyslice, nzslice = nzslice)
        if reference:
            listdir = glob.glob(currentDir + "/tmp/output_*")
            listdir.sort()
            firstnum = int(listdir[0].split('_')[-1])
            lastnum  = int(listdir[-1].split('_')[-1])
            st, out  = cmd("mv " + currentDir + "/tmp/output_%06d " %lastnum \
                         + currentDir + "/fig/reference/shock_tube/z/")
        else:
            plotAll(dumsesRoot, pro="Shock tube", save=True, vminv=-0.05 \
                  , vmaxv=0.9, datdim=3, direction='z', io_type=eval(io_type) \
                  , error=error)
        
    if(test == 'all' or test == 'orszag_tang'):
        print(bold + "----------| Launching Orszag Tang test... |----------" + reset)
        if not (mpi == '/' or mpi == ''):
            nxslice, nyslice, nzslice = "2", "2", "1"
        launchProblem(currentDir, dumsesRoot, ndim=2, fcc=fcc, io_type=io_type \
                  , hdf=hdf, cdf=cdf, mpi=mpi, oacc=oacc \
                  , pro="orszag_tang", tlim=".51", nx="64", ny="64", nz="1" \
                  , xmin="-.5d0", xmax=".5d0", ymin="-.5d0", ymax=".5d0" \
                  , zmin="-.5d0", zmax=".5d0", Omega0="0.d0", ciso="0.d0" \
                  , courant="0.8", dtdump=".05" \
                  , bdtypex="'periodic'", bdtypey="'periodic'" \
                  , nxslice = nxslice, nyslice = nyslice, nzslice = nzslice)
        if reference:
            listdir = glob.glob(currentDir + "/tmp/output_*")
            listdir.sort()
            firstnum = int(listdir[0].split('_')[-1])
            lastnum  = int(listdir[-1].split('_')[-1])
            st, out  = cmd("mv " + currentDir + "/tmp/output_%06d " %lastnum \
                         + currentDir + "/fig/reference/orszag_tang/x/")
        else:
            plotAll(dumsesRoot, pro="Orszag Tang", ndim=2, save=True, vmind=0.1\
                  , vmaxd=0.5, vminE=.15, vmaxE=.8, vminv=-1.3, vmaxv=1.3 \
                  , vminB=-.6, vmaxB=.6, io_type=eval(io_type), error=error)
        
        if not (mpi == '/' or mpi == ''):
            nxslice, nyslice, nzslice = "1", "2", "2"
        launchProblem(currentDir, dumsesRoot, ndim=3, fcc=fcc, io_type=io_type \
                  , hdf=hdf, cdf=cdf, mpi=mpi, oacc=oacc \
                  , pro="orszag_tang", tlim=".51", nx="1", ny="64", nz="64" \
                  , xmin="-.5d0", xmax=".5d0", ymin="-.5d0", ymax=".5d0" \
                  , zmin="-.5d0", zmax=".5d0", Omega0="0.d0", ciso="0.d0" \
                  , courant="0.8", dtdump=".05", direction='y' \
                  , bdtypex="'periodic'", bdtypey="'periodic'" \
                  , bdtypez="'periodic'" \
                  , nxslice = nxslice, nyslice = nyslice, nzslice = nzslice)
        if reference:
            listdir = glob.glob(currentDir + "/tmp/output_*")
            listdir.sort()
            firstnum = int(listdir[0].split('_')[-1])
            lastnum  = int(listdir[-1].split('_')[-1])
            st, out  = cmd("mv " + currentDir + "/tmp/output_%06d " %lastnum \
                         + currentDir + "/fig/reference/orszag_tang/y/")
        else:
            plotAll(dumsesRoot, pro="Orszag Tang", ndim=2, datdim=3 \
                  , save=True, vmind=0.1, vmaxd=0.5, vminE=.15, vmaxE=.8 \
                  , vminv=-1.3, vmaxv=1.3, vminB=-.6, vmaxB=.6 \
                  , io_type=eval(io_type), error=error, direction='y')
         
        if not (mpi == '/' or mpi == ''):
            nxslice, nyslice, nzslice = "2", "1", "2"
        launchProblem(currentDir, dumsesRoot, ndim=3, fcc=fcc, io_type=io_type \
                  , hdf=hdf, cdf=cdf, mpi=mpi, oacc=oacc \
                  , pro="orszag_tang", tlim=".51", nx="64", ny="1", nz="64" \
                  , xmin="-.5d0", xmax=".5d0", ymin="-.5d0", ymax=".5d0" \
                  , zmin="-.5d0", zmax=".5d0", Omega0="0.d0", ciso="0.d0" \
                  , courant="0.8", dtdump=".05", direction='z' \
                  , bdtypex="'periodic'", bdtypey="'periodic'" \
                  , bdtypez="'periodic'" \
                  , nxslice = nxslice, nyslice = nyslice, nzslice = nzslice)
        if reference:
            listdir = glob.glob(currentDir + "/tmp/output_*")
            listdir.sort()
            firstnum = int(listdir[0].split('_')[-1])
            lastnum  = int(listdir[-1].split('_')[-1])
            st, out  = cmd("mv " + currentDir + "/tmp/output_%06d " %lastnum \
                         + currentDir + "/fig/reference/orszag_tang/z/")
        else:
            plotAll(dumsesRoot, pro="Orszag Tang", ndim=2, datdim=3 \
                  , save=True, vmind=0.1, vmaxd=0.5, vminE=.15, vmaxE=.8 \
                  , vminv=-1.3, vmaxv=1.3, vminB=-.6, vmaxB=.6 \
                  , io_type=eval(io_type), error=error, direction='z')
        
    if(test == 'all' or test == 'magnetic_loop'):
        print(bold + "----------| Launching Magnetic loop test... |----------" + reset)
        if not (mpi == '/' or mpi == ''):
            nxslice, nyslice, nzslice = "2", "2", "1"
        launchProblem(currentDir, dumsesRoot \
                  , ndim=2, iso=1, fcc=fcc, io_type=io_type \
                  , hdf=hdf, cdf=cdf, mpi=mpi, oacc=oacc \
                  , pro="magnetic_loop", tlim="1.", nx="64", ny="64", nz="1" \
                  , xmin="-2.d0", xmax="2.d0", ymin="-2.d0", ymax="2.d0" \
                  , zmin="-2.d0", zmax="2.d0", Omega0="0.d0", ciso="1.d0" \
                  , courant="0.8", dtdump=".1" \
                  , bdtypex="'periodic'", bdtypey="'periodic'" \
                  , nxslice = nxslice, nyslice = nyslice, nzslice = nzslice)
        if reference:
            listdir = glob.glob(currentDir + "/tmp/output_*")
            listdir.sort()
            firstnum = int(listdir[0].split('_')[-1])
            lastnum  = int(listdir[-1].split('_')[-1])
            st, out  = cmd("mv " + currentDir + "/tmp/output_%06d " %lastnum \
                         + currentDir + "/fig/reference/magnetic_loop/x/")
        else:
            plotAll(dumsesRoot, pro="Magnetic loop", ndim=2, save=True \
                  , io_type=eval(io_type), vmind=1+.01e-12, vmaxd=1+2e-12 \
                  , vminE=-.5e-9, vmaxE=1.8e-9, vminv=2+.01e-13 \
                  , vmaxv=2+6.e-13, vminB=-1e-6, vmaxB=1e-6, error=error)
        
        if not (mpi == '/' or mpi == ''):
            nxslice, nyslice, nzslice = "1", "2", "2"
        launchProblem(currentDir, dumsesRoot \
                  , ndim=3, iso=1, fcc=fcc, io_type=io_type \
                  , hdf=hdf, cdf=cdf, mpi=mpi, oacc=oacc \
                  , pro="magnetic_loop", tlim="1.", nx="1", ny="64", nz="64" \
                  , xmin="-2.d0", xmax="2.d0", ymin="-2.d0", ymax="2.d0" \
                  , zmin="-2.d0", zmax="2.d0", Omega0="0.d0", ciso="1.d0" \
                  , courant="0.8", dtdump=".1", direction='y' \
                  , bdtypex="'periodic'", bdtypey="'periodic'" \
                  , bdtypez="'periodic'" \
                  , nxslice = nxslice, nyslice = nyslice, nzslice = nzslice)
        if reference:
            listdir = glob.glob(currentDir + "/tmp/output_*")
            listdir.sort()
            firstnum = int(listdir[0].split('_')[-1])
            lastnum  = int(listdir[-1].split('_')[-1])
            st, out  = cmd("mv " + currentDir + "/tmp/output_%06d " %lastnum \
                         + currentDir + "/fig/reference/magnetic_loop/y/")
        else:
            plotAll(dumsesRoot, pro="Magnetic loop", ndim=2, save=True \
                  , direction='y', io_type=eval(io_type), datdim=3 \
                  , vmind=1+.01e-12, vmaxd=1+2e-12, vminE=-.5e-9, vmaxE=1.8e-9 \
                  , vminv=2+.01e-13, vmaxv=2+6.e-13, vminB=-1e-6, vmaxB=1e-6 \
                  , error=error)
        
        if not (mpi == '/' or mpi == ''):
            nxslice, nyslice, nzslice = "2", "1", "2"
        launchProblem(currentDir, dumsesRoot \
                  , ndim=3, iso=1, fcc=fcc, io_type=io_type \
                  , hdf=hdf, cdf=cdf, mpi=mpi, oacc=oacc \
                  , pro="magnetic_loop", tlim="1.", nx="64", ny="1", nz="64" \
                  , xmin="-2.d0", xmax="2.d0", ymin="-2.d0", ymax="2.d0" \
                  , zmin="-2.d0", zmax="2.d0", Omega0="0.d0", ciso="1.d0" \
                  , courant="0.8", dtdump=".1", direction='z' \
                  , bdtypex="'periodic'", bdtypey="'periodic'" \
                  , bdtypez="'periodic'" \
                  , nxslice = nxslice, nyslice = nyslice, nzslice = nzslice)
        if reference:
            listdir = glob.glob(currentDir + "/tmp/output_*")
            listdir.sort()
            firstnum = int(listdir[0].split('_')[-1])
            lastnum  = int(listdir[-1].split('_')[-1])
            st, out  = cmd("mv " + currentDir + "/tmp/output_%06d " %lastnum \
                         + currentDir + "/fig/reference/magnetic_loop/z/")
        else:
            plotAll(dumsesRoot, pro="Magnetic loop", ndim=2, save=True \
                  , direction='z', io_type=eval(io_type), datdim=3 \
                  , vmind=1+.01e-12, vmaxd=1+2e-12, vminE=-.5e-9, vmaxE=1.8e-9 \
                  , vminv=2+.01e-13, vmaxv=2+6.e-13, vminB=-1e-6, vmaxB=1e-6 \
                  , error=error)
        
    if(test == 'all' or test == 'wind_tunnel'):
        print(bold + "----------| Launching Wind tunnel test... |----------" + reset)
        if not (mpi == '/' or mpi == ''):
            nxslice, nyslice, nzslice = "2", "2", "1"
        launchProblem(currentDir, dumsesRoot \
                  , ndim=2, iso=0, fcc=fcc, io_type=io_type \
                  , hdf=hdf, cdf=cdf, mpi=mpi, oacc=oacc \
                  , pro="wind_tunnel", tlim="100.", nx="64", ny="32", nz="1" \
                  , xmin="-2.d0", xmax="2.d0", ymin="-2.d0", ymax="2.d0" \
                  , zmin="-2.d0", zmax="2.d0", Omega0="0.d0", ciso="1.d0" \
                  , courant="0.8", dtdump="10." \
                  , bdtypex="'periodic'", bdtypey="'periodic'" \
                  , bdtypez="'periodic'" \
                  , nxslice = nxslice, nyslice = nyslice, nzslice = nzslice)
        if reference:
            listdir = glob.glob(currentDir + "/tmp/output_*")
            listdir.sort()
            firstnum = int(listdir[0].split('_')[-1])
            lastnum  = int(listdir[-1].split('_')[-1])
            st, out  = cmd("mv " + currentDir + "/tmp/output_%06d " %lastnum \
                         + currentDir + "/fig/reference/wind_tunnel/x/")
        else:
            plotAll(dumsesRoot, pro="Wind tunnel", ndim=2, save=True \
                  , vmind=0, vmaxd=30, vminE=0, vmaxE=1, vminv=-.03, vmaxv=.03 \
                  , vminB=0, vmaxB=0, io_type=eval(io_type), error=error)
        
        if not (mpi == '/' or mpi == ''):
            nxslice, nyslice, nzslice = "1", "2", "2"
        launchProblem(currentDir, dumsesRoot \
                  , ndim=3, iso=0, fcc=fcc, io_type=io_type \
                  , hdf=hdf, cdf=cdf, mpi=mpi, oacc=oacc \
                  , pro="wind_tunnel", tlim="100.", nx="1", ny="64", nz="32" \
                  , xmin="-2.d0", xmax="2.d0", ymin="-2.d0", ymax="2.d0" \
                  , zmin="-2.d0", zmax="2.d0", Omega0="0.d0", ciso="1.d0" \
                  , courant="0.8", dtdump="10.", direction='y' \
                  , bdtypex="'periodic'", bdtypey="'periodic'" \
                  , bdtypez="'periodic'" \
                  , nxslice = nxslice, nyslice = nyslice, nzslice = nzslice)
        if reference:
            listdir = glob.glob(currentDir + "/tmp/output_*")
            listdir.sort()
            firstnum = int(listdir[0].split('_')[-1])
            lastnum  = int(listdir[-1].split('_')[-1])
            st, out  = cmd("mv " + currentDir + "/tmp/output_%06d " %lastnum \
                         + currentDir + "/fig/reference/wind_tunnel/y/")
        else:
            plotAll(dumsesRoot, pro="Wind tunnel", ndim=2, save=True, datdim=3 \
                  , direction='y', vmind=0, vmaxd=30, vminE=0, vmaxE=1 \
                  , vminv=-.03, vmaxv=.03, vminB=0, vmaxB=0 \
                  , io_type=eval(io_type), error=error)
        
        if not (mpi == '/' or mpi == ''):
            nxslice, nyslice, nzslice = "2", "1", "2"
        launchProblem(currentDir, dumsesRoot \
                  , ndim=3, iso=0, fcc=fcc, io_type=io_type \
                  , hdf=hdf, cdf=cdf, mpi=mpi, oacc=oacc \
                  , pro="wind_tunnel", tlim="100.", nx="32", ny="1", nz="64" \
                  , xmin="-2.d0", xmax="2.d0", ymin="-2.d0", ymax="2.d0" \
                  , zmin="-2.d0", zmax="2.d0", Omega0="0.d0", ciso="1.d0" \
                  , courant="0.8", dtdump="10.", direction='z' \
                  , bdtypex="'periodic'", bdtypey="'periodic'" \
                  , bdtypez="'periodic'" \
                  , nxslice = nxslice, nyslice = nyslice, nzslice = nzslice)
        if reference:
            listdir = glob.glob(currentDir + "/tmp/output_*")
            listdir.sort()
            firstnum = int(listdir[0].split('_')[-1])
            lastnum  = int(listdir[-1].split('_')[-1])
            st, out  = cmd("mv " + currentDir + "/tmp/output_%06d " %lastnum \
                         + currentDir + "/fig/reference/wind_tunnel/z/")
        else:
            plotAll(dumsesRoot, pro="Wind tunnel", ndim=2, save=True, datdim=3 \
                  , direction='z', vmind=0, vmaxd=30, vminE=0, vmaxE=1 \
                  , vminv=-.03, vmaxv=.03, vminB=0, vmaxB=0 \
                  , io_type=eval(io_type), error=error)
    
    if(test == 'all' or test == 'shearing'):
        print(bold + "----------| Launching shearing box tests... |----------" + reset)
        if not (mpi == '/' or mpi == ''):
            nxslice, nyslice, nzslice = "2", "2", "1"
        launchProblem(currentDir, dumsesRoot, ndim=3, fcc=fcc, io_type=io_type \
                  , hdf=hdf, cdf=cdf, mpi=mpi, oacc=oacc, iso=1 \
                  , pro="mri", tlim="628.3", nx="32", ny="32", nz="32" \
                  , xmin="-.5d0", xmax=".5d0", ymin="0.d0", ymax="4.d0" \
                  , zmin="-.5d0", zmax=".5d0", Omega0="1.d-3", ciso="1.d-3" \
                  , courant="0.8", dtdump="62." \
                  , bdtypex="'shearingbox'", bdtypey="'periodic'" \
                  , bdtypez="'periodic'" \
                  , nxslice = nxslice, nyslice = nyslice, nzslice = nzslice \
                  , beta="400.d0", type_="Radial", amp="0.d0")
        if reference:
            listdir = glob.glob(currentDir + "/tmp/output_*")
            listdir.sort()
            firstnum = int(listdir[0].split('_')[-1])
            lastnum  = int(listdir[-1].split('_')[-1])
            st, out  = cmd("mv " + currentDir + "/tmp/output_%06d " %lastnum \
                         + currentDir + "/fig/reference/shearing/mri/")
        else:
            plotAll(dumsesRoot, pro="Radial shearing", ndim=2, save=True, direction="x"
                  , vmind=9.976e-1, vmaxd=1.0024, vminE=-2.4e-6, vmaxE=2.4e-6 \
                  , vminv=-1.6e-6, vmaxv=1.6e-6, vminB=-.6e-4, vmaxB=.6e-4 \
                  , io_type=eval(io_type), error=error)

        if not (mpi == '/' or mpi == ''):
            nxslice, nyslice, nzslice = "2", "2", "1"
        launchProblem(currentDir, dumsesRoot, ndim=3, fcc=fcc, io_type=io_type \
                  , hdf=hdf, cdf=cdf, mpi=mpi, oacc=oacc, iso=0 \
                  , pro="incomp_shearwave", tlim="1000" \
                  , nx="64", ny="64", nz="2" \
                  , xmin="-2.5d-1", xmax="2.5d-1" \
                  , ymin="-2.5d-1", ymax="2.5d-1" \
                  , zmin="-2.5d-1", zmax="2.5d-1" \
                  , Omega0="1.d-3", ciso="1.d-3" \
                  , courant="0.8", dtdump="100." \
                  , bdtypex="'shearingbox'", bdtypey="'periodic'" \
                  , bdtypez="'periodic'" \
                  , nxslice = nxslice, nyslice = nyslice, nzslice = nzslice)
        if reference:
            listdir = glob.glob(currentDir + "/tmp/output_*")
            listdir.sort()
            firstnum = int(listdir[0].split('_')[-1])
            lastnum  = int(listdir[-1].split('_')[-1])
            st, out  = cmd("mv " + currentDir + "/tmp/output_%06d " %lastnum \
                         + currentDir + "/fig/reference/shearing/ishwave/")
            print st, out
        else:
            plotAll(dumsesRoot, pro="Incompressible shearing wave", ndim=2 \
                  , save=True, direction="x"
                  , vmind=0.97539, vmaxd=1.03128, vminE=-1.6e-13, vmaxE=1.6e-13\
                  , vminv=-7.e-7, vmaxv=7.e-7, io_type=eval(io_type) \
                  , error=error)

        if not (mpi == '/' or mpi == ''):
            nxslice, nyslice, nzslice = "2", "2", "1"
        launchProblem(currentDir, dumsesRoot, ndim=3, fcc=fcc, io_type=io_type \
                  , hdf=hdf, cdf=cdf, mpi=mpi, oacc=oacc, iso=0 \
                  , pro="shearwave", tlim="1000" \
                  , nx="128", ny="64", nz="2" \
                  , xmin="-2.d-1", xmax="2.d-1" \
                  , ymin="-2.d-1", ymax="2.d-1" \
                  , zmin="-2.5d-1", zmax="2.5d-1" \
                  , Omega0="1.d-3", ciso="1.d-3" \
                  , courant="0.8", dtdump="100." \
                  , bdtypex="'shearingbox'", bdtypey="'periodic'" \
                  , bdtypez="'periodic'" \
                  , nxslice = nxslice, nyslice = nyslice, nzslice = nzslice)
        if reference:
            listdir = glob.glob(currentDir + "/tmp/output_*")
            listdir.sort()
            firstnum = int(listdir[0].split('_')[-1])
            lastnum  = int(listdir[-1].split('_')[-1])
            st, out  = cmd("mv " + currentDir + "/tmp/output_%06d " %lastnum \
                         + currentDir + "/fig/reference/shearing/shwave")
        else:
            plotAll(dumsesRoot, pro="Shearing wave", ndim=2 \
                  , save=True, direction="x"
                  , vmind=0.9842, vmaxd=1.016, vminE=-5.5e-14, vmaxE=5.5e-14\
                  , vminv=-4.13e-8, vmaxv=4.13e-8, io_type=eval(io_type) \
                  , error=error)

    if(test == 'all' or test == 'shearing' or test == 'mri'):
        print(bold + "----------| Launching MRI test... |----------" + reset)
        if not (mpi == '/' or mpi == ''):
            nxslice, nyslice, nzslice = "1", "2", "1"
        launchProblem(currentDir, dumsesRoot, ndim=3, fcc=fcc, io_type=io_type \
                  , hdf=hdf, cdf=cdf, mpi=mpi, oacc=oacc, iso=1 \
                  , pro="mri", tlim="60000" \
                  , nx="32", ny="64", nz="32" \
                  , xmin="-0.5d0", xmax="0.5d0" \
                  , ymin="0.d0", ymax="4.d0" \
                  , zmin="-0.5d0", zmax="0.5d0" \
                  , Omega0="1.d-3", ciso="1.d-3" \
                  , courant="0.8", dtdump="6000.", dthist="100.d0" \
                  , bdtypex="'shearingbox'", bdtypey="'periodic'" \
                  , bdtypez="'periodic'" \
                  , nxslice = nxslice, nyslice = nyslice, nzslice = nzslice \
                  , beta="1000.d0", type_="Test", amp="0.d0")
        if reference:
            listdir = glob.glob(currentDir + "/tmp/output_*")
            listdir.sort()
            firstnum = int(listdir[0].split('_')[-1])
            lastnum  = int(listdir[-1].split('_')[-1])
            st, out  = cmd("mv " + currentDir + "/tmp/history.txt " \
                         + currentDir + "/fig/reference/shearing/mri/")
            print st, out
        else:
            dp.compHistory(['maxwell', 'reynolds', 'dt', 'magp'], currentDir + '/fig/reference/shearing/mri/', 'tmp/', loc='upper right', tmax=60000, save=True)
            st, out  = cmd("mv " + currentDir + "/compHist_MaxReyDtMag.* " \
                         + currentDir + "/fig/.")
            print st, out
            
    cleanTmp(currentDir)

## Create PDF file from produced output figures
def makePDF(currentDir, error=False):
    os.chdir(currentDir + '/fig/')
    if error:
        st, out = cmd('pdflatex test_err_results.tex')
        st, out = cmd('rm test_err_results.{aux,log}')
    else:
        st, out = cmd('pdflatex test_results.tex')
        st, out = cmd('rm test_results.{aux,log}')
    os.chdir(currentDir)

## This is the main function
def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--test=", "-t", dest="test", type=str, default="all" \
                      , help="Test to run in [shock_tube, wind_tunnel, magnetic_loop, orszag_tang, shearing, all]. If 'all', run all the tests. Default: 'all'")
    parser.add_argument("--with-fortran-compiler=", "-f", dest="fcc", type=str \
                      , default="gfortran", help="Fortran compiler to compile the code. Default: 'gfortran'")
    parser.add_argument("--with-openacc", "-o", dest="oacc", action="store_true" \
                      , help="To compile the code for accelerators. Default: deactivated")
    parser.add_argument("--io-type=", "-i", dest="io_type", type=str \
                      , default="'binary'", help="I/O type to use, in [binary, hdf5, phdf5, pnetcdf]. Default: 'binary'")
    parser.add_argument("--with-phdf5=", "-H", dest="hdf", type=str \
                      , default="/", help="Path of the (Parallel)HDF5 library. Default: '/'")
    parser.add_argument("--with-pnetcdf=", "-c", dest="cdf", type=str \
                      , default="/", help="Path of the ParallelNetCDF library. Default: '/'")
    parser.add_argument("--with-mpi=", "-m", dest="mpi", type=str \
                      , default="", help="Path of the MPI library. Default: None")
    parser.add_argument("--with-openmp", "-O", dest="openmp" \
                      , action="store_true", help="Run with OpenMP")
    parser.add_argument("--error=", "-e", dest="error", action="store_true" \
                      , help="Print error(current run, reference run)")
    parser.add_argument("--reference=", "-r", dest="reference" \
                      , action="store_true", help="Generate reference output")
    parser.add_argument("--dumses-root", "-d", dest="dumsesroot", default=None \
                      , help="Define the absolute path of your current DUMSES-Hybrid installation. Default: None")
    args = parser.parse_args()
    test, fcc, oacc, io_type = args.test, args.fcc, args.oacc, args.io_type
    hdf, cdf, mpi, openmp    = args.hdf, args.cdf, args.mpi, args.openmp
    error, reference         = args.error, args.reference
    dumsesRoot               = args.dumsesroot

    if dumsesRoot:
        currentDir = dumsesRoot + "/utils/test/"
        os.chdir(currentDir)
    else:
        st, currentDir = cmd('pwd')
        try:
            assert("/".join(currentDir.split('/')[-2:] ) == "utils/test")
        except AssertionError:
            print("You are not launching test.py script in the right directory, nor you defined the DUMSES-Hybrid directory (with -d option).")
            sys.exit(2)
        dumsesRoot = currentDir + '/../../'
    
    try:
        os.makedirs(currentDir + './tmp')
    except OSError:
        pass

    if mpi and fcc == 'gfortran':
        fcc = ''

    if openmp:
        os.environ['OMP_NUM_THREADS'] = '4'
    else:
        os.environ['OMP_NUM_THREADS'] = '1'

    oacc = (0 if not(oacc) else 1)

    runTest(currentDir, dumsesRoot, test=test, fcc=fcc, io_type=io_type \
                , oacc=oacc, hdf=hdf, cdf=cdf, mpi=mpi, error=error \
                , reference=reference)
    if test == 'all':
        makePDF(currentDir, error=error)

if __name__ == "__main__":
    main()
