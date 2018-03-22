#!/usr/bin/python
#===============================================================================
## \file submit.py
# \brief
# \b DUMSES-Hybrid:
# This is submition script for BlueGene/Q (with LoadLeveler)
# \author
# Marc Joos <marc.joos@cea.fr>
# \copyright
# Copyrights 2013-2015, CEA.
# This file is distributed under the CeCILL-A & GNU/GPL licenses, see
# <http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html> and
# <http://www.gnu.org/licenses/>
# \date
# \b created:          05-23-2014 
# \b last \b modified: 05-26-2014

#===============================================================================
import os, sys, glob
if sys.version_info[0] <= 2:
    from commands import getstatusoutput as cmd
else:
    from subprocess import getstatusoutput as cmd
if sys.version_info[1] >= 7:
    import argparse
else:
    st, home = cmd("echo $HOME")
    sys.path.append(home + '/dumses_hybrid/utils/pyutils/')
    import argparse

bold    = "\033[1m"
reset   = "\033[0;0m"

## Create input for the given problem
def createInput(inputfile="input", problem="magnetic_loop" \
              , tlim="1.", verbose=".false.", debug=".false." \
              , bdtypex="'zero_gradient'", bdtypey="'zero_gradient'" \
              , bdtypez="'zero_gradient'", riemann="'hlld'", riemann2d="'hlld'"\
              , slope_type="2", courant="0.7", nx="32", ny="32", nz="32" \
              , xmin="-1.d0", xmax="1.d0", ymin="-1.d0", ymax="1.d0" \
              , zmin="0.d0", zmax="1.d0", Omega0="0.d0", ciso="0.d0" \
              , gamma="1.001d0", dtdump=".1", io_type="'binary'" \
              , nxslice="1", nyslice="1", nzslice="1", direction='x'):
    st, out = cmd('cp input.template ' + inputfile)
    template = open(inputfile, 'rt').read()
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
            "io_type": io_type,
            "nxslice": nxslice,
            "nyslice": nyslice,
            "nzslice": nzslice
            }

    with open(inputfile, 'wt') as output:
        output.write(template %data)

    lproDirect = ['shock_tube', 'orszag_tang', 'magnetic_loop', 'wind_tunnel']
    if problem in lproDirect:
        f = open(inputfile, 'a')
        initparam="\n&init_params\n        direction = '%s'\n        /\n" \
            %direction
        f.write(initparam)
        f.close()

## Check if the current executable is compatible with the input file
def checkExe(problem="magnetic_loop", exedir="$HOME/dumses_hybrid/bin/", force=False):
    st, out = cmd("grep PRBDIR " + exedir + "../Makefile")
    probMakefile = out.split()[2].split("/")[-1]
    if not force:
        if probMakefile != problem:
            print(bold + "Warning! " + reset + "The problem specified is different from the problem in " + exedir + "../Makefile.\n If you know what you are doing, please reexecute the script with the '-f, --force' option")
            sys.exit(1)

## Create jobfile for the given problem
def createJobLL(jobfile="job.ll", inputfile="input", jobname="test_hybrid" \
                  , wallclock="30:00", rundir="$WORKDIR/loop" \
                  , exedir="$HOME/dumses_hybrid/bin/", bg_size="64" \
                  , rankpernode="16", omp_num_threads="1"):
    st, out = cmd('cp job_ll.template ' + jobfile)
    template = open(jobfile, 'rt').read()
    data = {"jobname": jobname,
            "wallclock": wallclock,
            "bg_size": bg_size,
            "rankpnode": rankpernode,
            "omp_num_threads": omp_num_threads,
            "rundir": rundir,
            "exedir": exedir,
            "inputfile": inputfile,
            "jobfile": jobfile
            }

    with open(jobfile, 'wt') as output:
        output.write(template %data)

## Custom types to format entries in Fortran format
def FortranFloat(value):
    vfloat = float(value)
    vfloat = str(vfloat) + "d0"
    return vfloat

def FortranInteger(value):
    vint = int(value)
    vint = str(vint)
    return vint

## This is the main function
def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--inputfile=", dest="inputfile", type=str, default="default")
    parser.add_argument("--problem=", dest="problem", type=str, default="magnetic_loop")
    parser.add_argument("--tlim=", dest="tlim", type=FortranFloat, default="0.101")
    parser.add_argument("--verbose=", dest="verbose", type=str, default=".false.")
    parser.add_argument("--debug=", dest="debug", type=str, default=".false.")
    parser.add_argument("--bdtypex=", dest="bdtypex", type=str, default="'periodic'")
    parser.add_argument("--bdtypey=", dest="bdtypey", type=str, default="'periodic'")
    parser.add_argument("--bdtypez=", dest="bdtypez", type=str, default="'periodic'")
    parser.add_argument("--riemann=", dest="riemann", type=str, default="'hlld'")
    parser.add_argument("--Riemann2d=", dest="riemann2d", type=str, default="'hlld'")
    parser.add_argument("--slope_type=", dest="slope_type", type=FortranInteger, default="2")
    parser.add_argument("--courant=", dest="courant", type=FortranFloat, default="0.8")
    parser.add_argument("--Nx=", dest="nx", type=FortranInteger, default="32")
    parser.add_argument("--Ny=", dest="ny", type=FortranInteger, default="32")
    parser.add_argument("--Nz=", dest="nz", type=FortranInteger, default="32")
    parser.add_argument("--xmin=", dest="xmin", type=FortranFloat, default="-2.0")
    parser.add_argument("--xmax=", dest="xmax", type=FortranFloat, default="2.0")
    parser.add_argument("--ymin=", dest="ymin", type=FortranFloat, default="-2.0")
    parser.add_argument("--ymax=", dest="ymax", type=FortranFloat, default="2.0")
    parser.add_argument("--zmin=", dest="zmin", type=FortranFloat, default="-2.0")
    parser.add_argument("--zmax=", dest="zmax", type=FortranFloat, default="2.0")
    parser.add_argument("--Omega0=", dest="Omega0", type=FortranFloat, default="0.")
    parser.add_argument("--ciso=", dest="ciso", type=FortranFloat, default="1.")
    parser.add_argument("--gamma=", dest="gamma", type=FortranFloat, default="1.001")
    parser.add_argument("--dtdump=", dest="dtdump", type=FortranFloat, default="0.1")
    parser.add_argument("--io_type=", dest="io_type", type=str, default="'phdf5'")
    parser.add_argument("--nxslice=", dest="nxslice", type=FortranInteger, default="8")
    parser.add_argument("--nyslice=", dest="nyslice", type=FortranInteger, default="16")
    parser.add_argument("--nzslice=", dest="nzslice", type=FortranInteger, default="16")
    parser.add_argument("--direction=", dest="direction", type=str, default="x")
    parser.add_argument("--scheduler", dest="scheduler", type=str, default="loadleveler")
    parser.add_argument("--jobfile=", dest="jobfile", type=str, default="default")
    parser.add_argument("--jobname=", dest="jobname", type=str, default="test_hybrid")
    parser.add_argument("--wallclock=", dest="wallclock", type=str, default="30:00")
    parser.add_argument("--rundir=", dest="rundir", type=str, default="$WORKDIR/loop")
    parser.add_argument("--exedir=", dest="exedir", type=str, default="$HOME/dumses_hybrid/bin/")
    parser.add_argument("--bg_size=", dest="bg_size", type=FortranInteger, default="64")
    parser.add_argument("--rankpernode=", dest="rankpernode", type=FortranInteger, default="16")
    parser.add_argument("--omp_num_threads=", dest="omp_num_threads", type=FortranInteger, default="1")
    parser.add_argument("-f", "--force", dest="force", action="store_true")
    args = parser.parse_args()

    inputfile, problem        = args.inputfile, args.problem
    tlim, verbose, debug      = args.tlim, args.verbose, args.debug
    bdtypex, bdtypey, bdtypez = args.bdtypex, args.bdtypey, args.bdtypez
    riemann, riemann2d        = args.riemann, args.riemann2d
    slope_type, courant       = args.slope_type, args.courant
    nx, ny, nz                = args.nx, args.ny, args.nz
    xmin, ymin, zmin          = args.xmin, args.ymin, args.zmin
    xmax, ymax, zmax          = args.xmax, args.ymax, args.zmax
    Omega0, ciso, gamma       = args.Omega0, args.ciso, args.gamma
    dtdump, io_type           = args.dtdump, args.io_type
    nxslice, nyslice, nzslice = args.nxslice, args.nyslice, args.nzslice
    direction, scheduler      = args.direction, args.scheduler
    jobfile, jobname          = args.jobfile, args.jobname
    wallclock, rundir, exedir = args.wallclock, args.rundir, args.exedir
    bg_size, rankpernode      = args.bg_size, args.rankpernode
    omp_num_threads, force    = args.omp_num_threads, args.force

    # Check domain decomposition & machine configuration
    try:
        assert int(nxslice)*int(nyslice)*int(nzslice) <= int(bg_size)*int(rankpernode)*int(omp_num_threads)
    except AssertionError:
        print(bold + "   Error!\n" + reset + "Domain decomposition does not fit machine configuration;\n nxslice*nyslice*nzslice != bg_size*rankpernode")
        sys.exit(1)

    # Create extension depending on scheduler
    if scheduler == "loadlever": ext = ".ll"

    # Create input and job file names
    if inputfile == "default":
        inputfile = "input_%sx%s_%sx%sx%s_%s" %(bg_size,rankpernode,nxslice,nyslice,nzslice,omp_num_threads)
    if jobfile == "default":
        jobfile = "job_%sx%s_%sx%sx%s_%s" %(bg_size,rankpernode,nxslice,nyslice,nzslice,omp_num_threads) + ext

    checkExe(problem=problem, exedir=exedir, force=force)

    # Create input
    createInput(inputfile=inputfile, problem=problem \
                    , tlim=tlim, verbose=verbose, debug=debug \
                    , bdtypex=bdtypex, bdtypey=bdtypey, bdtypez=bdtypez \
                    , riemann=riemann, riemann2d=riemann2d \
                    , slope_type=slope_type, courant=courant \
                    , nx=nx, ny=ny, nz=nz, xmin=xmin, xmax=xmax \
                    , ymin=ymin, ymax=ymax, zmin=zmin, zmax=zmax \
                    , Omega0=Omega0, ciso=ciso, gamma=gamma \
                    , dtdump=dtdump, io_type=io_type \
                    , nxslice=nxslice, nyslice=nyslice, nzslice=nzslice \
                    , direction=direction)

    # Create job
    if scheduler == "loadleveler":
        createJobLL(jobfile=jobfile, inputfile=inputfile, jobname=jobname \
                  , wallclock=wallclock, rundir=rundir \
                  , exedir=exedir, bg_size=bg_size \
                  , rankpernode=rankpernode, omp_num_threads=omp_num_threads)

    # Submit job
    if scheduler == "loadleveler":
        st, out = cmd("llsubmit " + jobfile)
        print(out)

if __name__ == "__main__":
    main()
