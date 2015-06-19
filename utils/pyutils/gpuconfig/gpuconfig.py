#!/usr/bin/env python
#===============================================================================
# DUMSES-Hybrid:
# DUMSES GPU configuration test script.
# This script can be used to check the best configuration (in term of number of
# registers and of kernel configuration) for the most critical DUMSES routines.
# 
# author:
# Marc Joos <marc.joos@cea.fr>, Sebastien Fromang, Romain Teyssier, 
# Patrick Hennebelle
# copyright:
# Copyrights 2013-2015, CEA
# This file is distributed under the CeCILL-A & GNU/GPL licenses, see
# <http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html> and
# <http://www.gnu.org/licenses/>
# date:
# created:       04-07-2015
# last modified: 04-21-2015
#===============================================================================
import re
import os
import sys
import shutil
import pickle
import operator
import pylab as pl
from functools import partial
from collections import OrderedDict
from commands import getstatusoutput as cmd

bold    = "\033[1m"
reset   = "\033[0;0m"

class GPUTest(object):
    """Class to manage GPU configuration test.
    Usage:
       - initialization:
       data = GPUTest([rootdir, moddir, filename, makefile, mrc, nxy])

    Optional arguments:
       rootdir:  directory of your local DUMSES sources
       mrc:      maximum number of registers used
       nxy:      number of configurations tested (ex: if nxy=4, the kernels confiugation will be in [4, 8, 16, 32])
       moddir:   directory of the DUMSES modules
       filename: filename of the DUMSES module to modify
       makefile: filename of the Makefile to modify

       You should probably not modify moddir, filename and makefile (unless you know what's you are doing).

       - launch test:
       data.launchTest()

       - if you have a logfile available but no picklefile with the results:
       data.retrieveFromLogFile()

       - if you want to plot the result of a given test:
       data.plotTimeConfig()
       
       """
    def __init__(self, rootdir="/local/home/mjoos/codes/dumses/dumses_hybrid/" \
                     , moddir="./src/modules/" \
                     , filename="./commons.f90" \
                     , makefile="./Makefile" \
                     , mrc=64 \
                     , nxy=4):
        self.rootdir    = rootdir
        self.moddir     = moddir
        self.filename   = filename
        self.makefile   = makefile
        self.picklefile = "save_time_rc{}.pkl".format(mrc)
        self.logfile    = self.picklefile.replace(".pkl", ".dat").replace('save_time', 'log')
        self.allfile    = self.picklefile.replace("save_time", "save_all")
        self.blockDict  = {'blockx_solver': 32,	   
                           'blocky_solver': 4,	   
                           'blockx_solver_mag': 32, 
                           'blocky_solver_mag': 4, 
                           'blockx_trace': 32,	   
                           'blocky_trace': 4,	   
                           'blockx_update': 32,	   
                           'blocky_update': 4,	   
                           'blockx_slope': 32,	   
                           'blocky_slope': 4,	   
                           'blockx_primitive': 32,  
                           'blocky_primitive': 4}
        self.listTime   = ['slope', 'trace', 'riemann', 'riemann magnetic',
                           'update', 'primitive', 'elptime', 'cputime']
        self.cresults   = {}
        self.timeDict   = {}
        self.mrc = mrc
        self.nxy = nxy

    def rewriteFile(self):
        """Rewrite block variables in a module file."""
        filename = self.rootdir + self.moddir + self.filename
        tmpfile  = filename.replace('f90', 'tmp')

        block = re.compile(r'(^.*)(block[xy]\_)(\w.*)(\s*=\s*)(\d\d*)(\s*)($)', re.IGNORECASE)
    
        with open(filename, 'r') as f, open(tmpfile, 'w') as t:
            for line in f:
                if block.match(line):
                    line = "".join(block.match(line).group(1,2,3,4)) \
                   + str(self.blockDict["".join(block.match(line).group(2,3))]) \
                   + block.match(line).group(6)
                t.write(line)
    
        shutil.move(tmpfile, filename)

    def rewriteInput(self, inputdir, nmdim=128, nmslice=1):
        "Rewrite input file"
        filename = inputdir + "/input"
        tmpfile = filename + "tmp"

        verbose = re.compile(r'(^\s*[,\s*]\s*)(verbose\s*\=\s*)(\.\w+\.)', re.IGNORECASE)
        ndim    = re.compile(r'(^\s*[,\s*]\s*)(n\w\s*\=\s*)(\d+)', re.IGNORECASE)
        nslice  = re.compile(r'(^\s*[,\s*]\s*)(n\wslice\s*\=\s*)(\d+)', re.IGNORECASE)

        with open(filename, 'r') as f, open(tmpfile, 'w') as t:
            for line in f:
                if verbose.match(line):
                    line = "".join(verbose.match(line).group(1,2)) + ".true.\n"
                if ndim.match(line):
                    line = "".join(ndim.match(line).group(1,2)) + str(nmdim) + "\n"
                if nslice.match(line):
                    line = "".join(nslice.match(line).group(1,2)) + str(nmslice) + "\n"
                t.write(line)
        shutil.move(tmpfile, filename)

    def replaceStringMRC(self, mo, maxregcount):
        """Replace 'maxregcount' value in a string."""
        return mo.group(1) + str(maxregcount) + mo.group(3)

    def rewriteMakefile(self):
        """Rewrite DUMSES Makefile to change 'maxregcount' CUDA variable."""
        filename = self.rootdir + self.makefile
        tmpfile  = filename + ".tmp"

        replaceString = partial(self.replaceStringMRC, maxregcount=self.mrc)
        mrc = re.compile(r'(^.*maxregcount:)(\d*)(.*$)', re.IGNORECASE)
    
        with open(filename, 'r') as f, open(tmpfile, 'w') as t:
            for line in f:
                if mrc.match(line):
                    line = mrc.sub(replaceString, line)
                t.write(line)
        
        shutil.move(tmpfile, filename)

    def geneDict(self, x=32, y=4):
        """Fill blockDict dictionary with (x,y) couple."""
        blockx = re.compile('(blockx_)(\w*)', re.IGNORECASE)
        blocky = re.compile('(blocky_)(\w*)', re.IGNORECASE)
        
        for key in self.blockDict.iterkeys():
            if blockx.match(key):
                self.blockDict[key] = x
            elif blocky.match(key):
                self.blockDict[key] = y

    def geneCouple(self):
        """Generate a list with (x,y) couple from 4 to 2^n."""
        try:
            assert(self.nxy >= 1)
        except AssertionError:
            print("n must be strictly larger than 1")
            sys.exit(2)
        x = [2**i for i in xrange(2,self.nxy+2)]
        y = x
        couple = []
        for xitem in x:
            for yitem in y:
                couple.append((xitem, yitem)) 
        return couple

    def getTimes(self, list_):
        """Retrieve elapsed time from a list of lines (from DUMSES run output)."""
        ElapsedTime = re.compile('(^\s*Elapsed time.*)(\d+\.\d+E\+\d+)(.*)')
        CPUTime = re.compile('(^\s*CPU time\s+:\s+)(\d+\.\d+E\+\d+)(.*)')
    
        slope = re.compile('(^slope.+)(\d+\.\d+)')
        trace = re.compile('(^trace.+)(\d+\.\d+)')
        riemann = re.compile('(^riemann:\s*)(?!magnetic)(\d+\.\d+)')
        riemannMagnetic = re.compile('(^riemann magnetic.+)(\d+\.\d+)')
        update = re.compile('(^update.+)(\d+\.\d+)')
        primitive = re.compile('(^primitive.+)(\d+\.\d+)')
        
        results = {key: [] for key in self.listTime}
        for line in list_:
            if ElapsedTime.match(line):
                results['elptime'].append(float(ElapsedTime.match(line).group(2)))
            if CPUTime.match(line):
                results['cputime'].append(float(CPUTime.match(line).group(2)))
            if slope.match(line):
                results['slope'].append(float(slope.match(line).group(2)))
            if trace.match(line):
                results['trace'].append(float(trace.match(line).group(2)))
            if riemann.match(line):
                results['riemann'].append(float(riemann.match(line).group(2)))
            if riemannMagnetic.match(line):
                results['riemann magnetic'].append(float(riemannMagnetic.match(line).group(2)))
            if update.match(line):
                results['update'].append(float(update.match(line).group(2)))
            if primitive.match(line):
                results['primitive'].append(float(primitive.match(line).group(2)))
        return results

    def cutLogFile(self):
        """Cut log file in independent runs."""
        cr = re.compile(r'^\(c\)')
        
        ltest = []
        with open(self.logfile, 'r') as f:
            for i, line in enumerate(f):
                if cr.match(line):
                    ltest.append(i)
        return ltest
    
    def extractFromLogFile(self):
        """Extract elapsed times from log file."""
        ltest = self.cutLogFile()
        f = open(self.logfile)
        lines = f.readlines()
        f.close()
    
        try:
            assert(len(ltest) == len(self.geneCouple()))
        except AssertionError:
            print("Something is wrong with the log file (not enough tests run)")
            print len(ltest), len(geneCouple(nxy))
            sys.exit(2)
    
        i = 0
        for couple in self.geneCouple():
            currentKey = '{}x{}'.format(couple[0], couple[1])
            self.cresults[currentKey] = {}
            if i < len(ltest)-1:
                self.cresults[currentKey] = self.getTimes(lines[ltest[i]:ltest[i+1]])
            else:
                self.cresults[currentKey] = self.getTimes(lines[ltest[i]:])
            i += 1 
    
            with open(self.allfile, 'a') as f:
                pickle.dump({currentKey: self.cresults[currentKey]}, f)
    
            for key in self.cresults[currentKey].iterkeys():
                self.cresults[currentKey][key] = sum(self.cresults[currentKey][key])/len(self.cresults[currentKey][key])
    
        with open(self.picklefile, 'w') as f:
            pickle.dump(self.cresults, f)
    
    def launchTest(self):
        """Launch a GPU configuration test for a given maxregcount and a maximum block dimension. Store the log of the run in 'log_rc[maxregcount].dat', all the elapsed times in 'save_all_rc[maxregcount].pkl' and the reduced elapsed times in 'save_time_rc[maxregcount].pkl'. Return a dictionary with the reduced results."""
        firstPass = True
        for couple in self.geneCouple():
            currentKey = '{}x{}'.format(couple[0], couple[1])
            self.cresults[currentKey] = {}
    
            print("Computing time for {} with maxregcount={}".format(currentKey, self.mrc))
    
            self.geneDict(x=couple[0], y=couple[1])
            self.rewriteFile()
            st, out = cmd('pwd')
    
            localdir = out
            os.chdir(self.rootdir)
            if firstPass:
                st, out = cmd('./configure --problem=magnetic_loop -o')
                self.rewriteMakefile()
                st, out = cmd('./make.py')
                st, out = cmd('cp src/problem/magnetic_loop/input ' + localdir + '/.')
                self.rewriteInput(localdir)
            else:
                st, out = cmd('make shallow_clean; make')
    
            st, out = cmd('cp bin/dumses ' + localdir + '/.')
            os.chdir(localdir)
            st, out = cmd('./dumses')
    
            if firstPass:
                with open(self.logfile, 'w') as f:
                    f.write(out)
            else:
                with open(self.logfile, "a") as f:
                    f.write("\n")
                    f.write(out)
    
            if firstPass:
                firstPass = False
    
            self.cresults[currentKey] = self.getTimes(out.split('\n'))
    
            with open(self.allfile, 'a') as f:
                pickle.dump({currentKey: self.cresults[currentKey]}, f)
    
            for key in self.cresults[currentKey].iterkeys():
                if self.cresults[currentKey][key]:
                    self.cresults[currentKey][key] = sum(self.cresults[currentKey][key])/len(self.cresults[currentKey][key])
                else:
                    self.cresults[currentKey][key] = None
    
        with open(self.picklefile, 'w') as f:
            pickle.dump(self.cresults, f)
    
    def revertTimeDict(self, results):
        """Reorder reduced results dictionary, so the first entry is the name of the timed subroutines and the second entry the kernel configuration."""
        ndict = {}
        timedSubroutines = results[results.keys()[0]].keys()
        for tSub in timedSubroutines:
            ndict[tSub] = OrderedDict()
            lkeys = []
            for couple in self.geneCouple():
                lkeys.append("{}x{}".format(couple[0], couple[1]))
            for key in lkeys:
                ndict[tSub][key] = results[key][tSub]
        return ndict
    
    def plotTimeConfig(self, *args, **kwargs):
        """Plot all the elapsed time for a reduced results dictionary"""
        rdict = self.revertTimeDict(self.cresults)
        for tkey in rdict.iterkeys():
            if tkey not in ['elptime', 'cputime']:
                pl.plot(list(rdict[tkey].values()), label=tkey)
                rlabel = rdict[tkey].keys()
                pl.xticks(range(len(rlabel)), rlabel)
        pl.legend()
        pl.show()
        

def retrieveResults(picklefile):
    """Retrieve reduced results from a pickle file."""
    with open(picklefile, 'r') as f:
        results = pickle.load(f)
    return results

def getBestConfig(results, verbose=True):
    """Get best GPU configuration for every timed subroutine from reduced results dictionary."""
    timedSubroutines = results[results.keys()[0]].keys()

    timeGetter = operator.itemgetter(1)
    timeDict   = {}
    for tSub in timedSubroutines:
        timeDict[tSub] = []
        for key in results.iterkeys():
            if results[key][tSub]:
                timeDict[tSub].append((key, results[key][tSub]))
        timeDict[tSub] = sorted(timeDict[tSub], key=timeGetter)
        if verbose:
            if 'time' not in tSub: print('For {} subroutine, the best configuration is {} with {:.5f} s'.format(tSub, timeDict[tSub][0][0], timeDict[tSub][0][1]))
    return timeDict

def getBestFromAll(picklefile, fmrc=32, nmrc=8):
    """Get the best configuration across all the tests ran with different maxregcounts. 'fmrc' corresponds to the first value of maxregcount, and 'nmrc' to the number of different maxregcounts used (with strides of 32)."""
    lm = [fmrc + i*32 for i in xrange(nmrc)]
    if lm[-1] == 256: lm[-1] -= 1
    allDict = {}
    allBest = {}
    for i in lm:
        key = "%d" %i
        allDict[key] = retrieveResults(picklefile.replace("64", key))
        allBest[key] = getBestConfig(allDict[key], verbose=False)

    timeGetter = operator.itemgetter(2)
    timedSubroutines = allDict[key][allDict[key].keys()[0]].keys()
    timeDict = {}
    for tSub in timedSubroutines:
        timeDict[tSub] = []
        for key in allBest.iterkeys():
            timeDict[tSub].append((key, allBest[key][tSub][0][0], allBest[key][tSub][0][1]))
        timeDict[tSub] = sorted(timeDict[tSub], key=timeGetter)
        if 'time' not in tSub: print (bold + 'For {} subroutine:\n'.format(tSub) + reset + 'the best configuration is {} with a maxregcount of {} with {:.5f} s\n'.format(timeDict[tSub][0][0], timeDict[tSub][0][1], timeDict[tSub][0][2]) + '(worst in {:.5f} s)'.format(timeDict[tSub][-1][2]))

    return timeDict

def plotDiff(results, key, labels=None, *args, **kwargs):
    """Plot the elapsed time for a given subroutine (given by 'key') for all the given reduced results. 'results' is a list of reduced results dictionaries. optional 'labels' is a list of labels for the maxregcount of the given results."""
    if not labels:
        labels = [None for i in xrange(len(results))]
    for i, res in enumerate(results):
        rres = revertTimeDict(res)
        rlabel = rres[key].keys()
        pl.title(key)
        pl.plot(list(rres[key].values()), label=labels[i])
        pl.xticks(range(len(rlabel)), rlabel)
    if labels[0]:
        pl.legend(title="maxregcount:")
    pl.show()
