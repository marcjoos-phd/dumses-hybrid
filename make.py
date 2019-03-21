#!/usr/bin/env python
#===============================================================================
# DUMSES-Hybrid:
# DUMSES make script.
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
# created:       12-15-2014
# last modified: 06-25-2015
#===============================================================================
import os, sys, re, shutil
sys.path.append('utils/pyutils/')
from preproc import FileTree
if sys.version_info[0] <= 2:
    from commands import getstatusoutput as cmd
else:
    from subprocess import getstatusoutput as cmd

def editMakefile(fname='Makefile'):
    tname  = fname + '.temp'
    tmpdir = re.compile('^TMPDIR')
    with open(fname, 'r') as f, open(tname, 'w') as t:
        for line in f:
            if tmpdir.match(line):
                line = "TMPDIR  = tmp\n"
            t.write(line)
    shutil.move(tname, fname)

def main():
    try:
        os.makedirs('./bin')
    except OSError:
        pass
    tree = FileTree('./src')
    tree.processAllFiles()
    editMakefile()
    st, out = cmd('make')
    print(out)
    
if __name__ == "__main__":
    main()
