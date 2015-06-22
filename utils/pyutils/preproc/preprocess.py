#!/usr/bin/env python
#===============================================================================
## \file preprocess.py
# \brief
# \b Preprocessor for Fortran sources
# This is a simple preprocessor for Fortran sources in Python
# \author
# Marc Joos <marc.joos@cea.fr>
# \copyright
# Copyrights 2014-2015, CEA.
# This file is distributed under the CeCILL-A & GNU/GPL licenses, see
# <http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html> and
# <http://www.gnu.org/licenses/>
# \date
# \b created:          12-04-2014
# \b last \b modified: 06-19-2015

#===============================================================================
import os
import re
import sys
import shutil
import filecmp
import tempfile
from datetime import datetime

# Compute the depth of a list
depth = lambda list_: isinstance(list_, list) and max(map(depth, list_)) + 1

# Flatten a list of lists
def flatten(list_):
    flattened = list_
    while isinstance(flattened[0], list):
        flattened = [item for sublist in flattened for item in sublist]
    return flattened

# Flatten a list of lists to a defined depth
def flattenToDepth(list_, depth):
    flattened = list_
    while depth > 1:
        flattened = [item for sublist in flattened for item in sublist]
        depth -= 1
    return flattened

class File:
    """Class to manage and preprocess files"""
    def __init__(self, fname, overwrite=False, markedForSolver=False \
                     , listSolver=None):
        """Initialization method."""
        self.fname = fname
        sroot = '/'.join(fname.split('/')[:-1])
        sname = fname.split('/')[-1]
        self.oname = sroot + '/o_' + sname
        self.tname = sroot + '/temp_' + sname
        self.overwrite = overwrite
        self.markedForSolver = markedForSolver
        self.regex = re.compile("\s*\!\$py", re.IGNORECASE)
        self.bufferLabel = ''
        self.numTime = [0, 0]
        self.listTime = []
        self.listBegin = []
        self.listSolver = listSolver if markedForSolver else []

    def _findProgram(self, fname):
        """Find the extent of a program.
        Take a filename, return two lists of indices, first list is the beginning of a program, second list is the end of a program."""
        program    = re.compile(r"^\s*program", re.IGNORECASE)
        endProgram = re.compile(r"^\s*end program", re.IGNORECASE)
        progLine    = []
        endProgLine = []
        with open(fname, 'r') as f:
            for i, line in enumerate(f):
                if program.match(line): progLine.append(i)
                if endProgram.match(line): endProgLine.append(i)
        return progLine, endProgLine

    def _findDeclaration(self, list_):
        """Find all variable declarations in a list.
        Take a list, return a list of indices, each index is the line of a declaration."""
        decl = re.compile(r"(^\s*)(type|integer|real|complex|logical|character|namelist)", re.IGNORECASE)
        decLine = []
        decContent = []
        for i, line in enumerate(list_):
            if decl.match(line): 
                decLine.append(i)
                decContent.append(line)
        return decLine, decContent

    def _findDeclarationInFile(self, fname):
        """Find all variable declarations in the source.
        Take a filename, return a list of indices, each index is the line of a declaration."""
        with open(fname, 'r') as f:
            decLine, decContent = self._findDeclaration(f)
        return decLine, decContent

    def _findUsestatement(self, list_):
        """Find all use statements in a list.
        Take a list, return three lists: a list of line indices, a list of module names, and a list of 'only' statements"""
        use = re.compile(r"([\t ]*use[\t ]+)(\w+)((\s*,\s*only\s*:\s*\w+)*(\s*[,]\s*\w+)*)", re.IGNORECASE)
        useLine = []
        module  = []
        only    = []
        for i, line in enumerate(list_):
            if use.match(line):
                useLine.append(i)
                module.append(use.match(line).group(2))
                only.append(use.match(line).group(3))
        return useLine, module, only

    def _findUsestatementInFile(self, fname):
        """Find all use statements in the source.
        Take a filename, return three lists: a list of line indices, a list of module names, and a list of 'only' statements"""
        with open(fname, 'r') as f:
            useLine, module, only = self._findUsestatement(f)
        return useLine, module, only

    def _findImplicitNone(self, list_):
        """Find all 'implicit none' statement in a list.
        Take a list, return a list of line index of the statement."""
        implicitNone = re.compile(r"^\s*implicit none", re.IGNORECASE)
        impLine      = []
        for i, line in enumerate(list_):
            if implicitNone.match(line): impLine.append(i)
        return impLine

    def _findImplicitNoneInFile(self, fname):
        """Find all 'implicit none' statement.
        Take a filename, return a list of line index of the statement."""
        with open(fname, 'r') as f:
            impLine = self._findImplicitNone(f)
        return impLine

    def _findReturn(self, list_):
        """Find 'return' instruction in a list.
        Take a list, return an index line."""
        return_ = re.compile(r"^\s*return", re.IGNORECASE)
        retLine = []
        for i, line in enumerate(list_):
            if return_.match(line): retLine.append(i)
        return retLine

    def _findOMPPragma(self, list_):
        """Find OpenMP pragma in a list.
        Take a lsit, return an list of indices and of OpenMP pragma."""
        OMP = re.compile(r"^\s*\!\$OMP", re.IGNORECASE)
        OMPLine = []
        OMPPragma = []
        for i, line in enumerate(list_):
            if OMP.match(line):
                OMPLine.append(i)
                OMPPragma.append(line)
        return OMPLine, OMPPragma

    def _findSubroutine(self, list_):
        """Find all subroutines in a list.
        Take a list, return two lists of indices, first list is the beginning of a subroutine, second list is the end of a subroutine."""
        subroutine    = re.compile(r"^\s*subroutine", re.IGNORECASE)
        endSubroutine = re.compile(r"^\s*end subroutine", re.IGNORECASE)
        subroutineName = re.compile(r"(?!end)([\t ]*subroutine[\t ]+)(\w+)(\(.*\))?", re.IGNORECASE)
        subLine    = []
        endSubLine = []
        subName    = []
        subArgs    = []
        for i, line in enumerate(list_):
            if subroutine.match(line): subLine.append(i)
            if endSubroutine.match(line): endSubLine.append(i)
            if subroutineName.match(line): 
                subName.append(subroutineName.match(line).group(2))
                subArgs.append(subroutineName.match(line).group(3))
        return [subLine, endSubLine, subName, subArgs]

    def _findSubroutineInFile(self, fname):
        """Find all subroutines in a file.
        Take a filename, return two lists of indices, first list is the beginning of a subroutine, second list is the end of a subroutine."""
        with open(fname, 'r') as f:
            subLine, endSubLine, subName, subArgs = self._findSubroutine(f)
        return [subLine, endSubLine, subName, subArgs]
                    
    def _structure(self, fname):
        """Infer the structure of the file, to learn where to include extra declarations in the preprocessing step.
        Take a filename, return a list of line indices."""
        decList  = self._findDeclarationInFile(fname)
        progList = self._findProgram(fname)
        impList  = self._findImplicitNoneInFile(fname)
        subList  = self._findSubroutineInFile(fname)

        addLineList = []
        listAll = decList[0] + impList + subList[0]
        listAll.sort()
        for i, line in enumerate(listAll):
            if i < len(listAll)-1:
                nextLine = listAll[i+1]
                if line in impList and nextLine in subList[0]:
                    addLineList.append(line)
                if line in decList[0] and nextLine in subList[0]:
                    addLineList.append(line)
                if line in subList[0] and nextLine in subList[0]:
                    addLineList.append(line)
            else:
                addLineList.append(line)
        return addLineList

    def _extractExcerpt(self, fname, ibegin, iend):
        """Extract an excerpt from a file.
        Take a filename, and two line indices (begin and end of the region to retrieve in the source file)."""
        extLine = []
        with open(fname, 'r') as f:
            for i, line in enumerate(f):
                if i >= ibegin and i <= iend:
                    extLine.append(line)
                elif i > iend:
                    break
        return extLine

    def _addStatement(self, addLineOrig, subListOrig, listElmt, variable \
                          , module):
        """Add a 'use module' statement if needed.
        Take a list retrieve from self._structure of the original file, a list retrieve from self._findSubroutine of the original file, a list of elements to add, and two strings (the name of the variable to add, and the name of the module to add)."""
        addLineNew = self._structure(self.tname)
        subListNew = self._findSubroutineInFile(self.tname)
        useList    = self._findUsestatementInFile(self.tname)
        decList    = self._findDeclarationInFile(self.tname)
        impList    = self._findImplicitNoneInFile(self.tname)
        addUse = [[], [], []]
        if listElmt:
            for i in listElmt:
                distance = []
                for j in addLineOrig:
                    distance.append(i - j)
                mindist = max(distance)
                for dist in distance:
                    if dist > 0:
                        mindist = min(mindist, dist)
                addLine = addLineOrig[distance.index(mindist)]
                # Find current subroutine
                for j in xrange(len(subListOrig[0])):
                    if addLine > subListOrig[0][j] and \
                            addLine < subListOrig[1][j]:
                        beginSub = subListNew[0][j]
                        endSub   = subListNew[1][j]
                        nameSub  = subListNew[2][j]
                        break
                # Find current use statements, to check if 'params' is defined
                listUse = [[],[],[]]
                for j in xrange(len(useList[0])):
                    if useList[0][j] > beginSub and useList[0][j] < endSub:
                        listUse[0].append(useList[0][j])
                        listUse[1].append(useList[1][j])
                        listUse[2].append(useList[2][j])
                # Find current declarations, to see if 'verbose' is defined
                listDec = [[],[]]
                for j in xrange(len(decList[0])):
                    if decList[0][j] > beginSub and decList[0][j] < endSub:
                        listDec[0].append(decList[0][j])
                        listDec[1].append(decList[1][j])
                # Find 'implicit none' statement, to determine where to include
                # the use statement
                implicitPos = []
                for j in xrange(len(impList)):
                    if impList[j] > beginSub and impList[j] < endSub:
                        implicitPos.append(impList[j])
                        break
                # Finally, check where to include the use statement
                # 0: include 'use params, only: verbose'
                # 1: do nothing (but 'use params' is already there)
                # 2: include ', verbose'
                if listDec[0]:
                    allDec = " ".join(listDec[1])
                    if not variable in allDec:
                        if not module in listUse[1]:
                            if implicitPos:
                                addUse[0].append(implicitPos-1)
                            else:
                                addUse[0].append(beginSub)
                        else:
                            indPar = listUse[1].index(module)
                            if listUse[2][indPar]:
                                if not variable in listUse[2][indPar]:
                                    addUse[2].append(listUse[0][indPar])
                            else:
                                addUse[1].append(listUse[0][indPar])
        
        # Write the file
        with open(self.tname, 'r') as t:
            lines = list(t)
        if addUse[0]:
            for i, addline in enumerate(addUse[0]):
                lines[addline] += "  use " + module + ", only: " + \
                    variable + "\n"
        if addUse[2]:
            for i, addline in enumerate(addUse[2]):
                lines[addline] = lines[addline][:-1] + ", " + variable + "\n"
        with open(self.tname, 'w') as t:
            t.writelines(lines)

    def _addDeclaration(self, addLineOrig, listElmt, commentLine, decLine):
        """Add declaration if needed.
        Take a list retrieve from self._structure of the original file, a list of line of elements to add, and two strings (a comment line to explicitly state the declaration, and the actual declaration in the form ' integer :: a')."""
        addLineNew = self._structure(self.tname)
        addDecl = []
        if listElmt:
            for i in listElmt:
                distance = []
                for j in addLineOrig:
                    distance.append(i - j)
                mindist = max(distance)
                for dist in distance:
                    if dist > 0:
                        mindist = min(mindist, dist)
                addDecl.append(addLineNew[distance.index(mindist)])

        # Write the file
        addDeclDict = dict((i, addDecl.count(i)) for i in addDecl)
        addDeclClean = addDeclDict.keys()
        with open(self.tname, 'r') as t:
            lines = list(t)
        for addline in addDeclClean:
            if commentLine:
                lines[addline] += commentLine
            lines[addline] += decLine
        with open(self.tname, 'w') as t:
            t.writelines(lines)

    def _parseLine(self, line):
        """Parse line to find if it contains a !$py pragma."""
        search = self.regex.match(line)
        ppline = []
        if search:
            ppline = search.string.strip().lower()
            ppline = ppline.split('!$py', 1)[1].strip()
        return ppline

    def _preprocess_start_timing(self, ppline, i, line):
        """Preprocess a 'start_timing' instruction."""
        self.bufferLabel = " ".join(ppline.split()[1:])
        line = "  if (verbose) call system_clock(count=t0, count_rate=irate)\n"
        self.numTime[0] += 1
        self.listTime.append(i)
        return line

    def _preprocess_end_timing(self, ppline, line):
        """Preprocess a 'end_timing' instruction."""
        currentLabel = " ".join(ppline.split()[1:])
        try:
            assert self.bufferLabel == currentLabel
        except AssertionError:
            print("In file: %s" %self.fname)
            print("Timing error: labels do not match: %s != %s" \
                      %(self.bufferLabel, currentLabel))
            sys.exit(0)
        line  = "  if (verbose) then\n"
        line += "    call system_clock(count=t1, count_rate=irate)\n"
        line += "    print '(\"" + currentLabel \
              + ":   \", F12.8, \" s\")', (t1 - t0)/(irate*1.d0)\n"
        line += "  endif\n"
        self.numTime[1] += 1
        return line

    def _preprocess_begin_statement(self, i, line, label):
        """Preprocess a 'begin_statement' instruction."""
        line = "  if (verbose) print*, '> Entering " + label + "'\n"
        self.listBegin.append(i)
        return line

    def _preprocess_call_solver(self, fname):
        """Preprocess a 'call_solver' instruction."""
        # Retrieve the list of solvers generated previously
        depthSolver   = depth(self.listSolver)
        flattenSolver = flattenToDepth(self.listSolver, depthSolver-1)
        solverName = []
        solverArgs = []
        for solver in flattenSolver:
            solverName.append(solver[0])
            solverArgs.append(solver[1])

        # Retrieve the position where to insert solvers
        callSolverLine = []
        callSolverName = []
        callSolverArgs = []
        with open(fname, 'r') as f:
            for i, line in enumerate(f):
                ppline = self._parseLine(line)
                if ppline:
                    if ppline.startswith("call_solver"):
                        currentSolverName = ppline.split()[1]
                        currentSolverArgs = ' '.join(ppline.split()[2:])
                        callSolverLine.append(i)
                        callSolverName.append(currentSolverName)
                        callSolverArgs.append(currentSolverArgs)
        listNameSolver   = list(set(callSolverName))
        listLengthSolver = map(len, listNameSolver)
        listNameSolver   = [item for (len_, item) in \
                   sorted(zip(listLengthSolver, listNameSolver), reverse=True)]

        # Create a dictionary for each type of solver
        dictSolver = {}
        for name in listNameSolver:
            dictSolver[name] = []

        # Fill with the name and arguments of each associated solver
        for name in listNameSolver:
            remove = False
            for i, currentName in enumerate(solverName):
                if name in currentName:
                    dictSolver[name].append([solverName[i], solverArgs[i]])
                    remove = True
            if remove:
                for item in dictSolver[name]:
                    solverName.remove(item[0])
                    solverArgs.remove(item[1])
                    
        # Generate the lines to insert into the file
        newLines = []  
        for i, name in enumerate(callSolverName):
            try:
                assert(name in ['riemann', 'riemann_magnetic'])
            except AssertionError:
                print("Problem in generating call to solver:")
                print("Unknown solver type.")
                sys.exit(0)

            dimension = '' if name == 'riemann' else '2d'
            solverSet = dictSolver[name]
            line = ''
            for currentSolver in solverSet:
                try:
                    assert(len(callSolverArgs[i].split(',')) == len(currentSolver[1].split(',')))
                except AssertionError:
                    print("Problem in generating call to solver:")
                    print("The number of arguments in the call and in the solver does not fit.")
                    sys.exit(0)
                
                solverType = currentSolver[0].replace(name + '_', '')

                line += "  if (iriemann" + dimension + " == i" + solverType 
                line += ") then\n"
                line += "     call " + currentSolver[0] + callSolverArgs[i] 
                line += "\n  endif\n"
            newLines.append(line)

        return callSolverLine, newLines

    def _preprocess_solver(self, fname):
        """Preprocess a solver file.
        Take a filename, return a list of lines of code."""
        solverGenericTemplate = []
        insertSolver = []
        solverTemplate = []
        # Retrieve solver template definition in the source file
        with open(fname, 'r') as f:
            for i, line in enumerate(f):
                ppline = self._parseLine(line)
                if ppline:
                    if ppline.startswith("global_template"): 
                        solverGenericTemplate.append([i, line])
                    if ppline.startswith("insert_solver"):
                        insertSolver.append(i)
                    if ppline.startswith("solver_template"):
                        solverTemplate.append([i, line])

        # Perform a serie of check to assure the coherence of the template
        try:
            assert(len(solverGenericTemplate)%2 == 0)
        except AssertionError:
            print("In file: %s" %fname)
            print("Generic template solver is badly defined.")
            sys.exit(0)

        try:
            assert(len(solverGenericTemplate) == 2)
        except AssertionError:
            print("In file: %s" %fname)
            print("Generic template solver is badly defined: only one general template can be defined by file.")
            sys.exit(0)

        try:
            assert(len(solverTemplate)%2 == 0)
        except AssertionError:
            print("In file: %s" %fname)
            print("Template solvers are badly defined.")
            sys.exit(0)

        for i in xrange(0,len(solverTemplate),2):
            try:
                assert(solverTemplate[i][-1].split()[-1] == solverTemplate[i+1][-1].split()[-1])
            except AssertionError:
                print("In file: %s" %fname)
                print("Template solver " + solverTemplate[i][-1].split()[-1] + " is badly defined.")
                sys.exit(0)

        try:
            assert(len(insertSolver) == len(solverGenericTemplate)/2)
        except AssertionError:
            print("In file: %s" %fname)
            print("Generic template solver is badly defined: the number of 'insert_solver' instructions does not fit the number of generic template.")
            sys.exit(0)

        for i in xrange(0,len(solverGenericTemplate),2):
            begin = solverGenericTemplate[i][0]
            end   = solverGenericTemplate[i+1][0]
            try:
                assert(begin < insertSolver[i/2] & insertSolver[i/2] < end)
            except AssertionError:
                print("In file: %s" %fname)
                print("Generic template solver is badly defined: it does not contain an 'insert_solver' instruction.")
                sys.exit(0)
            
        # Build the solvers:
        #  - start by recovering informations on the global template
        #  - retrieve OpenMP instructions to update pragma for each solver
        #  - loop over all the solvers to generate the given solvers
        fullSolver = []
        # labelGeneric   = solverGenericTemplate[0][-1].split()[-1]
        beginGeneric   = solverGenericTemplate[0][0] + 1
        endGeneric     = solverGenericTemplate[1][0] - 1
        genericSolver  = self._extractExcerpt(fname, beginGeneric, endGeneric)
        gsUsestatement = self._findUsestatement(genericSolver)
        gsDeclaration  = self._findDeclaration(genericSolver)
        gsSubroutine   = self._findSubroutine(genericSolver)
        gsOMPPragma    = self._findOMPPragma(genericSolver)
        labelGeneric   = gsSubroutine[2][0]
        argsGeneric    = gsSubroutine[3][0]
        insertSolver   = insertSolver[0] - beginGeneric

        # Retrieve the position of the 'private' OpenMP pragma prior to the 
        # solver insertion
        OMPPrivate = gsOMPPragma[0][0]
        i = 0
        while (OMPPrivate < insertSolver):
            OMPPrivate = gsOMPPragma[0][i+1]
            i += 1
        lastPragma = gsOMPPragma[0][i-1]
        i -= 1
        private = re.compile(r'(^\s*)(\!\$OMP.*PRIVATE)', re.IGNORECASE)
        while (lastPragma > 0):
            if private.match(gsOMPPragma[1][i]):
                blankSpace = private.match(gsOMPPragma[1][i]).group(1)
                break
            else:
                i -= 1
                lastPragma = gsOMPPragma[0][i]
        privatePos = gsOMPPragma[0][i]

        for i in xrange(0,len(solverTemplate),2):
            # Retrieve solver structure
            beginSolver   = solverTemplate[i][0] + 1
            endSolver     = solverTemplate[i+1][0] - 1
            currentSolver = self._extractExcerpt(fname, beginSolver, endSolver)
            sUsestatement = self._findUsestatement(currentSolver)
            sDeclaration  = self._findDeclaration(currentSolver)
            sReturn       = self._findReturn(currentSolver)
            sSubroutine   = self._findSubroutine(currentSolver)

            labelSolver = sSubroutine[2][0]
            labelFull   = labelGeneric + '_' + labelSolver
            self.listSolver.append([labelFull, argsGeneric])

            beginBody = sDeclaration[0][-1] + 1
            endBody   = sSubroutine[1][0]-1 if not sReturn else sReturn[0]-1

            # Create new solver and change its name
            newSolver = list(genericSolver)
            newSolver[gsSubroutine[0][0]] = newSolver[gsSubroutine[0][0]].replace(labelGeneric, labelFull)
            newSolver[gsSubroutine[1][0]] = newSolver[gsSubroutine[1][0]].replace(labelGeneric, labelFull)

            # Check if new 'use' statement is needed
            for i in xrange(len(sUsestatement[0])):
                line     = sUsestatement[0][i] 
                module   = sUsestatement[1][i]
                variable = sUsestatement[2][i]
                if module not in gsUsestatement[1]:
                    newSolver[gsUsestatement[0][-1]] = newSolver[gsUsestatement[0][-1]] + '! Automatically generated by preprocessor >\n  use ' + module + variable + '\n! < Automatically generated by preprocessor\n'

            # Add variable declarations
            newSolver[gsDeclaration[0][-1]] = newSolver[gsDeclaration[0][-1]] +  '! Automatically generated by preprocessor >\n'
            newSolver[gsDeclaration[0][-1]] += ''.join(sDeclaration[1])
            newSolver[gsDeclaration[0][-1]] = newSolver[gsDeclaration[0][-1]] + '! < Automatically generated by preprocessor\n'

            # Retrieve internal variables of the solver; they need to be private
            # in OpenMP
            privateVar = []
            for lineOfDec in sDeclaration[1]:
                privateVar.append(lineOfDec.split('::')[-1][:-1].split(','))
            privateVar = flatten(privateVar)

            # Update the OpenMP 'DO PRIVATE' pragma
            if privateVar:
                blank = re.compile(r"(^\s*)")
                blankSpace = len(blank.match(newSolver[privatePos]).group(1))
                blankSpace = " "*blankSpace
                npvar = len(privateVar)
                nvarBatch = npvar/10
                nvarRemains = npvar%10
                i = 0
                newSolver[privatePos] = '! Automatically modified by preprocessor >\n' + newSolver[privatePos]
                while i < nvarBatch:
                    newSolver[privatePos] = newSolver[privatePos][:-1] + " &\n"
                    newSolver[privatePos] += blankSpace + "!$OMP PRIVATE(" \
                        + ",".join(privateVar[i*10:(i+1)*10]) + ")\n"
                    i += 1
                if nvarRemains: 
                    newSolver[privatePos] = newSolver[privatePos][:-1] \
                       + " &\n" + blankSpace \
                       + "!$OMP PRIVATE(" + ",".join(privateVar[i*10:]) + ")\n"
                newSolver[privatePos] = newSolver[privatePos] + "! < Automatically modified by preprocessor\n"

            # Add the main body of the solver
            blank = re.compile(r"(^\s*)")
            blankSpace = len(blank.match(newSolver[insertSolver]).group(1))
            blankCurrent = len(blank.match(currentSolver[beginBody]).group(1))
            blankSpace = " "*(blankSpace - blankCurrent - 1)
            newSolver[insertSolver] = '! Automatically generated by preprocessor >\n' + blankSpace.join(currentSolver[beginBody:endBody+1])
            newSolver[insertSolver] = re.sub(r"\n\s*#", r"\n#", newSolver[insertSolver])
            newSolver[insertSolver] = newSolver[insertSolver] + '! < Automatically generated by preprocessor\n'

            fullSolver.append(newSolver)

        allSolver = flatten(fullSolver)

        # Finally, add a warning at the beginning
        allSolver[0] = "! This source file was automatically generated by DUMSES preprocessor; do not modify it, unless you know what you are doing!\n" + allSolver[0]
        return allSolver

    def _preprocess(self, i, line, label):
        """Preprocess a given line.
        For now, the only preprocessing instructions used are 'start_timing' and 'end_timing' to time portions of the code, and 'begin_statement' to print (in verbose mode) when it enters a given subroutine."""
        ppline = self._parseLine(line)
        if ppline:
            if ppline.startswith("start_timing"): 
                line = self._preprocess_start_timing(ppline, i, line)
            if ppline.startswith("end_timing"):
                line = self._preprocess_end_timing(ppline, line)
            if ppline.startswith("begin_statement"):
                line = self._preprocess_begin_statement(i, line, label)
        return line

    def preprocess(self):
        """Preprocess a file."""
        # First, check if there are solvers to generate
        isSolver = False
        with open(self.fname, 'r') as f:
            for i, line in enumerate(f):
                ppline = self._parseLine(line)
                if 'solver_template' in ppline:
                    isSolver = True
                    break

        if isSolver:
            allSolver = self._preprocess_solver(self.fname)
            tmpfile    = self.tname
            self.tname = self.fname
            with open(self.tname, 'w') as t:
                t.writelines(allSolver)
            self.tname = tmpfile

        # Second, read and preprocess the file
        addLineOrig = self._structure(self.fname)
        subList = self._findSubroutineInFile(self.fname)
        try:
            assert(len(subList[0]) == len(subList[1]))
        except AssertionError:
            print("Error: something went wrong processing " + self.fname \
                      + "; the subroutine definitions cannot be understood.")
            print("You probably want to check your 'subroutine/end subroutine' declarations.")
            sys.exit(42)

        with open(self.fname, 'r') as f:
            with open(self.tname, 'w') as t:
                isub = 0
                for i, line in enumerate(f):
                    # Find the name of the current subroutine
                    if subList[0]:
                        if subList[1][isub] > i:
                            label = subList[2][isub]
                        else:
                            if isub < len(subList[0]) - 1: isub += 1
                            label = subList[2][isub]
                    else:
                        label = ''
                    line = self._preprocess(i, line, label)
                    t.write(line)

        # Check if the preprocessing directives were OK
        try:
            assert self.numTime[0] == self.numTime[1]
        except AssertionError:
            problem = "'start_timing'" \
                if (self.numTime[1] > self.numTime[0]) else "'end_timing'"
            print("Something went wrong: you forgot one or more " + problem)
            sys.exit(42)
                                
        # Add use statement if needed
        self._addStatement(addLineOrig, subList, self.listBegin, "verbose" \
                               , "params")

        # Declaration inclusion, if needed
        comment = "\n  ! Timing variables\n"
        decl    = "  integer :: t0, t1, irate\n"
        self._addDeclaration(addLineOrig, self.listTime, comment, decl)

        # Lastly, check if there are call to solvers to generate
        # This step has to be done in a second pass
        if not self.markedForSolver:
            with open(self.fname, 'r') as f:
                for i, line in enumerate(f):
                    ppline = self._parseLine(line)
                    if 'call_solver' in ppline:
                        self.markedForSolver = True
                        break
        else: 
            solverLine, newLines = self._preprocess_call_solver(self.fname)
            lineIndex = 0
            with open(self.fname, 'r') as f:
                with open(self.tname, 'w') as t:
                    for i, line in enumerate(f):
                        if lineIndex < len(solverLine):
                            if i == solverLine[lineIndex]:
                                line = newLines[lineIndex]
                                lineIndex += 1
                        t.write(line)
            self.markedForSolver = False

        # Overwrite the original file
        if self.overwrite:
            if os.path.isfile(self.oname):
                if not filecmp.cmp(self.tname, self.oname):
                    if not self.markedForSolver:
                        os.remove(self.oname)
                    shutil.move(self.tname, self.fname)
                else:
                    os.remove(self.tname)
                    shutil.move(self.oname, self.fname)
            else:
                shutil.move(self.tname, self.fname)

        return self.markedForSolver, self.listSolver

class FileTree:
    """Class to manage directory tree"""
    def __init__(self, inputDir="./", tmpRoot='./tmp'): #tempfile.mkdtemp()):
        """Initialization method."""
        self.dir = inputDir
        self.tmpRoot = tmpRoot
        if self.tmpRoot:
            self.tmp = True
        else:
            self.tmp = False

    # def __del__(self):
    #     """Deletion method."""
    #     shutil.rmtree(self.tmpRoot, ignore_errors=True)

    def listFiles(self):
        """List directory and files."""
        self.dirList = []
        self.fileList = []
        self._listIgnore = []
        dot = re.compile(r'\.[^\/]')
        for root, dirs, files in os.walk(self.dir):
            isdot = [dot.match(root.split('/')[-n]) for n in xrange(len(root.split('/')))]
            if not any(isdot):
                self._listFiles(root, files)
        self._cleanFiles()

    def _listFiles(self, root, files):
        """List files and files to ignore."""
        root = (root[:-1] if root[-1] == '/' else root)
        if root.split('/')[-1][0] != '.' or len(root.split('/')[-1]) == 1:
            if self._listIgnore:
                for ignore in self._listIgnore:
                    if root[:len(ignore)] == ignore:
                        return
                    else:
                        if files:
                            self.dirList.append(root)
                            self.fileList.append(files)
            else:
                if files:
                    self.dirList.append(root)
                    self.fileList.append(files)
        else:
            self._listIgnore.append(root)

    def _cleanFiles(self):
        """Clean directories and files to ignore."""
        if self.dirList:
            for i, fileList in enumerate(self.fileList):
                iterList = fileList.__iter__()
                filesToRemove = []
                for file_ in iterList:
                    if not(file_[-4:] in ['.f90', '.F90'] or \
                       file_[-2:] in ['.f', '.F']):
                        filesToRemove.append(file_)
                    if file_[0] == '#' or file_[:2] == '.#':
                        print('Warning: you probably forgot to save ' + file_ \
                                  + ' in ' + self.dirList[i])
                        if file_ not in filesToRemove:
                            filesToRemove.append(file_)
                if filesToRemove:
                    for file_ in filesToRemove:
                        self.fileList[i].remove(file_)
            for i, fileList in enumerate(self.fileList):
                if not fileList: 
                    self.fileList.remove(self.fileList[i])
                    self.dirList.remove(self.dirList[i])

    def _removeTrailingDot(self, dirList):
        """Remove trailing dots in a given path."""
        cleanDirs = []
        trailingDot = re.compile(r'\.\.\/')
        for dir_ in dirList:
            cleanDirs.append(trailingDot.sub('', dir_))
        return cleanDirs

    def createTmpTree(self):
        """Create temporary tree containing all the Fortran sources."""
        self.tmpTree = self._removeTrailingDot(self.dirList)
        for dir_ in self.tmpTree:
            try:
                dir_ = (dir_[1:] if dir_[0] == "." else dir_)
                os.makedirs(self.tmpRoot + dir_)
            except OSError:
                pass
        for i, dir_ in enumerate(self.tmpTree):
            dir_ = self.tmpRoot + (dir_[1:] if dir_[0] == "." else dir_)
            self.tmpTree[i] = dir_
            for file_ in self.fileList[i]:
                if os.path.isfile(dir_ + '/' + file_):
                    shutil.move(dir_ +  '/' + file_, dir_ + '/' + 'o_' + file_)
                shutil.copy(self.dirList[i] + '/' + file_, dir_ + '/' + file_)

    def processAllFiles(self):
        """Preprocess all the files in the temporary tree."""
        self.listFiles()
        self.createTmpTree()
        listSolver = []
        listMarked = []
        tree = self.tmpTree
        for i, dir_ in enumerate(tree):
            for file_ in self.fileList[i]:
                ftp = File(dir_ + '/' + file_, overwrite=True)
                currentMarked, currentSolver = ftp.preprocess()
                if currentSolver: listSolver.append(currentSolver)
                if currentMarked: listMarked.append(dir_ + '/' + file_)
        # Second pass to generate calls to solver
        if listMarked:
            for file_ in listMarked:
                ftp = File(file_, overwrite=True, markedForSolver=True \
                                , listSolver=listSolver)
                currentMarked, currentSolver = ftp.preprocess()
        with open(self.tmpRoot + '/log.dat', 'w') as f:
            gStr ="Generated on: " + datetime.now().strftime("%Y-%m-%d %H:%M")
            f.write(gStr)
