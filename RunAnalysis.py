# -*- coding: utf-8 -*-
"""
Created on Tue Jun 07 09:51:48 2016

@author: lbignell
"""

import oneTtools.Analyse
from optparse import OptionParser
import sys
import os
import numpy as np

parser = OptionParser()
parser.add_option("-N","--Nevents",default=1E4,type=int,
                  help="Number of events per input file to process. [default %default]")
parser.add_option("-r","--RunRange",default=None,type=str,
                  help="The range of run numbers to match when selecting Runs.\n\
                        Include the first and last run separated by a comma.\n\
                        Don't include spaces.")
parser.add_option('-i', '--inputfilepath', default=None, type=str,
                  help='path to the HDF5 file')
parser.add_option('-o', '--outputfilepath', default=None, type=str,
                  help='path to the output files')
(options, args) = parser.parse_args(args=sys.argv)

if options.RunRange is not None:
    runends = [int(run) for run in options.RunRange.split(',')]
    runs = ['{:06d}'.format(int(run)) 
            for run in np.linspace(runends[0], 
                                   runends[-1], runends[-1] - runends[0] + 1)]
    print('Runs to be processed: {0}'.format(runs))
else:
     print('WARNING: No runs selected!')
     runs = []
print('runs = {0}'.format(runs))

origpath = os.getcwd()
os.chdir(options.outputfilepath)

for run in runs:
    cmd = sys.executable + ' ' + oneTtools.Analyse.__file__ +\
            ' -N {0}'.format(options.Nevents) + ' -r ' + run +\
            ' -p {0}'.format(options.inputfilepath)
    os.system(cmd)
    
os.chdir(origpath)
