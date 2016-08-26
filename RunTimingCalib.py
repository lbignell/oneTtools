# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 09:47:33 2016

@author: lbignell

This code requires that the output HDF5 file already be processed into shelf
files containing the pandas dataframes, one for each run.

"""

import oneTtools.makeplots_save
from argparse import ArgumentParser
import sys
import os
import json
import subprocess

parser = ArgumentParser()
parser.add_argument("-N","--Nevents",default=1E4,type=int,
                  help="Number of events per input file to process. [default %default]")
parser.add_argument("-r","--Runs",default=None,type=str,
                  help="A semicolon-separated string of comma-separated run \
                  numbers.\nDon't use spaces. The colons separate LED \
                  reference positions, and the commas separate runs at that \
                  position.")
parser.add_argument('-c', '--refchans', default=None, type=str,
                  help="A comma-separated string of the reference channels. Don't use spaces")
parser.add_argument('-d', '--datafolder', default=None, type=str,
                  help='path to the data folder')
parser.add_argument('-o', '--outputfolder', default=None, type=str,
                  help='folder in which to put the output files')
parser.add_argument('-t', '--timecutlow', default=0, type=int,
                  help="The early cut on absolute pulse arrival time (ns) after trig.")
parser.add_argument('-T', '--timecuthigh', default=2560, type=int,
                  help="The late cut on absolute pulse arrival time (ns) after trig.")
parser.add_argument('-a', '--areacutlow', default=-100000, type=int,
                  help="The low cut on pulse area (AU).")
parser.add_argument('-A', '--areacuthigh', default=0, type=int,
                  help="The high cut on pulse area (AU).")
parser.add_argument('-s', '--reltimecutlow', default=-200, type=int,
                  help="The early cut on pulse arrival time (ns) relative to ref chan.")
parser.add_argument('-S', '--reltimecuthigh', default=200, type=int,
                  help="The late cut on pulse arrival time (ns) relative to ref chan.")
#(options, args) = parser.parse_args(args=sys.argv)
options = parser.parse_args()

error = False

if options.Runs is not None:
    runlist = [runs for runs in options.Runs.split(';')]
    print('Runs to be processed: {0}'.format(runlist))
else:
     print('ERROR: No runs selected!')
     error = True
     
if options.refchans is not None:
    refchans = [ch for ch in options.refchans.split(',')]
    print('Ref chans: {0}'.format(refchans))
else:
    print('ERROR: No reference channels selected')
    error = True
    
if len(runlist) != len(refchans):
    print("ERROR: the length of the ref chan list ({0}) doesn't match the \
            length of the run list ({1})".format(len(refchans), len(runlist)))
    error = True

if not os.path.exists(options.datafolder):
    print("ERROR: the path for the data folder ({0}) doesn't exist!".format(
            options.datafolder))
    error = True

tl = options.timecutlow
th = options.timecuthigh
al = options.areacutlow
ah = options.areacuthigh
rtl = options.reltimecutlow
rth = options.reltimecuthigh

cuts = {"timecutlow": tl, "timecuthigh": th, \
        "areacutlow": al, "areacuthigh": ah, \
        "reltimecutlow": rtl, "reltimecuthigh": rth}
print("cuts = {0}".format(cuts))
strcuts = r'''{{"timecutlow":{0}, "timecuthigh":{1}, "areacutlow":{2}, "areacuthigh":{3}, "reltimecutlow":{4}, "reltimecuthigh":{5}}}'''\
        .format(tl, th, 
                al, ah,
                rtl, rth)
#strcuts = json.dumps(cuts)
print("string cuts = {0}".format(strcuts))

#make the output folder if it doesn't exist
os.makedirs(options.outputfolder, exist_ok=True)    

if not error:
    for i,runs in enumerate(runlist):
        #opts = oneTtools.makeplots_save.__file__ +\        
        cmd = sys.executable + ' ' + oneTtools.makeplots_save.__file__ +\
                 ' -N {0}'.format(options.Nevents) + ' -r ' + runs +\
                ' -c ' + refchans[i] + ' -l ' + refchans[i] +\
                ' -d ' + options.datafolder + ' -o ' + options.outputfolder +\
                " --cutdict " + json.dumps(strcuts)
        print('command: {0}'.format(cmd))
        #print('options: {0}'.format(opts))
        os.system(cmd)
        #subprocess.check_output(['python3', opts])
