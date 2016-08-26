# -*- coding: utf-8 -*-
"""
Created on Thu May 26 14:25:37 2016

@author: lbignell
"""

from oneTtools.PandasAnalysis import GetData
from optparse import OptionParser
import sys
import psutil
import shelve
import os
import pandas
import numpy as np

class rundata():
    '''
    This is an object cointaining the run data.
    '''
    def __init__(self, filepath, Nevts, runmask):
        print('file: {0}'.format(filepath))
        print('runs: {0}'.format(runmask))
        print('memory used prior: {0} %'.format(psutil.virtual_memory().percent))
        if os.path.exists(filepath):
            print(filepath)
            self.QDCdata, self.TDCdata, self.WFDdata, self.runkeys =\
            GetData(filepath, maxevts=Nevts, runmask=runmask)

        else:
            print('ERROR: {0} does not exist'.format(filepath))

        print('memory used after: {0} %'.format(psutil.virtual_memory().percent))
        
        if os.path.dirname(filepath):
            basename = os.path.dirname(filepath) + os.path.sep

        else:
            basename = ''

        self.shelfpath = basename + 'runs' + \
                        str(self.runkeys[0]) + '--' + str(self.runkeys[-1]) + \
                        '_' + str(Nevts) + 'evts.shelf'
        theshelf = shelve.open(self.shelfpath, writeback=False)
        theshelf['QDCdata'] = self.QDCdata
        theshelf['TDCdata'] = self.TDCdata
        theshelf['WFDdata'] = self.WFDdata
        theshelf['runkeys'] = self.runkeys
        theshelf.close()    

    def gettrigtypes(self):
        '''
        Get event masks for the trigger types.
        '''
        self.LEDtrigs = np.all([pandas.isnull(self.TDCdata['CT']), 
                                pandas.isnull(self.TDCdata['M'])], axis=0)
        
    
        return

if  __name__ == '__main__' :
    parser = OptionParser()
    parser.add_option("-N","--Nevents",default=1E4,type=int,
                      help="Number of events per input file to process. [default %default]")
    parser.add_option("-r","--RunMask",default=None,type=str,
                      help="String to match when selecting Runs. Multiple runs are comma separated, don't include spaces.")
    parser.add_option('-p', '--filepath', default=None, type=str,
                      help='path to the HDF5 file')
    (options, args) = parser.parse_args(args=sys.argv)

    if options.RunMask is not None:
        runs = options.RunMask.split(',')
    else:
        runs = []
    
    if runs:
        runobj = rundata(options.filepath, options.Nevents, runs)
