# -*- coding: utf-8 -*-
"""
Created on Wed May 04 15:46:23 2016

@author: lbignell
"""
import h5py
import pandas
import matplotlib.pyplot as plt
import numpy as np
import gc
import sys

def myenumerate(x):
    if type(x) is np.float64 or type(x) is int:
        return enumerate([x])
    elif np.ndim(x)>1 and type(x) is np.ndarray:
        return enumerate(x[0])
    else:
        return enumerate(x)

def get_first(iterable, default=None):
    if iterable:
        for item in iterable:
            return item
    return default
        
def loopruns(rungroup, runs, maxevts=10000):
    emptyWFDevt = {'area': np.nan, 'npulse': np.nan, 'nsubp': np.nan,
                   'ped': np.nan, 'pedsd': np.nan, 'time': np.nan}
    QDCrunlist = []
    TDCrunlist = []
    WFDrunlist = []
    runkeys = []
    for runkey in runs:
        print('Run = {0}'.format(runkey))
        thisrun = rungroup[runkey]
        rnum = int(runkey)
        evtgroup = thisrun['Event']
        QDCevtlist = []
        TDCevtlist = []
        WFDevtlist = []
        QDCdict = {}
        TDCdict = {}
        WFDdict = {}
        evtcount = 0
        for evtkey in evtgroup.keys():
            if evtcount > maxevts:
                break
            evtcount += 1
            #if int(evtkey)%100==0: print('key = {0}'.format(evtkey))
            QDC = evtgroup[evtkey]['QDC']
            TDC = evtgroup[evtkey]['TDC']
            hazWFD = True
            try:
                WFD = evtgroup[evtkey]['WFD']
            except KeyError:
                hazWFD = False
            if sys.version_info[0] == 3:
                QDCdict = {}
            else:
                QDCalldata = {QDC[ch].name.split('/')[-1] : QDC[ch].value for ch in QDC.keys()}
                QDCdict = {key : QDCalldata[key]['QDC Data'] for key in QDCalldata.keys()}
            TDCdict = {TDC[ch].name.split('/')[-1] : [TDC[ch].value[1]] for ch in TDC.keys()}
            WFDdict = {}
            if hazWFD:
                WFDfields = WFD.keys()
                WFDchs = WFD[get_first(WFDfields)].keys()
                WFDdict = {ch : 
                           {field : 
                            pandas.DataFrame.from_dict({ID : 
                            [val] for (ID,val) in myenumerate(WFD[field][ch].value)}, 
                            orient='index').rename(columns={0:'entries'}) 
                             for field in WFDfields} 
                             for ch in WFDchs}
                WFDevtlist += [pandas.DataFrame(WFDdict)]
            else:
                WFDevtlist += [pandas.DataFrame({'NoPulses': emptyWFDevt})]
            QDCevtlist += [pandas.DataFrame(QDCdict)]
            TDCevtlist += [pandas.DataFrame(TDCdict)]
            #del QDC, TDC, WFD, QDCalldata, QDCdict, TDCdict, WFDdict
            #del QDCdict, TDCdict, WFDdict
            #gc.collect()
        try:
            QDCrun = pandas.concat(QDCevtlist, keys=range(len(QDCevtlist)))
            TDCrun = pandas.concat(TDCevtlist, keys=range(len(TDCevtlist)))
            WFDrun = pandas.concat(WFDevtlist, keys=range(len(WFDevtlist)))
            QDCrunlist += [QDCrun]
            TDCrunlist += [TDCrun]
            WFDrunlist += [WFDrun]
            runkeys += [rnum]
        except MemoryError:
            break
    return QDCrunlist, TDCrunlist, WFDrunlist, runkeys

def GetData(fname, maxevts=10000, runmask=[]):
    '''
    This code prepares the data from David's HDF5 output format for Pandas analysis.

    QDC is a data frame indexed by [run #][evt #][channel name]
    TDC is a data frame indexed by [run #][channel name]
    WFD is a data frame indexed by [run #][evt #][pulse parameter][channel name]
    '''
    h5file = h5py.File(fname, 'r')
    rungroup = h5file['Run']
    QDCrunlist = []
    TDCrunlist = []
    WFDrunlist = []
    runkeys = []
    if runmask:
        runs = [run for run in rungroup.keys() if run in runmask]
    else:
        runs = rungroup.keys()
    QDCrunlist, TDCrunlist, WFDrunlist, runkeys = loopruns(rungroup, runs, 
                                                           maxevts=maxevts)
    if len(runs) > 1:
        QDCdata = pandas.concat(QDCrunlist, keys=runkeys)
        TDCdata = pandas.concat(TDCrunlist, keys=runkeys)
        WFDdata = pandas.concat(WFDrunlist, keys=runkeys)
    else:
        QDCdata = QDCrunlist
        TDCdata = TDCrunlist
        WFDdata = WFDrunlist
    return QDCdata, TDCdata, WFDdata, runkeys

def processWFDtimediff(fname, maxevts=10000, **kwargs):
    h5file = h5py.File(fname)
    rungroup = h5file['Run']
    histdir_LED = {'S{0} - S2'.format(i) : [] for i in range(6)}
    histdir_CT = {'S{0} - S2'.format(i) : [] for i in range(6)}
    histdir_M = {'S{0} - S2'.format(i) : [] for i in range(6)}
    for run in rungroup.keys():
        QDCdata, TDCdata, WFDdata = GetData(fname, runmask=run, maxevts=maxevts)
        #WFD_S0mS2 = WFDdata.S0 - WFDdata.S2
        #WFD_S1mS2 = WFDdata.S1 - WFDdata.S2
        #WFD_S3mS2 = WFDdata.S3 - WFDdata.S2
        #WFD_S4mS2 = WFDdata.S4 - WFDdata.S2
        #WFD_S5mS2 = WFDdata.S5 - WFDdata.S2
        diffdict = {'S{0} - S2'.format(i) : 
                    getattr(WFDdata, 'S{0}'.format(i)) - getattr(WFDdata, 'S2')
                    for i in range(6)}
        LEDevts = getLEDevts(TDCdata)
        CTevts = getCTevts(TDCdata)
        Mevts = getMevts(TDCdata)
        LEDevts_WFD = WFDLEDevts(LEDevts)
        CTevts_WFD = WFDLEDevts(CTevts)
        Mevts_WFD = WFDLEDevts(Mevts)
        for i in range(6):
            thekey = 'S{0} - S2'.format(i)
            LEDhist = plotreltimeWFD(diffdict[thekey], LEDevts_WFD,
                                     thekey+', LED')
            CThist = plotreltimeWFD(diffdict[thekey], CTevts_WFD,
                                    thekey+', CT')
            Mhist = plotreltimeWFD(diffdict[thekey], Mevts_WFD,
                                   thekey+', M')
            if i is 0:
                xbins = LEDhist[1][0:-1]
            histdir_LED[thekey].append(LEDhist[0])
            histdir_CT[thekey].append(CThist[0])
            histdir_M[thekey].append(Mhist[0])

    LEDfig = plt.figure()
    LEDaxes = [plt.step(xbins, np.sum(histdir_LED[thekey], axis=0), label=thekey)
                for thekey in sorted(histdir_LED.keys())]
    plt.xlabel('Time Difference (ns)', fontsize=16)
    plt.ylabel('Counts', fontsize=16)
    plt.title('LED triggers', fontsize=18)
    CTfig = plt.figure()
    CTaxes = [plt.step(xbins, np.sum(histdir_CT[thekey], axis=0), label=thekey)
                for thekey in sorted(histdir_CT.keys())]
    plt.xlabel('Time Difference (ns)', fontsize=16)
    plt.ylabel('Counts', fontsize=16)
    plt.title('Cosmic triggers', fontsize=18)
    Mfig = plt.figure()
    Maxes = [plt.step(xbins, np.sum(histdir_M[thekey], axis=0), label=thekey)
                for thekey in sorted(histdir_M.keys())]
    plt.xlabel('Time Difference (ns)', fontsize=16)
    plt.ylabel('Counts', fontsize=16)
    plt.title('Multiplicity triggers', fontsize=18)
    return histdir_LED, histdir_CT, histdir_M
    
def plotreltimeWFD(DiffSeries, evtcut, label, paramcut='time',
                   nbins=200, binrange=(-50,50), newfig=False):
    tablelist = DiffSeries.loc[pandas.IndexSlice[:, evtcut, paramcut]].dropna().values
    vallist = [tablelist[i].entries.values for i in range(len(tablelist))]
    if newfig:
        plt.figure()
    thehist = plt.hist(np.hstack(vallist), bins=nbins, range=binrange,
                       histtype='step', label=label)
    return thehist

def getmeanstd_reltimeWFD(DiffSeries, evtcut, paramcut='time', lowcut=-100, highcut=100):
    tablelist = DiffSeries.loc[pandas.IndexSlice[:, evtcut, paramcut]].dropna().values
    vallist = [item.entries.values for item in tablelist]
    flatvals = np.hstack(vallist)
    cutvals = [val for val in flatvals if val<lowcut and val>highcut]
    return np.mean(cutvals), np.std(cutvals)
    
def plot_TorQ(TDCdata, nbins=100, binrange=None):
    axes = TDCdata.hist(bins=nbins, range=binrange, histtype='step')
    for ax in plt.gcf().get_axes():#np.nditer(axes):
        ax.set_yscale('log')
    return axes

def plot_TorQ_wSlice(thedata, chanslice=slice(None), evtslice=slice(None),
                     runslice=slice(None), nbins=100, binrange=None,
                     setlog=False, plottitle='', subplots=True, xlabel='',
                     ylim=None, xlim=None, **kwargs):
    '''
    Histogram all TDC or QDC data for a given slice.
    
    theslice should be an instance of pandas.IndexSlice
    '''
    if subplots:
        axes = thedata.loc[pandas.IndexSlice[runslice,
                                             evtslice],
                                             chanslice].hist(bins=nbins,
                                                             range=binrange,
                                                             histtype='step',
                                                             **kwargs)
    else:
        axes = thedata.loc[pandas.IndexSlice[runslice,
                                             evtslice],
                                             chanslice].plot(kind='hist',
                                                             bins=nbins,
                                                             range=binrange,
                                                             histtype='step',
                                                             **kwargs)        
    for ax in plt.gcf().get_axes():#np.nditer(axes):
        if setlog:
            try:
                ax.set_yscale('log')
            except ValueError:
                pass
        ax.set_ylabel('counts', fontsize=16)
        ax.set_xlabel(xlabel, fontsize=16)
        ax.set_ylim(ylim)
        ax.set_xlim(xlim)
    if plottitle: plt.gcf().suptitle(plottitle, fontsize=18)
    return axes

def plot_Series(thedata, evtslice=slice(None),
                     runslice=slice(None), nbins=100, binrange=None,
                     setlog=False, plottitle='', subplots=True, xlabel='',
                     ylim=None, xlim=None, **kwargs):
    '''
    Histogram all TDC or QDC data for a given slice.
    
    theslice should be an instance of pandas.IndexSlice
    '''
    if subplots:
        axes = thedata.loc[pandas.IndexSlice[runslice,
                                             evtslice]].hist(bins=nbins,
                                                             range=binrange,
                                                             histtype='step',
                                                             **kwargs)
    else:
        axes = thedata.loc[pandas.IndexSlice[runslice,
                                             evtslice]].plot(kind='hist',
                                                             bins=nbins,
                                                             range=binrange,
                                                             histtype='step',
                                                             **kwargs)        
    for ax in plt.gcf().get_axes():#np.nditer(axes):
        if setlog:
            try:
                ax.set_yscale('log')
            except ValueError:
                pass
        ax.set_ylabel('counts')
        ax.set_xlabel(xlabel)
        ax.set_ylim(ylim)
        ax.set_xlim(xlim)
    if plottitle: plt.gcf().suptitle(plottitle)
    return axes

def plot_Series_WFD(thedata, evtslice=slice(None),
                     runslice=slice(None), paramslice='time', nbins=100, binrange=None,
                     setlog=False, plottitle='', subplots=True, xlabel='',
                     ylim=None, xlim=None, **kwargs):
    '''
    Histogram all TDC or QDC data for a given slice.
    
    theslice should be an instance of pandas.IndexSlice
    '''
    if subplots:
        axes = thedata.loc[pandas.IndexSlice[runslice,
                                             evtslice, paramslice]].hist(bins=nbins,
                                                             range=binrange,
                                                             histtype='step',
                                                             **kwargs)
    else:
        axes = thedata.loc[pandas.IndexSlice[runslice,
                                             evtslice, paramslice]].plot(kind='hist',
                                                             bins=nbins,
                                                             range=binrange,
                                                             histtype='step',
                                                             **kwargs)        
    for ax in plt.gcf().get_axes():#np.nditer(axes):
        if setlog:
            try:
                ax.set_yscale('log')
            except ValueError:
                pass
        ax.set_ylabel('counts')
        ax.set_xlabel(xlabel)
        ax.set_ylim(ylim)
        ax.set_xlim(xlim)
    if plottitle: plt.gcf().suptitle(plottitle)
    return axes

def plotWFD(WFDdata, ch, param, nbins=100, binrange=None):
    '''
    Simply histogram all of the WFD data of a given channel, parameter
    '''
    vals = WFDdata.loc[pandas.IndexSlice[:,:,param], ch].dropna().values
    fig = plt.figure()
    if binrange is not None:
        plt.hist(np.hstack(vals), bins=nbins, range=binrange, histtype='step')
    else:
        plt.hist(np.hstack(vals), bins=nbins)
    fig.gca().set_yscale('log')
    return

def getLEDevts(TDCdata):
    ''''Return a boolean list of events with LED triggers'''
    return np.all([pandas.isnull(TDCdata['CT']), pandas.isnull(TDCdata['M'])], axis=0)

def getCTevts(TDCdata):
    '''Return a boolean list of events with cosmic triggers'''
    return pandas.notnull(TDCdata['CT'])

def getMevts(TDCdata):
    '''Return a boolean list of events with cosmic triggers'''
    return pandas.notnull(TDCdata['M'])

def getCTandMevts(TDCdata):
    '''Return a boolean list of events with cosmic triggers'''
    return np.all([pandas.notnull(TDCdata['CT']), pandas.notnull(TDCdata['M'])], axis=0)

def WFDLEDevts(LEDevts):
    '''Formats the output of getLEDevents for use by the WFD dataframe'''
    return [val for val in LEDevts for i in range(6)]

def plotWFDwSlice(WFDdata, ch, paramslice, evtslice=slice(None),
                  runslice=slice(None), nbins=100, binrange=None,
                  figref=None, plottitle='', pltlabel=None, legend=False):
    '''
    Slice on the WFD channels, parameters, events, and runs.
    '''
    tablelist = WFDdata.loc[pandas.IndexSlice[runslice,evtslice,paramslice], ch].dropna().values
    vals = [tablelist[i].entries.values for i in range(len(tablelist))]
    #if figref is None:
    #    fig = plt.figure()
    #    fig.add_axes()
    #    ax = fig.get_axes()
    #else:
    #    ax = figref.get_axes()
    plt.hist(np.hstack(vals), bins=nbins, range=binrange,
             histtype='step', label=pltlabel)
    ax = plt.gca()
    ax.set_yscale('log')
    if plottitle: plt.gcf().suptitle(plottitle)
    if legend: plt.legend()
    return vals

def plotAllWFDwSlice(WFDdata, paramslice, chanlist=False, **kwargs):
    if not chanlist:
        chanlist = WFDdata.keys()
    fig = plt.figure()
    for ch in [key for key in WFDdata.keys()
               if 'NoPulses' not in key and key in chanlist]:
        plotWFDwSlice(WFDdata, ch, paramslice=paramslice,
                      figref=fig, pltlabel=ch, **kwargs)
    return