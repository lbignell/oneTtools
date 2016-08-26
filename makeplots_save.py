# -*- coding: utf-8 -*-
"""
Created on Tue May 31 09:55:54 2016

@author: lbignell
"""
import shelve
import pandas
import oneTtools.PandasAnalysis
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from scipy.optimize import curve_fit
from scipy.optimize import leastsq
import itertools
import re
from argparse import ArgumentParser
import os
import json
import sys
from scipy import signal

def makeplots_LED_save(path, runnums, nevts, relch, savepath='', 
                       LEDloc='unspecified', cutdict={}, formats=['.svg'],
                        excludeMultiplicity=False, **kwargs):
    '''
    make plots of TDC and WFD timing relative to relch, then save data using shelve.

    Args:
    - Path to data folder (string).
    - Run numbers of data to plot (string list).
    - Number of events in the analyis.
    - String name of channel to be made relative to.
    - String title plot.
    - Path to save results.
    '''
    TDClist = []
    WFDlist = []
    runkeys = []
    for run in runnums:
        fname = path + 'runs{0}--{0}_{1}evts.shelf'.format(run, nevts)
        print('file: {0}'.format(fname))
        theshelf = shelve.open(fname, flag='r')
        TDClist += theshelf['TDCdata']
        WFDlist += theshelf['WFDdata']
        runkeys += theshelf['runkeys']
        theshelf.close()
    
    TDCdata = pandas.concat(TDClist, keys=runkeys)
    WFDdata = pandas.concat(WFDlist, keys=runkeys)
    LEDevts = oneTtools.PandasAnalysis.getLEDevts(TDCdata)
    CTevts = oneTtools.PandasAnalysis.getCTevts(TDCdata)
    Mevts = oneTtools.PandasAnalysis.getMevts(TDCdata)
    LEDevts_WFD = oneTtools.PandasAnalysis.WFDLEDevts(LEDevts)
    CTevts_WFD = oneTtools.PandasAnalysis.WFDLEDevts(CTevts)
    notCTevts_WFD = [not(evt) for evt in CTevts_WFD]
    Mevts_WFD = oneTtools.PandasAnalysis.WFDLEDevts(Mevts)
    #rel_TDC = {'S{0}'.format(i): getattr(TDCdata, 'S{0}'.format(i)) -
    #            getattr(TDCdata, relch) for i in range(6)}
    #rel_WFD = {'S{0}'.format(i): getattr(WFDdata, 'S{0}'.format(i)) -
    #            getattr(WFDdata, relch) for i in range(6)}
    #fig = plt.figure()
    #axes_LED = [oneTtools.PandasAnalysis.plot_Series(rel_TDC[key], 
    #                                                 evtslice=LEDevts,
    #                                                 nbins=400,
    #                                                 binrange=(-100,100),
    #                                                 subplots=False,
    #                                                 xlabel='TDC time (ns)',
    #                                                 label=key) for key in
    #                                                 sorted(rel_TDC.keys())]
    #fig = plt.figure()
    #axes_WFD = [oneTtools.PandasAnalysis.plotreltimeWFD(rel_WFD[key],
    #                                                    evtcut=LEDevts_WFD,
    #                                                    label=key,
    #                                                    nbins=400,
    #                                                    binrange=(-100,100),
    #                                                    newfig=False)
    #                                                    for key in
    #                                                    sorted(rel_WFD.keys())]
    
    chans = [key for key in TDCdata.keys() if re.match('S[0-5]', key)]
    fitresults_TDC = diffplots_TDC(TDCdata, chans, relch, LEDloc=LEDloc, 
                                   nbins=1000, binrange=(-100,100), xlim=(-40,40),
                                   savepath=savepath, formats=formats, 
                                   extraname='_TDC')

    #Since the timing calibration run uses a very large pulse, it isn't wise
    #to exclude multiplicity triggers.
    if excludeMultiplicity:
        evtcuts = LEDevts_WFD
    else:
        evtcuts = notCTevts_WFD
        
    fitresults_WFD = mkplts_WFD(WFDdata, evtcuts, relch, LEDloc=LEDloc,
                                cutdict=cutdict, savepath=savepath, 
                                formats=formats, nbins=2000)

    return fitresults_TDC, fitresults_WFD

def mkplts_WFD(WFDdata, LEDevts_WFD, refch, LEDloc='', cutdict={}, savepath='',
               formats=['.svg'], extraname='', cutsintitle=False, nbins=2000):
    '''
    This function will make some time vs charge plots for each channel in the
    given WFDdata dataset, then ask the user to supply some cuts based on the 
    plots (if no cuts were given).
    Then it will plot the charge, time, and time difference distributions for
    the cut data, relative to refch.
    '''
    fitresults = {}
    chans = [key for key in WFDdata.keys() if re.match('S[0-5]', key)]

    for ch in chans:
        timevschargeWFD(WFDdata, ch, LEDloc=LEDloc, evtslice=LEDevts_WFD, 
                        xlim=(0,500), savepath=savepath, formats=formats,
                        extraname=extraname)

    if not cutdict:
        print('No cuts have been selected. You can enter them now. If you \
                enter nothing, no cuts will be applied.\n')
        cutdict['timecutlow'] = float(input('\nearly absolute time cut (ns): ')
                                      or -np.inf)
        cutdict['timecuthigh'] = float(input('\nlate absolute time cut (ns): ')
                                       or np.inf)
        cutdict['areacutlow'] = float(input('\nlow charge cut (same arb units as plot): ')
                                      or -np.inf)
        cutdict['areacuthigh'] = float(input('\nhigh charge cut (same arb units as plot): ')
                                       or np.inf)
        cutdict['reltimecutlow'] = float(input('\nearly time cut relative to reference channel (ns): ')
                                         or -np.inf)
        cutdict['reltimecuthigh'] = float(input('\nlate time cut relative to reference channel (ns): ')
                                          or np.inf)
    
    for ch in chans:
        ch_time_cut, ch_area_cut, ref_time_cut, ref_area_cut, timediff_cut,\
        cutdict_applied = ApplyCutsWFD(WFDdata, ch, refch, LEDloc='', 
                                       evtslice=LEDevts_WFD, 
                                       timecutlow = cutdict['timecutlow'],
                                       timecuthigh = cutdict['timecuthigh'],
                                       areacutlow = cutdict['areacutlow'],
                                       areacuthigh = cutdict['areacuthigh'],
                                       reltimecutlow = cutdict['reltimecutlow'],
                                       reltimecuthigh = cutdict['reltimecuthigh'])
        #plot the cut time distribution.
        #plot_hist(ch_time_cut, ch, refch, LEDloc=LEDloc, cutdict=cutdict, 
        #          nbins=2000, binrange=(0,500), xlim=False, 
        #          savepath=savepath, module='WFD', plttype = 'time',
        #          units='ns', formats=formats, extraname=extraname)
        #plot_hist(ch_area_cut, ch, refch, LEDloc=LEDloc, cutdict=cutdict, 
        #          nbins=2000, binrange=None, xlim=False, 
        #          savepath=savepath, module='WFD', plttype = 'area',
        #          units='AU', formats=formats, extraname=extraname)
        try:
            thedict = diffplot_WFD(timediff_cut, ch, refch, LEDloc=LEDloc,
                                   xlim=(cutdict['reltimecutlow'], cutdict['reltimecuthigh']),
                                   cutdict={}, nbins=nbins, binrange=(-100,100), 
                                   savepath=savepath, formats=formats, 
                                   extraname=extraname)
            fitresults[ch] = thedict
        except ValueError:
            print("mkplotsWFD: Couldn't plot channel {0} WFD data".format(ch))
            print("ErrMsg:", sys.exc_info()[0])
            raise

    return fitresults

def CorrelationPlotsTDC(TDCdata, chlist, refch, nbins=100, absbinrange=(40,140),
                        abstitle='', xlim=False, ylim=False, savepath='',
                        formats=['.svg'], extraname=''):
    LEDevts = oneTtools.PandasAnalysis.getLEDevts(TDCdata)
    if not abstitle:
        abstitle = 'LED next to {0}'.format(refch)
    axes = oneTtools.PandasAnalysis.plot_TorQ_wSlice(TDCdata, chanslice=chlist,
                                                     evtslice=LEDevts,
                                                     nbins=nbins, binrange=absbinrange,
                                                     subplots=False, 
                                                     xlabel='TDC time (ns)',
                                                     plottitle=abstitle)
    if savepath:
        fname = savepath + 'AbsTDC_LEDNextTo{0}'.format(refch)
        if extraname:
            fname += extraname
        fig = plt.gcf()
        for fmt in formats:
            fig.savefig(fname+fmt)
    for ch in chlist:
        fig = plt.figure()
        axes = plt.hist2d(TDCdata.loc[pandas.IndexSlice[slice(None), LEDevts], refch],
                          TDCdata.loc[pandas.IndexSlice[slice(None), LEDevts], ch],
                          bins=[200,200], range=((40,140), (40,140)),
                          norm=mpl.colors.LogNorm())
        plt.xlabel('{0} TDC Time (ns)'.format(refch), fontsize=16)
        plt.ylabel('{0} TDC Time (ns)'.format(ch), fontsize=16)
        ax = plt.gca()        
        ax.xaxis.set_tick_params(labelsize=16)
        ax.yaxis.set_tick_params(labelsize=16)
        plt.title('{0} - {1} correlation, TDC'.format(refch, ch), fontsize=18)
        cbar = plt.colorbar()
        cbar.set_label('Counts', fontsize=16)
        plt.tight_layout()
        if xlim: plt.xlim(xlim)
        if ylim: plt.ylim(ylim)
        if savepath:
            fname = savepath + 'Corr{0}vs{1}'.format(refch, ch)
            if extraname:
                fname += extraname
            for fmt in formats:
                fig.savefig(fname+fmt)
    
    return

def diffplots_TDC(TDCdata, chlist, refch, LEDloc='unspecified', nbins=1000, 
                  binrange=(-100,100), xlim=False, savepath='', 
                  formats=['.svg'], extraname=''):
    '''
    Make TDC plots relative to some reference channel.
    returns a dictionary containing mean and standard deviation of time differences.
    
    Arguments:
    - The TDC data, as generated by PandasAnalysis.py
    - A list of (string) channels.
    - The reference channel name (string).
    - The number of bins for the histogram.
    - The histogram bin range.
    - If xlim is specified as a 2-element tuple range, fit the histogram to a
      Gaussian within that range and plot the fitted result.
    - If a save path is given, save to this folder.
    - A list of formats to save in.
    - An extra name to be appended to the saved files.
    '''
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    LEDevts = oneTtools.PandasAnalysis.getLEDevts(TDCdata)
    fitresults = {}
    for ch in chlist:
        fig = plt.figure()
        try:
            thehist = plt.hist(TDCdata.loc[pandas.IndexSlice[slice(None), LEDevts], 
                                           ch] - 
                                           TDCdata.loc[pandas.IndexSlice[slice(None), LEDevts],
                                                       refch], bins=nbins, range=binrange, 
                                           histtype='step')
            plt.xlabel('{0} - {1} TDC Time (ns)'.format(ch, refch), fontsize=16)
            plt.ylabel('Counts', fontsize=16)
            ax = plt.gca()
            ax.xaxis.set_tick_params(labelsize=16)
            ax.yaxis.set_tick_params(labelsize=16)
            plt.xlim(binrange[0], binrange[1])
            plt.title('{0} - {1} TDC time, LED location: {2}'.format(ch, refch, LEDloc), fontsize=18)
            if savepath:
                fname = savepath + '{0}rel{1}_LED{2}'.format(ch, refch, LEDloc)
                if extraname:
                    fname += extraname
                for fmt in formats:
                    fig.savefig(fname+fmt)
            if xlim:
                minidx = list(thehist[1]).index(xlim[0])
                maxidx = list(thehist[1]).index(xlim[1])
                params = fit_function([max(thehist[0][minidx:maxidx]),
                                           np.mean(thehist[0][minidx:maxidx]),
                                           np.std(thehist[0][minidx:maxidx])],
                                      GetBinCentres(thehist)[minidx:maxidx],
                                      thehist[0][minidx:maxidx], Gauss, 
                                      sigma=np.sqrt(thehist[0][minidx:maxidx]))
                paramdict = {'mean': params[0][1], 'Umean': params[1][1],
                             'std': params[0][2], 'Ustd': params[1][2]}
                xvals=GetBinCentres(thehist)
                yvals=Gauss(xvals,params[0])
                plt.plot(xvals,yvals,'r', linewidth=2)
                ax.text(0.05, 0.95, 
                        'mean = {:.3f} +/- {:.3f}\nstd = {:.3f} +/- {:.3f}'.format(paramdict['mean'], 
                                                                                   paramdict['Umean'], 
                                                                                   paramdict['std'],
                                                                                   paramdict['Ustd']), 
                        fontsize=14,verticalalignment='top', horizontalalignment='left',
                        bbox=props, transform=ax.transAxes)
                plt.xlim(xvals[xlim[0]], xvals[xlim[1]])
                fitresults[ch] = paramdict
                if savepath:
                    fname = savepath + '{0}rel{1}_LED{2}_fitted'.format(ch, refch, LEDloc)
                    if extraname:
                        fname += extraname
                    for fmt in formats:
                        fig.savefig(fname+fmt)
        except UnboundLocalError:
            print("{0}: Couldn't analyse channel {1} TDC data".format(__file__, ch))
            print("ErrMsg:", sys.exc_info()[0])

    return fitresults

#def CorrelationPlotWFD(WFDdata, ch, refch, paramslice='time', evtslice=slice(None), 
#                       runslice=slice(None), nbins= [2560,2560], pltlabel='',
#                       hrange=None, xlim=False, ylim=False):
#    tablelist = WFDdata.loc[pandas.IndexSlice[runslice,evtslice,paramslice], ch].dropna().values
#    chvals = [tablelist[i].entries.values for i in range(len(tablelist))]
#    tablelist = WFDdata.loc[pandas.IndexSlice[runslice,evtslice,paramslice], refch].dropna().values
#    refchvals = [tablelist[i].entries.values for i in range(len(tablelist))]
#    fig = plt.figure()
#    ret = plt.hist2d(np.hstack(refchvals), np.hstack(chvals), bins=nbins, range=hrange,
#                     norm=mpl.colors.LogNorm())
#    plt.xlabel('{0} WFD time (ns)'.format(refch), fontsize=16)
#    plt.ylabel('{0} WFD time (ns)'.format(ch), fontsize=16)
#    plt.title('{0} - {1} correlation, WFD'.format(refch, ch), fontsize=18)
#    cbar = plt.colorbar()
#    cbar.set_label('Counts', fontsize=16)
#    plt.tight_layout()
#    if xlim: plt.xlim(xlim)
#    if ylim: plt.ylim(ylim)
#    return

def plot_hist(values, ch, refch, LEDloc='', cutdict={}, 
              nbins=2000, binrange=(-100,100), xlim=False, 
              savepath='', module='module?', plttype = 'unspecified',
              units='units?', formats=['.svg'], extraname=''):
    fig = plt.figure()
    ret = plt.hist(values, bins=nbins, range=binrange, 
                   histtype='step')
    plt.xlabel('{0} - {1} {2} {3} ({4})'.format(ch, refch, module, plttype, units), 
               fontsize=16)
    plt.ylabel('Counts', fontsize=16)
    ax = plt.gca()        
    ax.xaxis.set_tick_params(labelsize=16)
    ax.yaxis.set_tick_params(labelsize=16)
    plttit = 'ch {0}, {1} {2}, LED location: {3}'.format(ch,module,plttype,LEDloc)
    if cutdict:
        cutstring = str(cutdict)[1:-1]
        plttit += '\n'+cutstring[1:40] + '\n' + cutstring[41:80] + '\n' + cutstring[81:-1]
    plt.title(plttit, fontsize=18)
    if xlim: plt.xlim(xlim)
    if savepath:
        fname = savepath + '{0}_{1}{2}_LED{3}'.format(ch, module, plttype, LEDloc)
        if extraname:
            fname += extraname
        for fmt in formats:
            fig.savefig(fname+fmt)

    return

def diffplot_WFD(timediff_cut, ch, refch, LEDloc='', cutdict={}, 
                 nbins=2000, binrange=(-100,100), xlim=False, 
                 savepath='', formats=['.svg'], extraname=''):
    fig = plt.figure()
    try:
        thehist = plt.hist(timediff_cut, bins=nbins, range=binrange, 
                           histtype='step')
    except ValueError:
        print("{0} diffplot_WFD: ValueError on hist command. sum(timediff_cut) = {1}"
            .format(__file__, sum(timediff_cut)))
        return
    plt.xlabel('{0} - {1} WFD Time (ns)'.format(ch, refch), fontsize=16)
    plt.ylabel('Counts', fontsize=16)
    ax = plt.gca()        
    ax.xaxis.set_tick_params(labelsize=16)
    ax.yaxis.set_tick_params(labelsize=16)
    plttit = '{0} - {1} WFD time, LED location: {2}'.format(ch,refch,LEDloc)
    paramdict = {}
    if cutdict:
        cutstring = str(cutdict)[1:-1]
        plttit += '\n'+cutstring
    plt.title(plttit, fontsize=18)
    if savepath:
        fname = savepath + '{0}rel{1}_LED{2}_WFD'.format(ch, refch, LEDloc)
        if extraname:
            fname += extraname
        for fmt in formats:
            fig.savefig(fname+fmt)
    if xlim:
        histy = np.array(thehist[0])
        histx = np.array(thehist[1])
        fitrange = (histx>xlim[0]) & (histx<xlim[1])
        fityvals = histy[fitrange]
        xvals=np.array(GetBinCentres(thehist))
        xvals = xvals[fitrange]
        rangemean = sum(xvals*fityvals/sum(fityvals))
        print('estimated mean = {0}'.format(rangemean))
        rangestd = (sum(((xvals**2)*fityvals)/sum(fityvals)) - rangemean**2)**0.5
        print('estimated std = {0}'.format(rangestd))
        paramestimate = [max(fityvals), rangemean, rangestd]
        params = fit_function(paramestimate, xvals, fityvals, Gauss, 
                              sigma=np.sqrt(fityvals))
        paramdict = {'mean': params[0][1], 'Umean': params[1][1],
                     'std': params[0][2], 'Ustd': params[1][2]}
        means, ampls = findpeaks(fityvals, xvals)
        yvals=Gauss(xvals,params[0])
        plt.plot(xvals,yvals,'r', linewidth=2)
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax.text(0.05, 0.95, 
                'mean = {:.3f} +/- {:.3f}\nstd = {:.3f} +/- {:.3f}'.format(paramdict['mean'], 
                                                                           paramdict['Umean'], 
                                                                           paramdict['std'],
                                                                           paramdict['Ustd']), 
                fontsize=14,verticalalignment='top', horizontalalignment='left',
                bbox=props, transform=ax.transAxes)
        pkfindfit=np.zeros(len(xvals))
        for i,mean in enumerate(means):
            if ampls[i]/max(ampls) > 0.05 and ampls[i]>10:
                thepeak = Gauss(xvals,[ampls[i], mean, 0.75])
                pkfindfit += np.array(thepeak)
        if sum(pkfindfit)>0:
            plt.plot(xvals, pkfindfit, 'g', linewidth=2)
        plt.xlim(xlim)
        print('{0} relative to {1}, peakfind results:\nmeans = {2}\namplitudes = {3}'\
            .format(ch, refch, means, ampls))
        if savepath:
            fname = savepath + '{0}rel{1}_LED{2}_fitted_WFD'.format(ch, refch, LEDloc)
            if extraname:
                fname += extraname
            for fmt in formats:
                fig.savefig(fname+fmt)

    return paramdict

def findpeaks(yvals, bincentres, minstd=0.5, maxstd=2.5):
    binwidth = bincentres[1] - bincentres[0]
    indicies = signal.find_peaks_cwt(yvals, np.linspace(0.5/binwidth, 2.5/binwidth, 10000))
    means = [bincentres[idx] for idx in indicies]
    ampls = [yvals[idx] for idx in indicies]
    return means, ampls

def timevschargeWFD(WFDdata, ch, LEDloc='', evtslice=slice(None), 
                    runslice=slice(None), nbins=[2560,2000], hrange=None, 
                    xlim=False, ylim=False,
                    savepath='', formats=['.svg'], extraname=''):
    '''
    Make a plot of time vs charge for a given channel.
    '''
    times = getvals(WFDdata, 'time', runslice, evtslice, ch)
    areas = getvals(WFDdata, 'area', runslice, evtslice, ch)
    fig = plt.figure()
    ret = plt.hist2d(np.hstack(times), np.hstack(areas), bins=nbins, range=hrange,
                     norm=mpl.colors.LogNorm())
    plt.xlabel('{0} WFD time (ns)'.format(ch), fontsize=16)
    plt.ylabel('{0} WFD charge (AU)'.format(ch), fontsize=16)
    ax = plt.gca()
    ax.xaxis.set_tick_params(labelsize=16)
    ax.yaxis.set_tick_params(labelsize=16)
    plt.title('LED location: {0}'.format(LEDloc), fontsize=18)
    cbar = plt.colorbar()
    cbar.set_label('Counts', fontsize=16)
    plt.tight_layout()
    if xlim: plt.xlim(xlim)
    if ylim: plt.ylim(ylim)
    if savepath:
        fname = savepath + '{0}_WFDChargeVsTime_LED{1}'.format(ch, LEDloc)
        if extraname:
            fname += extraname
        for fmt in formats:
            fig.savefig(fname+fmt)
    return

def ApplyCutsWFD(WFDdata, ch, refch, LEDloc='', evtslice=slice(None), 
                 timecutlow = -np.inf, timecuthigh = np.inf,
                 areacutlow = -np.inf, areacuthigh = np.inf,
                 reltimecutlow = -np.inf, reltimecuthigh = np.inf,                    
                 runslice=slice(None)):
    '''
    Apply cuts to the WFD data for ch, relative to ref ch.
    returns: ch time with cuts, ch area with cuts,
             ref time with cuts, ref area with cuts,
             ch time - ref time with cuts.
    '''
    chvals_time = getvals(WFDdata, 'time', runslice, evtslice, ch)
    chvals_area = getvals(WFDdata, 'area', runslice, evtslice, ch)
    refchvals_time = getvals(WFDdata, 'time', runslice, evtslice, refch)
    refchvals_area = getvals(WFDdata, 'area', runslice, evtslice, refch)
    #I need to cut the values, keeping the event info (that is encoded by the
    #2D nature of the chvals arrays). I then need to subtract every combination
    #on an event-by-event basis and throw away the data where the subtracted value
    #is unphysically big (say, > 20 ns).
    ch_time_cut = applycuts(chvals_time, chvals_area, 
                                timecutlow, timecuthigh,
                                areacutlow, areacuthigh)
    ch_area_cut = applycuts(chvals_area, chvals_time,
                                areacutlow, areacuthigh,
                                timecutlow, timecuthigh)
    ref_time_cut = applycuts(refchvals_time, refchvals_area,
                                   timecutlow, timecuthigh,
                                   areacutlow, areacuthigh)
    ref_area_cut = applycuts(refchvals_area, refchvals_time,
                                   areacutlow, areacuthigh,
                                   timecutlow, timecuthigh)

    timediff_cut = applydiffcut(ch_time_cut, ref_time_cut,
                                reltimecutlow, reltimecuthigh)
    return ch_time_cut, ch_area_cut, ref_time_cut, ref_area_cut, timediff_cut,\
        {'timecutlow': timecutlow, 'timecuthigh': timecuthigh,
         'areacutlow': areacutlow, 'areacuthigh': areacuthigh,
         'reltimecutlow': reltimecutlow, 'reltimecuthigh': reltimecuthigh}

def applydiffcut(vals1, vals2, lowdiffcut, highdiffcut):
    '''
    Return an array of the difference between the first and second argument,
    where the difference is between the cut levels.
    '''
    return [pair[0] - pair[1] for evt in range(len(vals1)) 
            for pair in itertools.product(vals1[evt],vals2[evt]) 
            if ((pair[0]-pair[1])>lowdiffcut) and
               ((pair[0]-pair[1])<highdiffcut)]    

def applycuts(vals1, vals2, lowcut1, highcut1, lowcut2, highcut2):
    '''
    This will cut two parameters that have one-to-one correspondence (for
    instance, the pulse area and time for a single digitizer channel).
    vals1 will be returned, being cut on low and high values of itself and
    vals2.
    '''
    return [vals1[evt][(vals1[evt] > lowcut1) & (vals1[evt] < highcut1) &
                       (vals2[evt] > lowcut2) & (vals2[evt] < highcut2)] 
            for evt in range(len(vals1))]

def getvals(WFDdata, param, runslice, evtslice, ch):
    chlist = WFDdata.loc[pandas.IndexSlice[runslice, evtslice, param], ch].values
    return [getentries_safe(frame) for frame in chlist]

def getentries_safe(frame):
    try:
        vals = frame.entries.values
    except AttributeError:
        if not np.isnan(frame):
            print('WARNING: unhandled frame. Frame: {0}'.frame)
        vals = np.array([], dtype=np.float64)
    return vals

def GetBinCentres(hist1D):
    return (hist1D[1][:-1] + hist1D[1][1:])/2

def Gauss(x, p):
        #A1, mu1, sigma1, A2, mu2, sigma2 = p
        return p[0]*np.exp(-(x-p[1])**2/(2.*p[2]**2))

def Gauss_cf(x, p0, p1, p2):
        #A1, mu1, sigma1, A2, mu2, sigma2 = p
        return p0*np.exp(-(x-p1)**2/(2.*p2**2))
    
def fit_function(p0, datax, datay, function, sigma=None,
                 cf_func=Gauss_cf, **kwargs):
        
        errfunc = lambda p, x, y,: function(x,p) - y
        #print "p0 = ", p0

        ##################################################
        ## 1. COMPUTE THE FIT AND FIT ERRORS USING leastsq
        ##################################################

        # If using optimize.leastsq, the covariance returned is the 
        # reduced covariance or fractional covariance, as explained
        # here :
        # http://stackoverflow.com/questions/14854339/in-scipy-how-and-why-does-curve-fit-calculate-the-covariance-of-the-parameter-es
        # One can multiply it by the reduced chi squared, s_sq, as 
        # it is done in the more recenly implemented scipy.curve_fit
        # The errors in the parameters are then the square root of the 
        # diagonal elements.   
        pfit, pcov, infodict, errmsg, success = \
            leastsq( errfunc, p0, args=(datax, datay), \
                              full_output=1)

        if (len(datay) > len(p0)) and pcov is not None:
            s_sq = (errfunc(pfit, datax, datay)**2).sum()/(len(datay)-len(p0))
            pcov = pcov * s_sq
        else:
            pcov = np.inf

        error = [] 
        for i in range(len(pfit)):
            try:
                error.append( np.absolute(pcov[i][i])**0.5)
            except:
                error.append( 0.00 )
        pfit_leastsq = pfit
        perr_leastsq = np.array(error) 
        #print "cov matrix = ", pcov


        ###################################################
        ## 2. COMPUTE THE FIT AND FIT ERRORS USING curvefit
        ###################################################
    
        # When you have an error associated with each dataY point you can use 
        # scipy.curve_fit to give relative weights in the least-squares problem. 
        datayerrors = sigma#kwargs.get('datayerrors', None)
        curve_fit_function = cf_func#kwargs.get('curve_fit_function', function)
        if datayerrors is None:
            try:
                pfit, pcov = \
                    curve_fit(curve_fit_function,datax,datay,p0=p0)
            except RuntimeError:
                pass#pfit and pcov will just be the leastsq values
        else:
            try:
                pfit, pcov = \
                     curve_fit(curve_fit_function,datax,datay,p0=p0,
                               sigma=datayerrors)
            except RuntimeError:
                pass#pfit and pcov will just be the leastsq values
        error = [] 
        for i in range(len(pfit)):
            try:
              error.append( np.absolute(pcov[i][i])**0.5)
            except:
              error.append( 0.00 )
        pfit_curvefit = pfit
        perr_curvefit = np.array(error)  


        ####################################################
        ## 3. COMPUTE THE FIT AND FIT ERRORS USING bootstrap
        ####################################################        

        # An issue arises with scipy.curve_fit when errors in the y data points
        # are given.  Only the relative errors are used as weights, so the fit
        # parameter errors, determined from the covariance do not depended on the
        # magnitude of the errors in the individual data points.  This is clearly wrong. 
        # 
        # To circumvent this problem I have implemented a simple bootstraping 
        # routine that uses some Monte-Carlo to determine the errors in the fit
        # parameters.  This routines generates random datay points starting from
        # the given datay plus a random variation. 
        #
        # The random variation is determined from average standard deviation of y
        # points in the case where no errors in the y data points are avaiable.
        #
        # If errors in the y data points are available, then the random variation 
        # in each point is determined from its given error. 
        # 
        # A large number of random data sets are produced, each one of the is fitted
        # an in the end the variance of the large number of fit results is used as 
        # the error for the fit parameters. 
    
        # Estimate the confidence interval of the fitted parameter using
        # the bootstrap Monte-Carlo method
        # http://phe.rockefeller.edu/LogletLab/whitepaper/node17.html
        #residuals = errfunc( pfit, datax, datay)
        #s_res = np.std(residuals)
        #ps = []
        # 100 random data sets are generated and fitted
        #for i in range(100):
        #    if datayerrors is None:
         #       randomDelta = np.random.normal(0., s_res, len(datay))
         #       randomdataY = datay + randomDelta
         #   else:
         #       randomDelta =  np.array( [ \
         #                        np.random.normal(0., derr + 1e-10,1)[0] \
         #                        for derr in datayerrors ] ) 
         #       randomdataY = datay + randomDelta
         #   randomfit, randomcov = \
         #       leastsq( errfunc, p0, args=(datax, randomdataY),\
         #               full_output=0)
         #   ps.append( randomfit ) 

        #ps = np.array(ps)
        #mean_pfit = np.mean(ps,0)
        #Nsigma = 1. # 1sigma gets approximately the same as methods above
                    # 1sigma corresponds to 68.3% confidence interval
                    # 2sigma corresponds to 95.44% confidence interval
        #err_pfit = Nsigma * np.std(ps,0) 

        #pfit_bootstrap = mean_pfit
        #perr_bootstrap = err_pfit


        # Print results 
        #print("\nlestsq method :")
        #print("pfit = ", pfit_leastsq)
        #print("perr = ", perr_leastsq)
        #print("\ncurvefit method :")
        #print("fit = ", pfit_curvefit)
        #print("err = ", perr_curvefit)
        #print("\nbootstrap method :")
        #print("pfit = ", pfit_bootstrap)
        #print("perr = ", perr_bootstrap)
        return pfit_curvefit, perr_curvefit#pfit_leastsq, perr_leastsq#_bootstrap

if  __name__ == '__main__' :
    parser = ArgumentParser()
    parser.add_argument("-N","--Nevents",default=1E4,type=int,
                      help="Number of events per input file to process. [default %default]")
    parser.add_argument("-r","--Runs",default=None,type=str,
                      help="String to match when selecting Runs.\nMultiple \
                      runs are comma separated, don't include spaces.")
    parser.add_argument("-c", "--refch", default="", type=str,
                      help="The reference channel (usually the LED location)")
    parser.add_argument("-l", "--LEDloc", default="unspecified", type=str,
                      help="The LED location")
    parser.add_argument('-d', '--datafolder', default=None, type=str,
                      help='path to the data folder')
    parser.add_argument('-o', '--outputfolder', default=None, type=str,
                      help='path to the output folder')
    parser.add_argument('-D', '--cutdict', default='{}', type=json.loads,
                      help='A string representation of the dictionary of cuts,\
                       to be interpreted by json.loads(cutdict).')
    #(options, args) = parser.parse_args(args=sys.argv)
    options = parser.parse_args()

    error = False

    if options.Runs is not None:
        #runs = ['{:06d}'.format(int(run)) for run in options.Runs.split(',')]
        runs = [int(run) for run in options.Runs.split(',')]
        print('Runs to be processed: {0}'.format(runs))
    else:
        print('makeplots_save ERROR: No runs selected!')
        error = True
        
    if not options.refch:
        print("makeplots_save ERROR: no reference channel specified")
        error = True

    os.makedirs(options.outputfolder, exist_ok=True)
    print('makeplots_save: options.outputfolder = {0}'.format(options.outputfolder))

    outfilename = options.outputfolder + \
        'FitResults_rel{0}_loc{1}.json'.format(options.refch, options.LEDloc)
    print('makeplots_save: outfilename = {0}'.format(outfilename))

    #cutdict = json.loads(options.cutdict)
    cutdict = options.cutdict
    print('cuts: {0}'.format(cutdict))

    fitresults_TDC = {}
    fitresults_WFD = {}

    if not error:
        if runs:
            fitresults_TDC[options.refch], fitresults_WFD[options.refch] = \
            makeplots_LED_save(options.datafolder, runs, options.Nevents, 
                               options.refch, savepath=options.outputfolder, 
                               LEDloc=options.LEDloc, cutdict=cutdict, 
                               formats=['.svg', '.png', '.pdf'])

        with open(outfilename, 'w') as fitfile:
            fitfile.write('cutdict='+json.dumps(options.cutdict)+'\n')
            fitfile.write('fitresults_TDC='+json.dumps(fitresults_TDC)+'\n')
            fitfile.write('fitresults_WFD='+json.dumps(fitresults_WFD)+'\n')
        
