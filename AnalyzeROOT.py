# -*- coding: utf-8 -*-
"""
Created on Tue Feb  23 15:28:05 2016

@author: lbignell
"""
import numpy as np
import scipy as sp
import ROOT
import matplotlib.pyplot as plt

def GetHists(rootfile, prefix, chlist, colorlist=None):
    '''
    Get a list of hists with name format prefix+chlist
    '''
    hists = []
    if colorlist is None:
        colorlist = [i+1 for i in range(len(chlist))]
    for i,ch in enumerate(chlist):
        thehist = getattr(rootfile, prefix+ch)
        thehist.SetLineColor(colorlist[i])
        hists.append(thehist)
    return hists

def SaveCanvas(canvas, filepath, basename, formats=['.svg', '.png', '.C']):
    '''
    Save a ROOT canvas in filepath with name basename.
    '''
    for ext in formats:
        canvas.SaveAs(filepath+basename+ext)
    return

def plotHists(histlist, xlabel="", ylabel="", newcanvas=True,
              legtitle="", leglabels=None):
    if newcanvas:
        thecanvas = ROOT.TCanvas()
    else:
        thecanvas = ROOT.gPad
    first = True
    for hist in histlist:
        if first:
            hist.Draw("")
            hist.SetStats(0)
            hist.GetXaxis().SetTitle(xlabel)
            hist.GetXaxis().SetTitleSize(0.05)
            hist.GetYaxis().SetTitle(ylabel)
            hist.GetYaxis().SetTitleSize(0.05)
            firsthist = hist
            first = False
        else:
            hist.Draw("SAME")
    leg = ROOT.TLegend(0.1, 0.7, 0.3, 0.9)
    if legtitle:
        leg.SetHeader(legtitle)
    if len(histlist) == len(leglabels):
        [leg.AddEntry(histlist[i], leglabels[i]) for i in range(len(histlist))]
        leg.Draw()
    thecanvas.Update()
    thetitle = thecanvas.FindObject('title')
    thetitle.Clear()
    thecanvas.Update()
    return histlist, firsthist, leg
    
def plot2DbyChannel(thehist, normed=True, chlist=[], newfig=True, xlim=False):
    numchans = thehist.GetNbinsY()
    numxbins = thehist.GetNbinsX()
    chanlabels = [thehist.GetYaxis().GetBinLabel(i+1) for i in range(numchans)]
    valuesdict = {label: [thehist.GetCellContent(xval, yval+1)
                          for xval in range(numxbins+2)]
                          for yval, label in enumerate(chanlabels)}
    chkeys = sorted(valuesdict.keys())
    if chlist:
        chkeys = sorted(set(chkeys).intersection(chlist))
    xvals = [thehist.GetXaxis().GetBinCenter(i) for i in range(numxbins+2)]
    if newfig:
        fig = plt.figure()
    else:
        fig = plt.gcf()
    if normed:
        normfactors = {chkey : sum(valuesdict[chkey]) for chkey in chkeys}
    else:
        normfactors = [1 for chkey in chkeys]
    axes = [plt.step(xvals,np.divide(valuesdict[chkey], normfactors[chkey]), label=chkey)
            for chkey in chkeys]
    plt.legend(loc='upper left', fontsize=16)
    if xlim:
        plt.xlim(xlim)
    plt.xlabel('Charge (AU)', fontsize=16)
    plt.ylabel('Counts', fontsize=16)
    return fig, axes