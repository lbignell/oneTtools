# -*- coding: utf-8 -*-
"""
Created on Mon Aug 15 10:33:08 2016

@author: lbignell
"""
import oneton.getPositions
import numpy as np
from scipy.spatial.distance import euclidean

def LED_TOF(LEDloc, rindex=1.3390):
    '''
    This function will return a dictionary containing the TOF in ns from the
    LED each PMT as a dictionary, for a given LED location (lengths are cm).
    Dispersion is ignored, so the refractive index is a constant.
    (The default is water at 400 nm).
    
    Note LEDloc should just be the channel name as a string.
    '''
    pos = oneton.getPositions.getPosition()
    c = 29.9792458 #cm/ns
    chans = ['S{0}'.format(i) for i in range(6)]
    LEDcoords = pos.getData('LEDat{0}'.format(LEDloc))
    #print("LED coords = {0}".format(LEDcoords))
    PMTpos = {ch:pos.getData(ch) for ch in chans}
    #print("PMTpos = {0}".format(PMTpos))
    return {ch:rindex*euclidean(LEDcoords, PMTpos[ch][0:3])/c for ch in chans}
