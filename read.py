# -*- coding: utf-8 -*-
"""
Created on Tue Feb  23 15:28:05 2016

@author: lbignell
"""
import numpy as np
import scipy as sp
import h5py

class datafile:
    '''
    Class to unpack HDF5 data file

    Limitations (to be fixed):


    Base units: ns, GHz, fC
    '''
    def __init__(self, fname=None):
        '''
        Instanciate an object and (if path is given) unpack some basic info.
        '''
        self.fname = fname
        if self.fname is not None:
            try:
                self.h5file = h5py.File(self.fname,mode='r')
            except IOError:
                self.h5file = None
                print('ERROR: file {0} does not exist or can\'t be opened'.format(
                    self.fname))
                return
            #register the calibration data
            self.reg_caldata()
            #Get some basic info from Run_Info
            self.Run_Info = self.h5file['Run_Info']
            self.Configuration_Settings = self.Run_Info['Configuration_Settings']
            self.Nevts = self.Run_Info['Num_Evts_When_Flushed'][()]
            self.fclock = self.Configuration_Settings['Digitizer_Clock_Freq_GHz'][()]
            #Trig time is given by total acquisition time - nbcols after trig * 1/Fp
            #I'm assuming fp = 100 MHz (period = 10 ns)
            self.lenWFD = max(np.shape(self.dig1cal))
            self.TotTimeWFD = self.lenWFD/self.fclock
            self.TrigTimeWFD = self.TotTimeWFD - \
                self.Configuration_Settings['Digitizer_Num_Cols_After_Trig'][()]*10
            self.binw_QDC = (self.Configuration_Settings['QDC_High_Range_Bin_Width_fC'][()],
                             self.Configuration_Settings['QDC_Low_Range_Bin_Width_fC'][()])
            self.threshBin_QDC = (self.Configuration_Settings['QDC_High_Range_Threshold_bin'][()],
                                  self.Configuration_Settings['QDC_Low_Range_Threshold_bin'][()])
            self.binw_TDC = self.Configuration_Settings['TDC_Bin_Width_ns'][()]
            self.threshbin_TDC = self.Configuration_Settings['TDC_Threshold_bin'][()]
            self.Run_Details = self.Run_Info['Run_Details']
            self.HVon = self.Run_Details['HV_On'][()]
            self.runnum = self.Run_Details['Run_Number'][()]
            self.runtype = self.Run_Details['Run_Type'][()]
            self.Time_Info = self.Run_Info['Time_Info']
            self.tstart = self.Time_Info['Start_Time_UNIX'][()]
        return

    def set_fname(self, fname):
        self.fname = fname

    def reg_caldata(self):        
        '''
        This method gets the digitizer calibration data and stores it in the
        digNcal arrays.
        '''
        self.dig1cal = self.h5file['Calibration']['Digitizer_1']['Pedestal'][:]
        self.dig2cal = self.h5file['Calibration']['Digitizer_2']['Pedestal'][:]
        return

    def setbranch(self, branch):
        self.branch = branch

    def _effmodel(self, KE, Mparticle, Zparticle, factor,verbose=0):
        return (1 - np.e**(-1*self.CollEff*self.QE*self.scint.LY*
            self.scint.quenched_en(KE, Mparticle, Zparticle,verbose-1)[0]*factor))

    def eff_nPMT_beta(self, n, factor=1, verbose=0):
        '''
        Calculate the efficiency of 1 PMT to the registered scintillator and
        beta decay branch.
        
        Arguments:

        - number of PMTs

        - Optional factor to model efficiency extrapolation.
        '''
        if self.scint is None:
            print('ERROR!! No scintillator has been registered.')
            return None
        elif self.branch is None:
            print('ERROR!! No branch has been registered.')
            return None
        eff = lambda KE, Mparticle, Zparticle, n:\
            self.branch.beta_spec(KE)*self._effmodel(KE,Mparticle,Zparticle,factor,verbose-1)**n
        #I should add a verbosity option later to plot the efficiency spectrum.
        unnorm = sp.integrate.quad(eff, 0, self.branch.Q, 
                                   args=(self.branch.Melec,self.branch.betacharge,n))
        if verbose>0:
            print('TDCR.eff_nPMT_beta: n = {0}, factor= {1}, integral = {2} +/- {3}'.format(
                n, factor, unnorm[0], unnorm[1]))
        return unnorm[0]/sp.integrate.quad(self.branch.beta_spec, 0, self.branch.Q)[0]
        
    def eff_nPMT_monoenergetic(self, KE, Mparticle, Zparticle, n, factor=1):
        '''
        Calculate the efficiency of 1 PMT to the registered scintillator and
        a monoenergetic emission.
        
        Arguments:

        - Particle energy (MeV)

        - Particle rest mass (MeV/c^2)

        - Particle charge (e)

        - number of PMTs

        - Optional factor to model efficiency extrapolation.
        '''
        if self.scint is None:
            print('ERROR!! No scintillator has been registered.')
            return None
        return self._effmodel(KE,Mparticle,Zparticle,factor)**n

    def eff_extrap_beta(self, factorvals, kBvals, verbose=0):
        '''
        Calculate the apparent activity for an efficiency extrapolation measurement,
        for various kB values and TDCRs.
        
        Note that the *true* kB value is that given in the registered scintillator.
        Arguments:
        
        - A list of values to multiply the detection efficiency by.

        - A list of kB values (cm/MeV)
        
        The beta spectrum and scintillator will be inferred from the
        currently-registered objects.
        '''
        if self.scint is None:
            print('ERROR!! No scintillator has been registered.')
            return None
        elif self.branch is None:
            print('ERROR!! No branch has been registered.')
            return None
        kBtrue = self.scint.kB
        #Calculate the detection efficiency vs factor for the true kB.
        #Array stores efficiency[#PMTs][factor]
        effn_true = [[self.eff_nPMT_beta(n+1,factor,verbose-1) for factor in factorvals]
                     for n in range(3)]
        TDCR_true = [effn_true[2][idx]/effn_true[1][idx] for idx in range(len(factorvals))]        
        #This bit needs careful thought.
        #The TDCR is observed, so is independent of the estimated kB.
        #So the problem is to find effn_true[TDCR_wrong], as this is what the
        #efficiency actually is, as opposed to what we calculate with a
        #wrong kB.
        #The ratio effn_wrong[TDCR_wrong]/effn_true[TDCR_wrong] is therefore
        #proportional to the apparently measured activity.
        activity_meas = []
        TDCR_meas = []
        for thiskB in kBvals:
            if verbose>0:
                print('TDCR.eff_extrap_beta: kB = {:0.3f}'.format(thiskB))
            self.scint.setkB(thiskB)
            effn_wrong = [[self.eff_nPMT_beta(n+1,factor,verbose-1) for factor in factorvals]
                          for n in range(3)]
            TDCR_wrong = [effn_wrong[2][idx]/effn_wrong[1][idx] for idx in range(len(factorvals))]
            activity_meas += [[effn_wrong[1][idx]/np.interp(TDCR_wrong[idx],
                              TDCR_true,effn_true[1][:], left=np.inf,right=np.inf) 
                              for idx in range(len(factorvals))]]
            TDCR_meas += [TDCR_wrong]
        #reset kB to its original value
        self.scint.setkB(kBtrue)        
        return activity_meas, TDCR_meas, TDCR_true, effn_true
