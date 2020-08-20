# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 11:11:57 2019

@author: Cota
"""

import Read_IrisMeet as ReadIris
import Scan2Signal as S2S
import FunctionsT as FUn
#import Models as Model
#import Signatures as Sign
#import Iso_Functions as IsoFu
#import AnisotropyDecomposition as ADeco

import time 
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize, least_squares, leastsq, curve_fit
#from scipy import linalg
from sympy import symbols, Matrix, lambdify, exp
from matplotlib.backends.backend_pdf import PdfPages
from scipy.special import erf
from scipy.interpolate import interp1d
#from scipy.interpolate import interp1d

linewide=2


###############################################################################
########################    EXPERIMENTAL DATA    ##############################
###############################################################################

#cd = ('/Users/Cota/Documents/PhD/UltraFast/Anisotropy_Meas/20190306/')

cd = ('/Users/Cota/Documents/PhD/UltraFast/ForsterTransferOH/NaOD/')

files = ['20190601_NaOD_1M_1.dat',
         '20190601_NaOD_2M_1.dat',
         '20190601_NaOD_3M_1.dat',
         '20190601_NaOD_4500mM_1.dat',
         '20190601_NaOD_6M_1.dat']


ScanG = [0,
         0,
         [0,2,3,4,6],
         [0,1,2,3],
         [1,2,3,4,5,6,7,8,10,12,13,14,15,16,17,18,19],
         [0,2,3,4,5,6],
         0]
#         [0,1,2,4,5,6,7,8,9,10]]


Sampler = ['1.0 M','2.0 M','3.0 M','4.5 M','6.0 M']


Par, ParE, Per, PerE, Iso, IsoE = [],[],[],[],[],[]
Time = []
Con = [1.0, 2.0, 3.0, 4.5, 6.0]

l=0
for sample in files:
    
#    concen = (sample.split('mM')[0]).split('_')
#    concen = (float(concen[len(concen)-1]))/1000.0
    
    ###Read the file given by the IrisMeet software
    dat = ReadIris.Read(sample)
    
    Sx = S2S.Scan_to_Signal(dat,ScanG[l])  

    ###Swap polarization
    SigA, ErrA, SigB, ErrB = Sx.S1, Sx.S1e, Sx.S2, Sx.S2e 
    Sx.S1, Sx.S1e, Sx.S2, Sx.S2e = SigB, ErrB, SigA, ErrA
                                                                                
                                                                          
    ###############################################################################        
    ####################      VISUALIZING THE RAW DATA      #######################
    
    w0 = 3380
    t_plot = np.asarray([-10, 0.45, 0.5, 0.6, 0.8, 1, 1.2, 1.5, 2, 3, 4, 8, 15, 50, 100])
    
    t_nod = FUn.Plot_Raw(Sx, w0, t_plot, 0)
    ###Shift the time to t where A_iso(w0,t)/A_iso^max(w0, t) = - 0.7  
    Sx.t = FUn.ShiftT(Sx,t_nod)
    ###Delete pixels 
    Sx.v, Sx.S1, Sx.S1e, Sx.S2, Sx.S2e = FUn.delete_pixel(Sx, [0,1,2,3,4,5,6,7,8,9,28,29,30,31]) 
    ###Show raw data
    FUn.Plot_Raw(Sx, w0, t_plot, 1)
    ###Calculate the isotropic transient absorption
    Ax_iso, Ax_isoe = FUn.Free_Rot(Sx) 
        


    ###############################################################################        
    ########      DEFINING THE THE SPECTRAL AND TIME RANGE ANALYSIS      ########## 
    #########################      FOR THE ANALYSIS      ##########################    
    ###############################################################################        
      
    t_indx = FUn.t2pos(Sx.t,0.44)
    t_indxF = FUn.t2pos(Sx.t,30.5)+1    
        
        
    t = Sx.t[t_indx:t_indxF]
    np.savetxt('Matrices/Time_%smM.dat' % int(1000*Con[l]), t, fmt='%s', delimiter='\t')    
        
    ###############################################################################        
    #######################       RAW DATA      ###################################      
    ###############################################################################        
    np.savetxt('Matrices/RawFrequency_%smM.dat' % int(1000*Con[l]), Sx.v, fmt='%s', delimiter='\t')    
    
    np.savetxt('Matrices/RawIsoExp_%smM.dat' % int(1000*Con[l]), Ax_iso[t_indx:t_indxF,:], fmt='%s', delimiter='\t')    
    np.savetxt('Matrices/RawIsoExpStDev_%smM.dat' % int(1000*Con[l]), Ax_isoe[t_indx:t_indxF,:], fmt='%s', delimiter='\t')    
    ###############################################################################        
    ###############################################################################        
        

    ###############################################################################        
    ##################        CUTTING UP DATA MATRICES        #####################      
    ###############################################################################        
    
    f_indx_max = FUn.v2p(Sx,3600)
    f_indx_min = FUn.v2p(Sx,3250)
    
    freq = Sx.v[f_indx_max:f_indx_min + 1]
    np.savetxt('Matrices/Frequency_%smM.dat' % int(1000*Con[l]), freq, fmt='%s', delimiter='\t')    
    
   
    print ('The analysis is made in the spectral range of %s - %s cm^1\nstarting from %s ps\n'% (int(min(freq)),int(max(freq)), round(t[0],2)))
    
    Para, Para_error = Sx.S1[t_indx:t_indxF,f_indx_max:f_indx_min + 1], Sx.S1e[t_indx:t_indxF,f_indx_max:f_indx_min + 1]
    Perp, Perp_error = Sx.S2[t_indx:t_indxF,f_indx_max:f_indx_min + 1], Sx.S2e[t_indx:t_indxF,f_indx_max:f_indx_min + 1]
    
    Iso_signal = Ax_iso[t_indx:t_indxF,f_indx_max:f_indx_min + 1]
    Iso_error = Ax_isoe[t_indx:t_indxF,f_indx_max:f_indx_min + 1]
    
#    ###############################################################################    
#    Del = [41, 42, 43, 44, 45, 46, 47, 48, 52]
#    
#    Para = np.delete(Para, Del, axis = 0)
#    Para_error = np.delete(Para_error, Del, axis = 0)
#    
#    Perp = np.delete(Perp, Del, axis = 0)
#    Perp_error = np.delete(Perp_error, Del, axis = 0)
#    
#    Iso_signal = np.delete(Iso_signal, Del, axis = 0)
#    Iso_error =np.delete(Iso_error, Del, axis = 0)
    
    
    
#    u, s, vh = np.linalg.svd(Iso_signal, full_matrices=True, compute_uv=True)
#    
#    print (100*(s[:4]/sum(s)))
    
#    A = np.arange(14)
#    
#    
#    plt.figure(figsize=(10,6))
#    for i in range (5):
#        plt.plot(freq, vh[i], ls = '-', lw = 2, label=r'$\mathrm{Level}$ %s' % (i+1))
#    
#    plt.legend(loc=1)
#    #plt.grid(True)
#    plt.show()
#    plt.close()
#    
#    
#    B = np.arange(len(t))
#    
#    plt.figure(figsize=(10,6))
#    for i in range (5):
#        plt.plot(B, u[:,i], ls = '-', lw = 2, label=r'$\mathrm{Level}$ %s' % (i+1))
#    
#    plt.legend(loc=1)
#    #plt.grid(True)
#    plt.show()
#    plt.close()
    
    
    
    
#    t = np.delete(t, Del)
#    np.savetxt('Matrices/Time_%smM.dat' % int(1000*Con[l]), t, fmt='%s', delimiter='\t')    
    
    ###############################################################################   
    Par.append(Para)
    ParE.append(Para_error)
    
    Per.append(Perp)
    PerE.append(Perp_error)
    
    Iso.append(Iso_signal)
    IsoE.append(Iso_error)
    
    Time.append(t)
    
    l+=1
    
#    Con.append(concen)
    ###############################################################################

###############################################################################
###############################################################################
###############################################################################



#Counts = [1,10]
#Counts = np.array(Counts, dtype='int32')
#np.savetxt('Counts.dat', Counts, fmt='%d', delimiter='\t')

#%%

def Survival(t,Molar,kf,a0):
    Gra = np.exp(-(4.0*np.pi*0.0006022*Molar/3.0)*(np.power(a0,3.0)*(np.exp(-(kf*t)/np.power(a0,6.0))-1.0) + np.sqrt(kf)*np.sqrt(np.pi)*np.sqrt(t)*erf(np.sqrt(kf*t)/np.power(a0,3.0))))
#    Gra = (4.0*np.pi*Rho/3.0)*(np.sqrt(kf)*np.sqrt(np.pi)*np.sqrt(t)*erf(np.sqrt(kf*t)/np.power(a0,3.0)))
    return Gra



def Global_Analysis(ParFit, Conchen, Time, freq, Iso, IsoE):
    
    def Minimo(P,Spec,Iso,IsoE):
        N = np.array(P)
        diff = np.power((Iso-N@Spec)/IsoE,2.0)
        return diff
        
#    iteration, cycle = np.loadtxt('Counts.dat',dtype='int32')

     
#    SignW = (np.loadtxt('Matrices/NeatHDO_Signature.dat', unpack=False))
#    SignW = (np.loadtxt('Matrices/Signatures_NaOH_Solved.dat', unpack=True))[0]
    Sign1 = np.loadtxt('/Users/Cota/Documents/PhD/UltraFast/Anisotropy_Meas/20190601/Matrices/SpectralDecomposition_1000mM.dat', unpack=True)[1]
    Sign2 = np.loadtxt('/Users/Cota/Documents/PhD/UltraFast/Anisotropy_Meas/20190601/Matrices/SpectralDecomposition_6000mM.dat', unpack=True)[2]
#    SignB = (np.loadtxt('Matrices/Signatures_NaOH_Solved.dat', unpack=True))[1]



    GlobalResidual = []   
    
    Scalar = [1.6,1.9,2.2,2.5,3.7]
    
    for idx in range(len(Conchen)):
        
        Iso_signal = Iso[idx]#.astype('float128')
        Iso_error = IsoE[idx]#.astype('float128')
        
        t = Time[idx]
        
        ###############################################################################   
        ###############################################################################   
        ###############################################################################   
                         
        Weights = 1.0/(Iso_error*Iso_error)
        Iso_Heat = np.average(Iso_signal[(len(t) - 4):,:], axis = 0, weights = Weights[(len(t) - 4):,:])
                
        
        SolveSpecs = np.vstack((Sign1,Sign2,Iso_Heat))
        
        np.savetxt('Matrices/Signatures_NaOD_%smM.dat' % int(1000*Conchen[idx]), SolveSpecs.T, fmt='%.18e', delimiter='\t')
        
        PopM = np.zeros((len(t),len(SolveSpecs)),dtype='float128')
        PopM_StDev = np.zeros((len(t),len(SolveSpecs)),dtype='float128')

        p0 = np.array([0.5,0.5,0.0])        
        
        for i in range (len(t)):
            params = least_squares(Minimo, p0, args=(SolveSpecs, Iso_signal[i,:], Iso_error[i,:]), bounds=(0.0, 1.8),  method='trf', loss='linear', max_nfev = 1000)
            
            PopM[i,:] =  params['x']
            
            PopM_StDev[i,:] = np.sqrt((sum(params['fun']))*np.diag(np.linalg.inv((params['jac'].T)@params['jac'])))
            
            p0 = params['x']


        np.savetxt('Matrices/Population_NaOD_%smM.dat' % int(1000*Conchen[idx]), np.insert(PopM,0,t,axis=1), fmt='%.18e', delimiter='\t')
        np.savetxt('Matrices/Population_StDev_NaOD_%smM.dat' % int(1000*Conchen[idx]), np.insert(PopM_StDev,0,t,axis=1), fmt='%.18e', delimiter='\t')        
        
        np.savetxt('Matrices/IsoFit_NaOD_%smM.dat' % int(1000*Conchen[idx]), PopM@SolveSpecs, fmt='%s', delimiter='\t')    
        
        
        plt.plot(t,Scalar[idx]*PopM[:,0])
        plt.errorbar(t,Scalar[idx]*PopM[:,0],yerr=Scalar[idx]*PopM_StDev[:,0])
#        plt.plot(t,PopM[:,1])
#        plt.plot(t,PopM[:,2])
        tt=np.arange(0,7,0.02)
        plt.plot(tt,np.exp(-tt/0.78))
        plt.xlim(0,3.5)
        plt.ylim(0.01,1)
        plt.yscale('log')
        plt.show()
        plt.close()        
        

        
#    if cycle == 10:
        
        if idx == 0:
            plt.plot(freq,Sign1,color='b')
            plt.plot(freq,Sign2,color='g')

#                plt.plot(freq,SignW,color='b')
#                plt.plot(freq,ParFit,color='g')

            plt.plot(freq,Iso_Heat,color='r')
#            plt.text(2520,-0.015,'%d' % iteration)
#            plt.show()
            plt.close()
        
        
        n = [16,7,4,2,1,0.42]
        for l in range(len(n)):
            colors = plt.cm.gnuplot(np.linspace(0,1,11))
            
            idxt = (np.abs(t - n[l])).argmin()

            plt.plot(freq,Iso_signal[idxt,:],'.',color=colors[l])
            plt.plot(freq, (PopM@SolveSpecs)[idxt,:], color = colors[l])
            
        plt.text(3330,-0.017,'%s' % idx, fontsize=15)
#        plt.show()
        plt.close() 
                     
        
        Res = np.power((Iso_signal - PopM@SolveSpecs)/Iso_error,2.0)
         
        GlobalResidual.extend(Res)  
    

    
#    if cycle==10:
#        cycle = 0
#            
#    iteration += 1
#    cycle += 1
#    
##    print (iteration, '\t', cycle)        
#    
#    
#    np.savetxt('Counts.dat', np.array([iteration,cycle]), fmt='%d', delimiter='\t')
    
    
    GlobalResidual = np.array(GlobalResidual)
        
    return GlobalResidual.reshape(-1)








#SignW = (np.loadtxt('Matrices/NeatHDO_Signature.dat', unpack=False))
#SignW = (np.loadtxt('Matrices/Signatures_NaOH_Solved.dat', unpack=True))[0]

#SignB = np.loadtxt('/Users/Cota/Documents/PhD/UltraFast/Anisotropy_Meas/20190306/Background/Signature_Background.dat', unpack=True)[1]



p0 = -58.0


time1 = time.time()
GValues = least_squares(Global_Analysis, p0, args=(Con, Time, freq, Iso, IsoE), bounds=(-np.inf,0.0),  method='trf', loss='linear', max_nfev = 2500)
time2 = time.time()
print ('Processing time: ', time2 - time1)


#Residual = sum(Global_Analysis(SignW, Con, Time, freq, Iso, IsoE))


#print (stop)


#%%


## Population Normalization
Scalar = [1.5,1.75,2.2,3.5,4.1]


with PdfPages('Survival_DelayCurves.pdf') as pdf:

    fig, ax1 = plt.subplots(1, 1, figsize=(9,6))

    colors = plt.cm.gnuplot(np.linspace(0,1,7))
#    colors = colors[::-1]
    
    for i in range (5):
        t, popW, popIon, popH = np.loadtxt(cd + 'Matrices/Population_NaOD_%smM.dat' % int(1000*Con[i]), unpack=True)
        tz, popW_StDev, popIon_StDev, popH_StDev = np.loadtxt(cd + 'Matrices/Population_StDev_NaOD_%smM.dat' % int(1000*Con[i]), unpack=True)

        Timet = np.arange(0,4,0.02)

#        ax1.plot(Timet, np.exp(-Timet/T1_opt)*Survival(Timet,Conc[i],rf_opt,1.5), lw=1.5, color = colors[i], label = r'$\mathrm{%s\/\/M}$' % Conc[i])
#
#        ax1.plot(t, scal_opt[i]*popW, '.', lw = 2.0, color = colors[i])
        
        plt.errorbar(t, Scalar[i]*popW, yerr=Scalar[i]*popW_StDev, color=colors[i+1], fmt='.',label='%s' % i)#,label='isotropic signal '+str(t[i]))
        
    tt=np.arange(0,4,0.02)
    ax1.plot(tt,np.exp(-tt/0.74),'k', lw=2)

    ax1.set_xlabel(r'$\mathrm{Time}$ $\mathrm{delay}$ $\mathrm{[ps]}$', size = 24)
    ax1.set_ylabel(r'$\mathrm{Population}$  $\mathrm{\mathbf{OD}} \cdots \mathrm{H_2O}$',size=24)

    plt.yscale('log')
    ax1.set_xlim(0.0,3.0)
    ax1.set_ylim(0.01,1.0)
    ax1.tick_params(labelsize=16)

    ax1.legend(loc=1,fontsize=20)
    
    
    pdf.savefig(bbox_inches='tight')
    plt.show()
    plt.close()

