# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 18:05:51 2019

@author: Cota
"""


from Read_IrisMeet import FileIRISmeet
from Scan2Signal import Scan2Sig
import Supporting_Functions as FUn

import time 
import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import minimize, least_squares, leastsq, curve_fit
from matplotlib.backends.backend_pdf import PdfPages
from scipy.special import erf
from scipy.interpolate import interp1d



######################################################################
############       IMPORT EXPERIMENTAL DATA       ####################
######################################################################

cd = 'Data/'

files = ['20190306_NaOH_1000mM.dat',
         '20190306_NaOH_1500mM.dat',
         '20190306_NaOH_2000mM.dat',
         '20190306_NaOH_2500mM.dat',
         '20190306_NaOH_3000mM.dat',
         '20190306_NaOH_4000mM.dat',
         '20190306_NaOH_5000mM.dat']


Time = []
Con = []
Par, ParE, Per, PerE, Iso, IsoE = [],[],[],[],[],[]

for sample in files:
  
    concen = (sample.split('mM')[0]).split('_')
    concen = (float(concen[len(concen)-1]))/1000.0
    
    print ('\n')    
    print ('{}\nConcentration: {} M\n{}'.format('-'*61,concen,'-'*20))    
    
    ###Read the file given by the IrisMeet software
    dat = FileIRISmeet(cd + sample)
    
    Sx = Scan2Sig(dat,0)                                                                                  
                                                                          
    ###############################################################################        
    ####################      VISUALIZING THE RAW DATA      #######################
    
    w0 = 2510
    t_plot = np.asarray([-10, 0.45, 0.5, 0.6, 0.8, 1, 1.2, 1.5, 2, 3, 4, 8, 15, 50, 100])
    
    t_nod = FUn.Plot_Raw(Sx, w0, t_plot, 0)
    ###Shift the time to t where A_iso(w0,t)/A_iso^max(w0, t) = - 0.7  
    Sx.t = FUn.ShiftT(Sx,t_nod)
    ###Delete pixels 
    Sx.v, Sx.S1, Sx.S1e, Sx.S2, Sx.S2e = FUn.delete_pixel(Sx, [0,1,18,29,30,31]) 
    ###Show raw data
    FUn.Plot_Raw(Sx, w0, t_plot, 0)
    ###Calculate the isotropic transient absorption
    Ax_iso, Ax_isoe = FUn.Free_Rot(Sx) 
        


    ###############################################################################        
    ########      DEFINING THE THE SPECTRAL AND TIME RANGE ANALYSIS      ########## 
    #########################      FOR THE ANALYSIS      ##########################    
    ###############################################################################        
      
    t_indx = FUn.t2pos(Sx.t,0.4)
    t_indxF = FUn.t2pos(Sx.t,75.0)    
        
        
    f_indx_max = FUn.v2p(Sx,2590)
    f_indx_min = FUn.v2p(Sx,2430)

    ###############################################################################        
    #######################       RAW DATA      ###################################      
    ###############################################################################        
    np.savetxt('Matrices/RawFrequency_%s.dat' % sample[9:-4], Sx.v, fmt='%s', delimiter='\t')    
    
    np.savetxt('Matrices/RawIsoExp_%s.dat' % sample[9:-4], Ax_iso[t_indx:t_indxF,:], fmt='%s', delimiter='\t')    
    np.savetxt('Matrices/RawIsoExpStDev_%s.dat' % sample[9:-4], Ax_isoe[t_indx:t_indxF,:], fmt='%s', delimiter='\t')    
    ###############################################################################        
    ###############################################################################        
        
    
    
    
    ###############################################################################        
    ##################        CUTTING UP DATA MATRICES        #####################      
    ###############################################################################        


    t = Sx.t[t_indx:t_indxF]
    np.savetxt('Matrices/RawTime_%s.dat' % sample[9:-4], t, fmt='%s', delimiter='\t')    
    
    freq = Sx.v[f_indx_max:f_indx_min + 1]
    np.savetxt('Matrices/Frequency_%s.dat' % sample[9:-4], freq, fmt='%s', delimiter='\t')    
    
   
    print ('\nThe analysis is set to the spectral range of %s - %s cm^1\nstarting from %s ps'% (int(min(freq)),int(max(freq)), round(t[0],2)))
    
    Para, Para_error = Sx.S1[t_indx:t_indxF,f_indx_max:f_indx_min + 1], Sx.S1e[t_indx:t_indxF,f_indx_max:f_indx_min + 1]
    Perp, Perp_error = Sx.S2[t_indx:t_indxF,f_indx_max:f_indx_min + 1], Sx.S2e[t_indx:t_indxF,f_indx_max:f_indx_min + 1]
    
    Iso_signal = Ax_iso[t_indx:t_indxF,f_indx_max:f_indx_min + 1]
    Iso_error = Ax_isoe[t_indx:t_indxF,f_indx_max:f_indx_min + 1]
    
    ###############################################################################    
    Del = [41, 42, 43, 44, 45, 46, 47, 48, 52]
    
    Para = np.delete(Para, Del, axis = 0)
    Para_error = np.delete(Para_error, Del, axis = 0)
    
    Perp = np.delete(Perp, Del, axis = 0)
    Perp_error = np.delete(Perp_error, Del, axis = 0)
    
    Iso_signal = np.delete(Iso_signal, Del, axis = 0)
    Iso_error =np.delete(Iso_error, Del, axis = 0)
    
    
    
    
    t = np.delete(t, Del)
    np.savetxt('Matrices/Time_%s.dat' % sample[9:-4], t, fmt='%s', delimiter='\t')    
    
    ###############################################################################   
    Par.append(Para)
    ParE.append(Para_error)
    
    Per.append(Perp)
    PerE.append(Perp_error)
    
    Iso.append(Iso_signal)
    IsoE.append(Iso_error)
    
    Time.append(t)
    
    Con.append(concen)
    ###############################################################################


#%%
###############################################################################        
###########            SINGULAR VALUE DECOMPOSITION               #############
###############################################################################        

# In some cases a SVD analysis can help to understand the undelying structure
# of the data.

print ('\n')
print ('Singular value decomposition (in %)\n{}'.format('-'*61))

print ('{:11}|{:11}|{:11}|{:11}|{:11}\n{}'.format('Conc', '1st Comp', '2nd Comp', '3rd Comp', '4th Comp', '-'*61))

l=0
for Iso_signal in Iso:
    u, s, vh = np.linalg.svd(Iso_signal, full_matrices=True, compute_uv=True)
    print( '{:11}|{:11.2f}|{:11.2f}|{:11.2f}|{:11.2f}'.format(Con[l],*((s/sum(s))[:4]*100)))
    l+=1
print ('{}'.format('-'*61))

# >99% of the structure is well represented with the first 3 components
# for all the samples under study here.

###############################################################################        



#%%
###############################################################################
###############################################################################
# Creates an external file to count the number of iterations 
# cycles to show the evolution of the fit
Counts = [1,10]
Counts = np.array(Counts, dtype='int32')
np.savetxt('Counts.dat', Counts, fmt='%d', delimiter='\t')


# Calculates the X2 from the difference between the experimental data
# and the predicted values N@Spec.
def Minimize_Iso(P,Spec,dat,datE):
    N = np.array(P)
    diff = np.power((dat-N@Spec)/datE,2.0)
    return diff


# In this model we assume that the spectra are conformed by the linear 
# combination of three spectral shapes: two originally excited states 
# and one thermal state.

# The shapes of these spectral shapes are assumed concentration independent
# and only the population matrix tell us the amplitude of each and every 
# time point. In pther words, the population matrix contains the information 
# of concentration.

def Global_Analysis(paramFit, concentration, Time, freq, Iso, IsoE):
    
    iteration, cycle = np.loadtxt('Counts.dat',dtype='int32')
     
     
    # Based on experiments of NaCl in pure H2O (no D2O), I measured
    # the shape of the ion-associated excited state (called Background here)
     
    # The first of the following two lines cointains the results of such
    # pure H2O experiments, thus teh backgroud signal can be used here.
    # The second line contains the results that were saved in a previous 
    # analysis. 
    # They will give exactly the same result, but the first will take 
    # something like 10 times longer to converge.
     
#    SignB = np.loadtxt('Matrices/Signature_Background.dat', unpack=True)[1]
    SignB = (np.loadtxt('Matrices/Signatures_NaOH_Solved.dat', unpack=True))[1]


    GlobalResidual = []   
    
    for idx in range(len(concentration)):
        
        Iso_signal = Iso[idx]
        Iso_error = IsoE[idx]
        
        t = Time[idx]
        
        ###############################################################################   
        
        # Calculate the spectral shape of the thermal signal from the average
        # of the last four observed transient absorption signal                
        Weights = 1.0/(Iso_error*Iso_error)
        Iso_Heat = np.average(Iso_signal[(len(t) - 4):,:], axis = 0, weights = Weights[(len(t) - 4):,:])
                
        # Create a matrix with dimensions len(spectral range)x3 that contains
        # the shapes of the three spectral components
        SolveSpecs = np.vstack((paramFit,SignB,Iso_Heat))
        
        # Save the spectral shapes. Note: the spectral shape of the two
        # excited components are concentration independent. The third is not,
        # but very well defined at long delay times for each concentration.
        np.savetxt('Matrices/Signatures_NaOH_%smM.dat' % int(1000*concentration[idx]), SolveSpecs.T, fmt='%.18e', delimiter='\t')
        
        # Create the population matrix
        PopM = np.zeros((len(t),len(SolveSpecs)),dtype='float128')
        PopM_StDev = np.zeros((len(t),len(SolveSpecs)),dtype='float128')

        # Calculates the extend to which each component contributes 
        # to the isotropic signal at each delay time t
        p0 = np.array([0.5,0.5,0.0])        
        for i in range (len(t)):
            params = least_squares(Minimize_Iso, p0, 
                                   args=(SolveSpecs, Iso_signal[i,:], Iso_error[i,:]), 
                                   bounds=(0.0, 1.8),  method='trf', loss='linear', max_nfev = 1000)
            
            PopM[i,:] =  params['x']
            
            PopM_StDev[i,:] = np.sqrt((sum(params['fun']))*np.diag(np.linalg.inv((params['jac'].T)@params['jac'])))
            
            p0 = params['x']

        # Save the population matric for each concentration
        np.savetxt('Matrices/Population_NaOH_%smM.dat' % int(1000*concentration[idx]), np.insert(PopM,0,t,axis=1), fmt='%.18e', delimiter='\t')
        np.savetxt('Matrices/Population_StDev_NaOH_%smM.dat' % int(1000*concentration[idx]), np.insert(PopM_StDev,0,t,axis=1), fmt='%.18e', delimiter='\t')        
        
        # Save the calculated isotropic transient signal for each concentration
        np.savetxt('Matrices/IsoFit_NaOH_%smM.dat' % int(1000*concentration[idx]), PopM@SolveSpecs, fmt='%s', delimiter='\t')    
             
        # Print the evolution of the global fit every 10 cycles.
        if cycle == 10:
            
            if idx == 0:
                plt.plot(freq,paramFit,color='b',label='Bulk')
                plt.plot(freq,SignB,color='g',label='Ion associated')
#                plt.plot(freq,SignW,color='b')
#                plt.plot(freq,paramFit,color='g')
                plt.plot(freq,Iso_Heat,color='r',label='Thermal')
                plt.text(2530,-0.018,'Iteration = %d' % iteration)
                plt.title('Specral components')
                plt.legend(loc=3)
                plt.show()
                plt.close()
            
            for l in range(10):
                colors = plt.cm.gnuplot(np.linspace(0,1,11))
                n = [0,2,4,8,12,16,22,30,38,50]
                plt.plot(freq,Iso_signal[n[l],:],'.',color=colors[l])
                plt.plot(freq, (PopM@SolveSpecs)[n[l],:], color = colors[l])
            plt.title('Isotropic transient absorption - Concentration %s M' % Con[idx], fontsize=15)
            plt.show()
            plt.close() 
                     
        
        Res = np.power((Iso_signal - PopM@SolveSpecs)/Iso_error,2.0)
         
        GlobalResidual.extend(Res)  
    

    
    if cycle==10:
        cycle = 0
            
    iteration += 1
    cycle += 1    
    
    np.savetxt('Counts.dat', np.array([iteration,cycle]), fmt='%d', delimiter='\t')
    
    GlobalResidual = np.array(GlobalResidual)
        
    return GlobalResidual.reshape(-1)




# I let the spectral shape of the excited molecules in the bulk
# change along the fit. Yet, the calculated shape is shared in
# all the ion concentrations

# In the following two line, you can extract the initial guess for 
# for the spectral shape of the excited molecules in the bulk.
# The first contains the spectral shape from a sample of pure water
# The second is the solution that was calculated here in a previous run
# of the program. Both will give the same results, but the first will 
# take at least 10 times longer. 

#SignW = (np.loadtxt('Matrices/NeatHDO_Signature.dat', unpack=False))[0]
SignW = (np.loadtxt('Matrices/Signatures_NaOH_Solved.dat', unpack=True))[0]



# Calculated the spectral shape of the excited molecules in the bulk
# that best fit the experimental data.
time1 = time.time()
GValues = least_squares(Global_Analysis, SignW, args=(Con, Time, freq, Iso, IsoE), bounds=(-np.inf,0.0),  method='trf', loss='linear', max_nfev = 2500)
time2 = time.time()
print ('Processing time: ', time2 - time1)


#Residual = sum(Global_Analysis(SignW, Con, Time, freq, Iso, IsoE))



#%%
# Compare the final result with the spectral shape of the excited molecules 
# in the bulk from a sample of pure water
SignWater = GValues['x']
SignWater_StDev = np.sqrt((sum(GValues['fun']))*np.diag(np.linalg.inv((GValues['jac'].T)@GValues['jac'])))

SignNeat = (np.loadtxt('Matrices/NeatHDO_Signature.dat', unpack=False))

plt.errorbar(freq,SignWater,yerr=SignWater_StDev)
plt.plot(freq,SignNeat)
plt.show()
plt.close()
#
#

#%%


# Forster survival probability.
# Reference: https://doi.org/10.1063/1.3432616
def Forster_Survival(t,Molar,kf,a0):
    Gra = np.exp(-(4.0*np.pi*0.0006022*Molar/3.0)*(np.power(a0,3.0)*(np.exp(-(kf*t)/np.power(a0,6.0))-1.0) + np.sqrt(kf)*np.sqrt(np.pi)*np.sqrt(t)*erf(np.sqrt(kf*t)/np.power(a0,3.0))))
#    Gra = (4.0*np.pi*Rho/3.0)*(np.sqrt(kf)*np.sqrt(np.pi)*np.sqrt(t)*erf(np.sqrt(kf*t)/np.power(a0,3.0)))
    return Gra


# Fit the decay of the excited stated as the convolution of an exponential 
# decay and a Forster energy transfer process
def residual_Forster(param,conc,Surv):
    
    T1, kf = param[0], param[1]
    Scalar = param[2:]
    
    Res = []    
    
    for i in range(len(conc)):
        t, popW = np.loadtxt('Matrices/Population_NaOH_%smM.dat' % int(1000*conc[i]), unpack=True)[:2]
        
        diff = np.power(popW*Scalar[i] - np.exp(-t/T1)*Forster_Survival(t,conc[i],kf,1.5),2.0)
        
        Res.extend(diff)
        
    return Res


# Initial values for the rate of the exponential decay and the rate of 
# the Forster energy tranfer process ([0,1]-elements). The extra elements
# are used to normalize the decay rate to 1 at t=0.
Par0 = np.array([1.8,300.0,1.82,2.1,2.6,2.95,3.1,5.4,6.24])    
    
Popt = least_squares(residual_Forster, Par0, args=(Con, Forster_Survival), bounds=(0.0,np.inf),  method='trf', loss='linear', max_nfev = 20000)

#print (Popt['x'])

ValOpt = Popt['x']

T1_opt, rf_opt = ValOpt[0], ValOpt[1]
normScalar = ValOpt[2:]


with PdfPages('Plots/Survival_DelayCurves.pdf') as pdf:

    fig, ax1 = plt.subplots(1, 1, figsize=(9,6))

    colors = plt.cm.gnuplot(np.linspace(0,1,8))
    
    for i in range (7):
        t, popW, popLoc, popH = np.loadtxt('Matrices/Population_NaOH_%smM.dat' % int(1000*Con[i]), unpack=True)
        tz, popW_StDev, popLoc_StDev, popH_StDev = np.loadtxt('Matrices/Population_StDev_NaOH_%smM.dat' % int(1000*Con[i]), unpack=True)

        Timet = np.arange(0,10,0.05)

        ax1.plot(Timet, np.exp(-Timet/T1_opt)*Forster_Survival(Timet,Con[i],rf_opt,1.5), lw=1.5, color = colors[i], label = r'$\mathrm{%s\/\/M}$' % Con[i])

        ax1.plot(t, normScalar[i]*popW, '.', lw = 2.0, color = colors[i])
        
        plt.errorbar(t, normScalar[i]*popW, yerr=normScalar[i]*popW_StDev, color=colors[i], fmt='.')#,label='isotropic signal '+str(t[i]))
        
        
#    ax1.plot(t,np.exp(-t/1.75),'k--')

    ax1.set_xlabel(r'$\mathrm{Time}$ $\mathrm{delay}$ $\mathrm{[ps]}$', size = 24)
    ax1.set_ylabel(r'$\mathrm{Population}$  $\mathrm{\mathbf{OD}} \cdots \mathrm{H_2O}$',size=24)

#    plt.yscale('log')
    ax1.set_xlim(0.0,5.0)
    ax1.set_ylim(0.02,1.0)
    ax1.tick_params(labelsize=16)

    ax1.legend(loc=1,fontsize=20)
    
    
    pdf.savefig(bbox_inches='tight')
    plt.show()
    plt.close()



#%%
# Some statistical analysis
# Import the data for the Forster distribution function (FDF) from a numerical solution
# using Wolfram Mathematica of equation 15 in Ref:https://doi.org/10.1063/1.3432616.
# The FDF was calculated both for hydroxide ions and for protons.
Time, ForsterD = np.loadtxt('ForsterDistance.dat', delimiter = ',', unpack=True)

# Separate the data for hydroxide ions and protons
Time1, Time2 = Time[:443], Time[443:]
ForsterD11, ForsterD22 = ForsterD[:443], ForsterD[443:]

Time1 = Time1[35:] # Hydroxide ions
Time2 = Time2[39:] # Protons

ForsterD1 = ForsterD11[35:] # Hydroxide ions
ForsterD2 = ForsterD22[39:] # Protons


f = interp1d(Time1, ForsterD1)
Distance = np.arange(1.514,49.9,0.005)
DisForster1 = f(Distance)
     
ForD1 = sum(DisForster1*Distance)/sum(DisForster1)    

print ('\n')     
print ('Mediam (hydroxide ions) = ', ForD1)
    
       
f2 = interp1d(Time2, ForsterD2)
Distance = np.arange(1.514,49.9,0.005)
DisForster2 = f2(Distance)
    
ForD2 = sum(DisForster2*Distance)/sum(DisForster2)    

print ('Median (protons:Rutger result) = ', ForD2)
print ('\n')     
 


# Plot the FDF and show the median
with PdfPages('Plots/ForsterDistribution.pdf') as pdf:

    fig, ax1 = plt.subplots(1, 1, figsize=(9,6))

    colors = plt.cm.gnuplot(np.linspace(0,1,8))

    ax1.axvline(x=ForD1, color = 'r', ls= '--')
    ax1.axvline(x=ForD2, color = 'b', ls= '--')
    
    ax1.plot(Time1, ForsterD1/sum(ForsterD1), lw=1.5, color = 'r', label = r'$\mathrm{OH^-}$')
    ax1.plot(Time2, ForsterD2/sum(ForsterD2), lw=1.5, color = 'b', label = r'$\mathrm{H^+}$')

    ax1.set_xlabel(r'$\mathrm{Distance}$ $\mathrm{[\AA]}$', size = 24)
    ax1.set_ylabel(r'$\mathrm{F\"orster}$ $\mathrm{transfer}$ $\mathrm{distribution}$',size=24)

    ax1.set_xlim(0.0,10.0)
    
    ax1.tick_params(labelsize=16)
    
    ax1.legend(loc=1,fontsize=20)
    
    pdf.savefig(bbox_inches='tight')
    plt.show()
    plt.close()

