# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 18:19:37 2019

@author: Cota
"""

#import Read_IrisMeet as ReadIris
#import Scan2Signal as S2S
#import FunctionsT as FUn
#import Models as Model
#import Signatures as Sign
#import Iso_Functions as IsoFu
#import AnisotropyDecomposition as ADeco

#import time 
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize, least_squares, leastsq, curve_fit
#from scipy import linalg
#from sympy import symbols, Matrix, lambdify, exp
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import OldScalarFormatter, ScalarFormatter, LogFormatter, NullFormatter, FuncFormatter
import matplotlib.ticker as mticker
#from scipy.interpolate import interp1d
from matplotlib.patches import Rectangle
from scipy.special import erf
from scipy.interpolate import interp1d

plt.style.use('classic')

#files = [   'Water.dpt',
#            'Phenol_840mM.dpt',
#            'Phenolate_840mM.dpt',
#            'Phenolate_1000mM.dpt',
#            'Phenolate_2000mM.dpt',
#            'Phenolate_3000mM.dpt']
#
#
#Names = [   '4% D/H',
#            '[PhOH] = 0.84 M', 
#            '[PhONa] = 0.84 M', 
#            '[PhONa] = 1 M',
#            '[PhONa] = 2 M',
#            '[PhONa] = 3 M']
            

########################################
##########        FTIR        ##########
########################################


def Frequency(file):
    data = open(file, 'r').read()
    freqs = data.split('\n')[2]
    Freqs = np.array(list(map(float, freqs.split()[1:])))
    return Freqs

Amarillo = '#ad5b50'
Rojo='#78b26d'
Purpura='#4986ae'



cd1 = '/Users/Cota/Documents/PhD/UltraFast/DifferentSalts/NaOH/H2O/FTIR/'

Abs_Spectrum01 = np.loadtxt(cd1 + 'FTIR_4pD2O_H2O_25um_001.dat',comments='#', unpack=True)
Abs_Spectrum02 = np.loadtxt(cd1 + 'FTIR_4pD2O_H2O_25um_002.dat',comments='#', unpack=True)

Abs_Spectrum11 = np.loadtxt(cd1 + 'FTIR_NaOH_1M_4pD2O_H2O_25um_001.dat',comments='#', unpack=True)
Abs_Spectrum12 = np.loadtxt(cd1 + 'FTIR_NaOH_1M_4pD2O_H2O_25um_002.dat',comments='#', unpack=True)

Abs_Spectrum21 = np.loadtxt(cd1 + 'FTIR_NaOH_2M_4pD2O_H2O_25um_001.dat',comments='#', unpack=True)
Abs_Spectrum22 = np.loadtxt(cd1 + 'FTIR_NaOH_2M_4pD2O_H2O_25um_002.dat',comments='#', unpack=True)


Spectrum0MH = (Abs_Spectrum01 + Abs_Spectrum02)/2.0
Spectrum15MH = (Abs_Spectrum11 + Abs_Spectrum12)/2.0
Spectrum3MH = (Abs_Spectrum21 + Abs_Spectrum22)/2.0






cd2 = '/Users/Cota/Documents/PhD/UltraFast/DifferentSalts/NaOH/D2O/FTIR/'

Abs_Spectrum01 = np.loadtxt(cd2 + 'FTIR_HeavyWater_4pH2O_25umSpacer_001.dat',comments='#', unpack=True)
Abs_Spectrum02 = np.loadtxt(cd2 + 'FTIR_HeavyWater_4pH2O_25umSpacer_002.dat',comments='#', unpack=True)

Abs_Spectrum11 = np.loadtxt(cd2 + 'FTIR_HeavyWater_4pH2O_NaOD_1M_25umSpacer_001.dat',comments='#', unpack=True)
Abs_Spectrum12 = np.loadtxt(cd2 + 'FTIR_HeavyWater_4pH2O_NaOD_1M_25umSpacer_002.dat',comments='#', unpack=True)

Abs_Spectrum21 = np.loadtxt(cd2 + 'FTIR_HeavyWater_4pH2O_NaOD_2M_25umSpacer_001.dat',comments='#', unpack=True)
Abs_Spectrum22 = np.loadtxt(cd2 + 'FTIR_HeavyWater_4pH2O_NaOD_2M_25umSpacer_002.dat',comments='#', unpack=True)

Abs_Spectrum41 = np.loadtxt(cd2 + 'FTIR_HeavyWater_4pH2O_NaOD_4M_25umSpacer_001.dat',comments='#', unpack=True)
Abs_Spectrum42 = np.loadtxt(cd2 + 'FTIR_HeavyWater_4pH2O_NaOD_4M_25umSpacer_002.dat',comments='#', unpack=True)


Spectrum0MD = (Abs_Spectrum01 + Abs_Spectrum02)/2.0
Spectrum15MD = Abs_Spectrum12
Spectrum3MD = (Abs_Spectrum21 + Abs_Spectrum22)/2.0
Spectrum4MD = (Abs_Spectrum41 + Abs_Spectrum42)/2.0







Freq1 = Frequency(cd1 + 'FTIR_4pD2O_H2O_25um_001.dat')

Freq2 = Frequency(cd2 + 'FTIR_HeavyWater_4pH2O_25umSpacer_001.dat')


#plt.style.use('default')




#with PdfPages('LinearSpectra.pdf') as pdf:
#    fig, [ax1,ax2] = plt.subplots(1, 2, figsize=(14,5))
##            fig.subplots_adjust(hspace=0.03)
#    
#    ax1.plot(Freq1, Spectrum3MH, ':', color = Rojo, lw=2, label = '3.0 M')
#    ax1.plot(Freq1, Spectrum15MH, '--', color = Purpura, lw=2, label = '1.5 M')
#    ax1.plot(Freq1, Spectrum0MH, 'k', lw=2, label = '0.0 M')
#    ax1.set_ylim(0.0,3.0)
#    ax1.set_xlim(2250, 2950)
#    ax1.legend(loc=1)
#
#    ax1.set_ylabel(r'$\mathrm{Absorbance}$', size = 18)
#    ax1.set_xlabel(r'$\mathrm{Frequency}$ $\mathrm{[cm^{-1}]}$', size = 18)
#
#   
##    left, bottom, width, height = 
#    ax11 = fig.add_axes([0.155, 0.58, 0.19, 0.29])    
#    
#    ax11.plot(Freq1, Spectrum3MH, ':', color = Rojo, lw=2)
#    ax11.plot(Freq1, Spectrum15MH, '--', color = Purpura, lw=2)
#    ax11.plot(Freq1, Spectrum0MH, 'k', lw=2)
#    ax11.yaxis.set_ticks(np.arange(0,4.1,1))
#    ax11.xaxis.set_ticks([1500,2500,3500])
#    ax11.set_xlim(1000,4000)
#    ax11.set_ylim(0.0,4.0)
#    
#    
#    ax2.plot(Freq2, 1.25*Spectrum4MD, ':', color = Amarillo, lw=2, label = '6.0 M')
#    ax2.plot(Freq2, 1.15*Spectrum3MD, '-.', color = Rojo, lw=2, label = '3.0 M')
#    ax2.plot(Freq2, 1.22*Spectrum15MD, '--', color = Purpura, lw=2, label = '1.5 M')
#    ax2.plot(Freq2, Spectrum0MD, 'k', lw=2, label = '0.0 M')
#    ax2.set_xlim(3000,3700)
#    ax2.set_ylim(0.0,2.0)
#    ax2.legend(loc=1)
#
#    ax2.set_ylabel(r'$\mathrm{Absorbance}$', size = 18)
#    ax2.set_xlabel(r'$\mathrm{Frequency}$ $\mathrm{[cm^{-1}]}$', size = 18)
#
#
#    ax22 = fig.add_axes([0.58, 0.58, 0.19, 0.29])    
#    
#    
#    ax22.plot(Freq2, Spectrum4MD, ':', color = Amarillo, lw=2, label = '6.0 M')
#    ax22.plot(Freq1, Spectrum3MD, '-.', color = Rojo, lw=2)
#    ax22.plot(Freq2, Spectrum15MD, '--', color = Purpura, lw=2, label = '1.5 M')
#    ax22.plot(Freq1, Spectrum0MD, 'k', lw=2)
#    ax22.yaxis.set_ticks(np.arange(0,4.1,1))
#    ax22.xaxis.set_ticks([1500,2500,3500])
#    ax22.set_xlim(1000,4000)
#    ax22.set_ylim(0.0,4.0)
#
#    
#    pdf.savefig(bbox_inches='tight')
##    plt.show()
#    plt.close()





with PdfPages('LinearSpectraOD.pdf') as pdf:
    fig, ax1 = plt.subplots(1, 1, figsize=(8,5))
#            fig.subplots_adjust(hspace=0.03)
    
    ax1.add_patch(Rectangle((2411, 0.0),178,3.5,fill=True, color = 'b', alpha=0.2))
    
    ax1.plot(Freq1, Spectrum3MH, '-', color = 'r', lw=3, label = r'$\mathrm{3.0}\/\/\mathrm{M}$')
    ax1.plot(Freq1, Spectrum15MH, '-', color = 'b', lw=3, label = r'$\mathrm{1.5}\/\/\mathrm{M}$')
    ax1.plot(Freq1, Spectrum0MH, 'k', lw=3, label = r'$\mathrm{0.0}\/\/\mathrm{M}$')
    ax1.set_ylim(0.0,3.5)
    ax1.yaxis.set_ticks([0,1,2,3])
    ax1.set_xlim(2250, 2950)
    ax1.legend(loc=1,fontsize=20,frameon=False)
    
#    ax1.text(2750,2,r'$\mathrm{3.0}\/\/\mathrm{M}$', size = 24)        
#    ax1.text(2750,1.18,r'$\mathrm{1.5}\/\/\mathrm{M}$', size = 24)
#    ax1.text(2750,0.38,r'$\mathrm{0.0}\/\/\mathrm{M}$', size = 24)

    ax1.set_ylabel(r'$\mathrm{Absorbance}$', size = 25)
    ax1.set_xlabel(r'$\mathrm{Frequency}$ $\mathrm{[cm^{-1}]}$', size = 25)
    ax1.tick_params(labelsize=16)
   
   
#    left, bottom, width, height = 
    ax11 = fig.add_axes([0.17, 0.57, 0.4, 0.29])    
    
    ax11.plot(Freq1, Spectrum3MH, '-', color = 'r', lw=2.5)
    ax11.plot(Freq1, Spectrum15MH, '-', color = 'b', lw=2.5)
    ax11.plot(Freq1, Spectrum0MH, 'k', lw=2.5)
    ax11.yaxis.set_ticks(np.arange(0,4.1,1))
    ax11.xaxis.set_ticks([1500,2500,3500])
    ax11.set_xlim(1000,4000)
    ax11.set_ylim(0.0,4.0)
    ax11.tick_params(labelsize=14)
    
    pdf.savefig(bbox_inches='tight')
    plt.show()
    plt.close()



#%%

########################################
##########     ISOTROPIC     ###########
########################################


cd = '/Users/Cota/Documents/PhD/UltraFast/ForsterTransferOH/NaOH/Matrices/'

#####20190306_Water.dat

Freq_Water = np.loadtxt(cd + 'Frequency_H2O.dat',comments='#', unpack=True)

Time_Water = np.loadtxt(cd + 'Time_H2O.dat',comments='#', unpack=True)
IsoFit_Water = np.loadtxt(cd + 'IsoFit_H2O.dat',comments='#', unpack=True)
#IsoExp_Water = np.loadtxt(cd + 'IsoExp_Water.dat',comments='#', unpack=True)
#IsoExpStDev_Water = np.loadtxt(cd + 'IsoExpStDev_Water.dat',comments='#', unpack=True)


RawFreq_Water = np.loadtxt(cd + 'RawFrequency_H2O.dat',comments='#', unpack=True)
RawIsoExp_Water = np.loadtxt(cd + 'RawIsoExp_H2O.dat',comments='#', unpack=True)
RawIsoExpStDev_Water = np.loadtxt(cd + 'RawIsoExpStDev_H2O.dat',comments='#', unpack=True)








Freq_NaOH2 = np.loadtxt(cd + 'Frequency_NaOH_2000mM.dat',comments='#', unpack=True)

Time_NaOH2 = np.loadtxt(cd + 'Time_NaOH_2000mM.dat',comments='#', unpack=True)
IsoFit_NaOH2 = np.loadtxt(cd + 'IsoFit_NaOH_2000mM.dat',comments='#', unpack=True)
#IsoExp_NaOH2 = np.loadtxt(cd + 'IsoExp_NaPh_2M.dat',comments='#', unpack=True)
#IsoExpStDev_NaOH2 = np.loadtxt(cd + 'IsoExpStDev_NaPh_2M.dat',comments='#', unpack=True)

RawTime_NaOH2 = np.loadtxt(cd + 'RawTime_NaOH_2000mM.dat',comments='#', unpack=True)
RawFreq_NaOH2 = np.loadtxt(cd + 'RawFrequency_NaOH_2000mM.dat',comments='#', unpack=True)
RawIsoExp_NaOH2 = np.loadtxt(cd + 'RawIsoExp_NaOH_2000mM.dat',comments='#', unpack=True)
RawIsoExpStDev_NaOH2 = np.loadtxt(cd + 'RawIsoExpStDev_NaOH_2000mM.dat',comments='#', unpack=True)


Population_2M = np.loadtxt(cd + 'Population_NaOH_2000mM.dat',comments='#', unpack=True)
Population_2M_stdev = np.loadtxt(cd + 'Population_StDev_NaOH_2000mM.dat',comments='#', unpack=True)

Signatures_2M = np.loadtxt(cd + 'Signatures_NaOH_2000mM.dat',comments='#', unpack=True)







Freq_NaOH5 = np.loadtxt(cd + 'Frequency_NaOH_5000mM_2.dat',comments='#', unpack=True)

Time_NaOH5 = np.loadtxt(cd + 'Time_NaOH_5000mM_2.dat',comments='#', unpack=True)
IsoFit_NaOH5 = np.loadtxt(cd + 'IsoFit_NaOH_5000mM.dat',comments='#', unpack=True)
#IsoExp_NaOH2 = np.loadtxt(cd + 'IsoExp_NaPh_2M.dat',comments='#', unpack=True)
#IsoExpStDev_NaOH2 = np.loadtxt(cd + 'IsoExpStDev_NaPh_2M.dat',comments='#', unpack=True)

RawTime_NaOH5 = np.loadtxt(cd + 'RawTime_NaOH_5000mM_2.dat',comments='#', unpack=True)
RawFreq_NaOH5 = np.loadtxt(cd + 'RawFrequency_NaOH_5000mM_2.dat',comments='#', unpack=True)
RawIsoExp_NaOH5 = np.loadtxt(cd + 'RawIsoExp_NaOH_5000mM_2.dat',comments='#', unpack=True)
RawIsoExpStDev_NaOH5 = np.loadtxt(cd + 'RawIsoExpStDev_NaOH_5000mM_2.dat',comments='#', unpack=True)






#
#
#
#with PdfPages('IsoTransientSpectra.pdf') as pdf:
#    fig, [ax1,ax2,ax3] = plt.subplots(3, 1, figsize=(7.36,12), sharex='col', sharey='row')
#    fig.subplots_adjust(hspace=0.03, wspace=0.03)
#
##    n = [0,2,4,8,12,20,30]
##    n = [100,25,16,7,4,2,1,0.4]
##    n=[0.4, 0.7, 1, 1.5, 2, 3.2, 7, 15, 100]
#    n=[0.36,0.6, 1, 2, 3.1, 4, 10, 50]
#    colors = plt.cm.gnuplot(np.linspace(0,1,len(n)+1))
#    colors = colors[::-1]     
#    
#    
#    
#    for i in range(len(n)):
#        idx = (np.abs(Time_Water - n[i])).argmin()
#        if i < 3:
#            ax1.plot(Freq_Water, 1000*IsoFit_Water[:,idx], color = colors[i+1], lw=2, label = r'$\mathrm{%s\/\/ps}$' % round(Time_Water[idx],1))
#            ax1.errorbar(RawFreq_Water, 1000*RawIsoExp_Water[:,idx], yerr = 1000*RawIsoExpStDev_Water[:,idx], color = colors[i+1], fmt='.',lw=2)
#
#        else:
#            ax1.plot(Freq_Water, 1000*IsoFit_Water[:,idx], color = colors[i+1], lw=2, label = r'$\mathrm{%s\/\/ps}$' % int(Time_Water[idx]))
#            ax1.errorbar(RawFreq_Water, 1000*RawIsoExp_Water[:,idx], yerr = 1000*RawIsoExpStDev_Water[:,idx], color = colors[i+1], fmt='.',lw=2)
#
#
#    ax1.arrow(2560, -45, 0, 48, head_width=2, head_length=3.0,fc='k', ec='k',lw=2)
#    ax1.text(2564,-45,r'$\mathrm{0.4}\/\/\mathrm{ps}$', size = 20)        
#    ax1.text(2534,3,r'$\mathrm{50}\/\/\mathrm{ps}$', size = 20)
#        
#    ax1.text(2411,-59,r'$\mathrm{[NaOH] = } \mathrm{0\/M}$', size = 24)        
#        
#    handles, labels = ax1.get_legend_handles_labels()  
#    ax1.plot(RawFreq_Water, 0.0*RawFreq_Water, 'k--', lw=1)
##    ax1.set_xlim(2400,2590)
#    ax1.yaxis.set_ticks([0,-20,-40,-60])
#    ax1.set_ylim(-65,12)
##    ax1.legend(handles[::-1], labels[::-1],loc=4,fontsize=14,frameon=False,labelspacing=0.1,borderpad=0)
#    ax1.set_ylabel(r'$\Delta \alpha_\mathrm{iso} \times 10^3$', size = 25)
##    ax1.set_xlabel(r'$\mathrm{Frequency}$ $\mathrm{[cm^{-1}]}$', size = 25)
#    ax1.tick_params(labelsize=16)
#    
#    
#    
#    
#    
#    
#    
#    for i in range(len(n)):
#        idx = (np.abs(Time_NaOH2 - n[i])).argmin()
#        Rawidx = (np.abs(RawTime_NaOH2 - n[i])).argmin()
#        if i < 3:
#            ax2.plot(Freq_NaOH2, 1000*IsoFit_NaOH2[:,idx], color = colors[i+1], lw=2, label = r'$\mathrm{%s\/\/ps}$' % round(Time_NaOH2[idx],1))
#            ax2.errorbar(RawFreq_NaOH2, 1000*RawIsoExp_NaOH2[:,Rawidx], yerr = 1000*RawIsoExpStDev_NaOH2[:,Rawidx], color = colors[i+1], fmt='.',lw=2)
#
#        else:
#            ax2.plot(Freq_NaOH2, 1000*IsoFit_NaOH2[:,idx], color = colors[i+1], lw=2, label = r'$\mathrm{%s\/\/ps}$' % int(Time_NaOH2[idx]))
#            ax2.errorbar(RawFreq_NaOH2, 1000*RawIsoExp_NaOH2[:,Rawidx], yerr = 1000*RawIsoExpStDev_NaOH2[:,Rawidx], color = colors[i+1], fmt='.',lw=2)
#
#    ax2.arrow(2560, -29, 0, 30, head_width=2, head_length=2.0,fc='k', ec='k',lw=2)
#    ax2.text(2564,-29,r'$\mathrm{0.4}\/\/\mathrm{ps}$', size = 20)        
#    ax2.text(2534,1,r'$\mathrm{50}\/\/\mathrm{ps}$', size = 20)
#
#    ax2.text(2411,-32,r'$\mathrm{[NaOH] = } \mathrm{2\/M}$', size = 24)        
#        
#    handles, labels = ax2.get_legend_handles_labels()  
#    ax2.plot(RawFreq_NaOH2, 0.0*RawFreq_NaOH2, 'k--', lw=1)
##    ax2.set_xlim(2400,2595)
#    ax2.yaxis.set_ticks([0,-10,-20,-30,-40])
#    ax2.set_ylim(-35,6)
##    ax2.legend(handles[::-1], labels[::-1],loc=4,fontsize=14,frameon=False,labelspacing=0.1,borderpad=0)
#    ax2.set_ylabel(r'$\Delta \alpha_\mathrm{iso} \times 10^3$', size = 30)
##    ax2.set_xlabel(r'$\mathrm{Frequency}$ $\mathrm{[cm^{-1}]}$', size = 25)
#    ax2.tick_params(labelsize=20)
#
#
#
#
#
#
#
#    for i in range(len(n)):
#        idx = (np.abs(Time_NaOH5 - n[i])).argmin()
#        Rawidx = (np.abs(RawTime_NaOH5 - n[i])).argmin()
#        if i < 3:
#            ax3.plot(Freq_NaOH5, 1000*IsoFit_NaOH5[:,idx], color = colors[i+1], lw=2, label = r'$\mathrm{%s\/\/ps}$' % round(Time_NaOH5[idx],1))
#            ax3.errorbar(RawFreq_NaOH5, 1000*RawIsoExp_NaOH5[:,Rawidx], yerr = 1000*RawIsoExpStDev_NaOH5[:,Rawidx], color = colors[i+1], fmt='.',lw=2)
#
#        else:
#            ax3.plot(Freq_NaOH5, 1000*IsoFit_NaOH5[:,idx], color = colors[i+1], lw=2, label = r'$\mathrm{%s\/\/ps}$' % int(Time_NaOH5[idx]))
#            ax3.errorbar(RawFreq_NaOH5, 1000*RawIsoExp_NaOH5[:,Rawidx], yerr = 1000*RawIsoExpStDev_NaOH5[:,Rawidx], color = colors[i+1], fmt='.',lw=2)
#
#    ax3.arrow(2560, -18, 0, 13, head_width=2, head_length=1.0,fc='k', ec='k',lw=2)
#    ax3.text(2564,-18,r'$\mathrm{0.4}\/\/\mathrm{ps}$', size = 20)        
#    ax3.text(2534,-5,r'$\mathrm{50}\/\/\mathrm{ps}$', size = 20)
#
#    ax3.text(2411,-17.3,r'$\mathrm{[NaOH] = } \mathrm{5\/M}$', size = 24)     
##    ax3.text(2415,-10,r'$\times 3$', size = 30)        
#        
#    handles, labels = ax3.get_legend_handles_labels()  
#    ax3.plot(RawFreq_NaOH2, 0.0*RawFreq_NaOH5, 'k--', lw=1)
#    ax3.xaxis.set_ticks([2420,2460,2500,2540, 2580])
#    ax3.set_xlim(2405,2595)
#    ax3.yaxis.set_ticks([0,-5,-10,-15,-20])
#    ax3.set_ylim(-19,2)
##    ax3.legend(handles[::-1], labels[::-1],loc=4,fontsize=14,frameon=False,labelspacing=0.1,borderpad=0)
#    ax3.set_ylabel(r'$\Delta \alpha_\mathrm{iso} \times 10^3$', size = 30)
#    ax3.set_xlabel(r'$\mathrm{Frequency}$ $\mathrm{[cm^{-1}]}$', size = 30)
#    ax3.tick_params(labelsize=20)
#
#
#
#
#
#    pdf.savefig(bbox_inches='tight')
#    pdf.savefig()
#    plt.show()
#    plt.close()
#
#
#
#
#
#
#
#
#
#
#with PdfPages('IsoTransientSpectra_2.pdf') as pdf:
#    fig, [ax1,ax3] = plt.subplots(2, 1, figsize=(7.36,8), sharex='col', sharey='row')
#    fig.subplots_adjust(hspace=0.03, wspace=0.03)
#
##    n = [0,2,4,8,12,20,30]
##    n = [100,25,16,7,4,2,1,0.4]
##    n=[0.4, 0.7, 1, 1.5, 2, 3.2, 7, 15, 100]
#    n=[0.36, 0.6, 1, 2, 4, 50]
#    colors = plt.cm.gnuplot(np.linspace(0,1,len(n)+1))
#    colors = colors[::-1]     
#    
#    markers = ['o','s','v','d','^','p']
#    LegOrder = [0,2,1,3,4,5]
#    
#    for i in range(len(n)):
#        idx = (np.abs(Time_Water - n[i])).argmin()
#        
#        if Time_Water[idx] < 1.0:
#            ax1.plot(Freq_Water, 1000*IsoFit_Water[:,idx], color = colors[i+1], lw=3, zorder=1)
##            ax1.errorbar(RawFreq_Water, 1000*RawIsoExp_Water[:,idx], yerr = 1000*RawIsoExpStDev_Water[:,idx], color = colors[i+1], fmt='.',lw=2)
#
#            ax1.scatter(RawFreq_Water, 1000*RawIsoExp_Water[:,idx], s=60, marker=markers[i], lw = 1.5, facecolors='w', edgecolors = colors[i+1], zorder=2, label = r'$\mathrm{%s\/\/ps}$' % round(Time_Water[idx],1))
#
#
#        else:
#            ax1.plot(Freq_Water, 1000*IsoFit_Water[:,idx], color = colors[i+1], lw=3, zorder=1)
##            ax1.errorbar(RawFreq_Water, 1000*RawIsoExp_Water[:,idx], yerr = 1000*RawIsoExpStDev_Water[:,idx], color = colors[i+1], fmt='.',lw=2)
#
#            ax1.scatter(RawFreq_Water, 1000*RawIsoExp_Water[:,idx], s=50, marker=markers[i], lw = 1.5, facecolors='w', edgecolors = colors[i+1], zorder=2, label = r'$\mathrm{%s\/\/ps}$' % int(Time_Water[idx]))
#
#
#
#
#
#
#
##    ax1.arrow(2560, -45, 0, 48, head_width=2, head_length=3.0,fc='k', ec='k',lw=2)
##    ax1.text(2564,-45,r'$\mathrm{0.4}\/\/\mathrm{ps}$', size = 20)        
##    ax1.text(2534,3,r'$\mathrm{50}\/\/\mathrm{ps}$', size = 20)
#        
#    ax1.text(2411,-59,r'$\mathrm{[NaOH] = } \mathrm{0\/M}$', size = 24)        
#        
#    handles, labels = ax1.get_legend_handles_labels()  
#    ax1.plot(RawFreq_Water, 0.0*RawFreq_Water, 'k--', lw=1)
##    ax1.set_xlim(2400,2590)
#    ax1.yaxis.set_ticks([0,-20,-40,-60])
#    ax1.set_ylim(-65,12)
##    ax1.legend(handles[::-1], labels[::-1],loc=4,fontsize=14,frameon=False,labelspacing=0.1,borderpad=0)
#    
#    ax1.legend(loc='upper center',ncol=3,bbox_to_anchor=(0.5, 1.31),fontsize=21.45,labelspacing=0.1,borderpad=0)
#    
#    ax1.set_ylabel(r'$\Delta \alpha_\mathrm{iso} \times 10^3$', size = 30)
##    ax1.set_xlabel(r'$\mathrm{Frequency}$ $\mathrm{[cm^{-1}]}$', size = 25)
#    ax1.tick_params(labelsize=20)
#    
#    
#    
#    
#    
#    
#    
##    for i in range(len(n)):
##        idx = (np.abs(Time_NaOH2 - n[i])).argmin()
##        Rawidx = (np.abs(RawTime_NaOH2 - n[i])).argmin()
##        if i < 3:
##            ax2.plot(Freq_NaOH2, 1000*IsoFit_NaOH2[:,idx], color = colors[i+1], lw=2, label = r'$\mathrm{%s\/\/ps}$' % round(Time_NaOH2[idx],1))
##            ax2.errorbar(RawFreq_NaOH2, 1000*RawIsoExp_NaOH2[:,Rawidx], yerr = 1000*RawIsoExpStDev_NaOH2[:,Rawidx], color = colors[i+1], fmt='.',lw=2)
##
##        else:
##            ax2.plot(Freq_NaOH2, 1000*IsoFit_NaOH2[:,idx], color = colors[i+1], lw=2, label = r'$\mathrm{%s\/\/ps}$' % int(Time_NaOH2[idx]))
##            ax2.errorbar(RawFreq_NaOH2, 1000*RawIsoExp_NaOH2[:,Rawidx], yerr = 1000*RawIsoExpStDev_NaOH2[:,Rawidx], color = colors[i+1], fmt='.',lw=2)
##
##    ax2.arrow(2560, -29, 0, 30, head_width=2, head_length=2.0,fc='k', ec='k',lw=2)
##    ax2.text(2564,-29,r'$\mathrm{0.4}\/\/\mathrm{ps}$', size = 20)        
##    ax2.text(2534,1,r'$\mathrm{50}\/\/\mathrm{ps}$', size = 20)
##
##    ax2.text(2411,-32,r'$\mathrm{[NaOH] = } \mathrm{2\/M}$', size = 24)        
##        
##    handles, labels = ax2.get_legend_handles_labels()  
##    ax2.plot(RawFreq_NaOH2, 0.0*RawFreq_NaOH2, 'k--', lw=1)
###    ax2.set_xlim(2400,2595)
##    ax2.yaxis.set_ticks([0,-10,-20,-30,-40])
##    ax2.set_ylim(-35,6)
###    ax2.legend(handles[::-1], labels[::-1],loc=4,fontsize=14,frameon=False,labelspacing=0.1,borderpad=0)
##    ax2.set_ylabel(r'$\Delta \alpha_\mathrm{iso} \times 10^3$', size = 30)
###    ax2.set_xlabel(r'$\mathrm{Frequency}$ $\mathrm{[cm^{-1}]}$', size = 25)
##    ax2.tick_params(labelsize=20)
#
#
#
#
#
#
#
#    for i in range(len(n)):
#        idx = (np.abs(Time_NaOH5 - n[i])).argmin()
#        Rawidx = (np.abs(RawTime_NaOH5 - n[i])).argmin()
#        if i < 2:
#            ax3.plot(Freq_NaOH5, 1000*IsoFit_NaOH5[:,idx], color = colors[i+1], lw=3, label = r'$\mathrm{%s\/\/ps}$' % round(Time_NaOH5[idx],1), zorder=1)
##            ax3.errorbar(RawFreq_NaOH5, 1000*RawIsoExp_NaOH5[:,Rawidx], yerr = 1000*RawIsoExpStDev_NaOH5[:,Rawidx], color = colors[i+1], fmt='.',lw=2)
#
#            ax3.scatter(RawFreq_NaOH5, 1000*RawIsoExp_NaOH5[:,Rawidx], s=60, marker=markers[i], lw = 1.5, facecolors='w', edgecolors = colors[i+1], zorder=2)
#
#
#        else:
#            ax3.plot(Freq_NaOH5, 1000*IsoFit_NaOH5[:,idx], color = colors[i+1], lw=3, label = r'$\mathrm{%s\/\/ps}$' % int(Time_NaOH5[idx]), zorder=1)
##            ax3.errorbar(RawFreq_NaOH5, 1000*RawIsoExp_NaOH5[:,Rawidx], yerr = 1000*RawIsoExpStDev_NaOH5[:,Rawidx], color = colors[i+1], fmt='.',lw=2)
#
#            ax3.scatter(RawFreq_NaOH5, 1000*RawIsoExp_NaOH5[:,Rawidx], s=60, marker=markers[i], lw = 1.5, facecolors='w', edgecolors = colors[i+1], zorder=2)
#
#
#
#
##    ax3.arrow(2560, -18, 0, 13, head_width=2, head_length=1.0,fc='k', ec='k',lw=2)
##    ax3.text(2564,-18,r'$\mathrm{0.4}\/\/\mathrm{ps}$', size = 20)        
##    ax3.text(2534,-5,r'$\mathrm{50}\/\/\mathrm{ps}$', size = 20)
#
#    ax3.text(2411,-17.3,r'$\mathrm{[NaOH] = } \mathrm{5\/M}$', size = 24)     
##    ax3.text(2415,-10,r'$\times 3$', size = 30)        
#        
#    handles, labels = ax3.get_legend_handles_labels()  
#    ax3.plot(RawFreq_NaOH2, 0.0*RawFreq_NaOH5, 'k--', lw=1)
#    ax3.xaxis.set_ticks([2420,2460,2500,2540, 2580])
#    ax3.set_xlim(2406,2595)
#    ax3.yaxis.set_ticks([0,-5,-10,-15,-20])
#    ax3.set_ylim(-19,2)
##    ax3.legend(handles[::-1], labels[::-1],loc=4,fontsize=14,frameon=False,labelspacing=0.1,borderpad=0)
#    ax3.set_ylabel(r'$\Delta \alpha_\mathrm{iso} \times 10^3$', size = 30)
#    ax3.set_xlabel(r'$\mathrm{Frequency}$ $\mathrm{[cm^{-1}]}$', size = 30)
#    ax3.tick_params(labelsize=20)
#
#
#
#
#
#    pdf.savefig(bbox_inches='tight')
##    pdf.savefig()
##    plt.show()
#    plt.close()
#
#
#
#
#




with PdfPages('IsoTransientSpectra_3.pdf') as pdf:
    fig, [ax1,ax3] = plt.subplots(2, 1, figsize=(7.50,8), sharex='col', sharey='row')
    fig.subplots_adjust(hspace=0.03, wspace=0.03)

    ax1.plot(RawFreq_Water, 0.0*RawFreq_Water, 'k--', lw=1, zorder=1)

    n=[0.36, 2, 0.6, 4, 1, 50]
    
    colors = plt.cm.gnuplot(np.linspace(0,1,len(n)+1))
    colors = colors[::-1]     
    
    markers = ['o','s','v','d','^','p']
    
    LegOrder = [0,3,1,4,2,5]
        
    for i in range(len(n)):
        idx = (np.abs(Time_Water - n[i])).argmin()
        
        if Time_Water[idx] < 1.0:
            ax1.plot(Freq_Water[:-2], 1000*IsoFit_Water[:-2,idx], color = colors[LegOrder[i]+1], lw=3, zorder=2)
#            ax1.errorbar(RawFreq_Water, 1000*RawIsoExp_Water[:,idx], yerr = 1000*RawIsoExpStDev_Water[:,idx], color = colors[i+1], fmt='.',lw=2)

            ax1.scatter(RawFreq_Water, 1000*RawIsoExp_Water[:,idx], s=60, marker=markers[LegOrder[i]], lw = 1.5, facecolors='w', edgecolors = colors[LegOrder[i]+1], zorder=3, label = r'$\mathrm{%s\/\/ps}$' % round(Time_Water[idx],1))


        else:
            ax1.plot(Freq_Water[:-2], 1000*IsoFit_Water[:-2,idx], color = colors[LegOrder[i]+1], lw=3, zorder=2)
#            ax1.errorbar(RawFreq_Water, 1000*RawIsoExp_Water[:,idx], yerr = 1000*RawIsoExpStDev_Water[:,idx], color = colors[i+1], fmt='.',lw=2)

            ax1.scatter(RawFreq_Water, 1000*RawIsoExp_Water[:,idx], s=60, marker=markers[LegOrder[i]], lw = 1.5, facecolors='w', edgecolors = colors[LegOrder[i]+1], zorder=3, label = r'$\mathrm{%s\/\/ps}$' % int(Time_Water[idx]))

    
    ax1.text(2411,-59,r'$\mathrm{[NaOH] = } \mathrm{0\/M}$', size = 22)        
        
    handles, labels = ax1.get_legend_handles_labels()  
#    ax1.set_xlim(2400,2590)
    ax1.yaxis.set_ticks([0,-20,-40,-60])
    ax1.set_ylim(-65,12)
#    ax1.legend(handles[::-1], labels[::-1],loc=4,fontsize=14,frameon=False,labelspacing=0.1,borderpad=0)
    
    ax1.legend(loc='upper center',ncol=3,bbox_to_anchor=(0.5, 1.295),fontsize=20,labelspacing=0.1,borderpad=0.19,scatterpoints = 2)
    
    ax1.set_ylabel(r'$\Delta \alpha_\mathrm{iso} \times 10^3$', size = 25)
#    ax1.set_xlabel(r'$\mathrm{Frequency}$ $\mathrm{[cm^{-1}]}$', size = 25)
    ax1.tick_params(labelsize=16)
    
    




    ax3.plot(RawFreq_NaOH2, 0.0*RawFreq_NaOH5, 'k--', lw=1, zorder=1)

    for i in range(len(n)):
        idx = (np.abs(Time_NaOH5 - n[i])).argmin()
        Rawidx = (np.abs(RawTime_NaOH5 - n[i])).argmin()
        if i < 2:
            ax3.plot(Freq_NaOH5, 1000*IsoFit_NaOH5[:,idx], color = colors[LegOrder[i]+1], lw=3, label = r'$\mathrm{%s\/\/ps}$' % round(Time_NaOH5[idx],1), zorder=2)
#            ax3.errorbar(RawFreq_NaOH5, 1000*RawIsoExp_NaOH5[:,Rawidx], yerr = 1000*RawIsoExpStDev_NaOH5[:,Rawidx], color = colors[i+1], fmt='.',lw=2)

            ax3.scatter(RawFreq_NaOH5, 1000*RawIsoExp_NaOH5[:,Rawidx], s=60, marker=markers[LegOrder[i]], lw = 1.5, facecolors='w', edgecolors = colors[LegOrder[i]+1], zorder=3)


        else:
            ax3.plot(Freq_NaOH5, 1000*IsoFit_NaOH5[:,idx], color = colors[LegOrder[i]+1], lw=3, label = r'$\mathrm{%s\/\/ps}$' % int(Time_NaOH5[idx]), zorder=2)
#            ax3.errorbar(RawFreq_NaOH5, 1000*RawIsoExp_NaOH5[:,Rawidx], yerr = 1000*RawIsoExpStDev_NaOH5[:,Rawidx], color = colors[i+1], fmt='.',lw=2)

            ax3.scatter(RawFreq_NaOH5, 1000*RawIsoExp_NaOH5[:,Rawidx], s=60, marker=markers[LegOrder[i]], lw = 1.5, facecolors='w', edgecolors = colors[LegOrder[i]+1], zorder=3)



    ax3.text(2411,-17.3,r'$\mathrm{[NaOH] = } \mathrm{5\/M}$', size = 22)     
#    ax3.text(2415,-10,r'$\times 3$', size = 30)        
        
    handles, labels = ax3.get_legend_handles_labels()  
    ax3.xaxis.set_ticks([2420,2460,2500,2540, 2580])
    ax3.set_xlim(2406,2595)
    ax3.yaxis.set_ticks([0,-5,-10,-15,-20])
    ax3.set_ylim(-19,2)
#    ax3.legend(handles[::-1], labels[::-1],loc=4,fontsize=14,frameon=False,labelspacing=0.1,borderpad=0)
    ax3.set_ylabel(r'$\Delta \alpha_\mathrm{iso} \times 10^3$', size = 25)
    ax3.set_xlabel(r'$\mathrm{Frequency}$ $\mathrm{[cm^{-1}]}$', size = 25)
    ax3.tick_params(labelsize=16)





    pdf.savefig(bbox_inches='tight')
#    pdf.savefig()
    plt.show()
    plt.close()







with PdfPages('IsoTransientSpectra_4.pdf') as pdf:
    fig, [ax1,ax2,ax3] = plt.subplots(3, 1, figsize=(7.50,10), sharex='col', sharey='row')
    fig.subplots_adjust(hspace=0.03, wspace=0.03)

    ax1.plot(RawFreq_Water, 0.0*RawFreq_Water, 'k--', lw=1, zorder=1)

    n=[0.36, 2, 0.6, 4, 1, 50]
    
    colors = plt.cm.gnuplot(np.linspace(0,1,len(n)+1))
    colors = colors[::-1]     
    
    markers = ['o','s','v','d','^','p']
    
    LegOrder = [0,3,1,4,2,5]
        
    for i in range(len(n)):
        idx = (np.abs(Time_Water - n[i])).argmin()
        
        if Time_Water[idx] < 1.0:
            ax1.plot(Freq_Water[:-2], 1000*IsoFit_Water[:-2,idx], color = colors[LegOrder[i]+1], lw=3, zorder=2)
#            ax1.errorbar(RawFreq_Water, 1000*RawIsoExp_Water[:,idx], yerr = 1000*RawIsoExpStDev_Water[:,idx], color = colors[i+1], fmt='.',lw=2)

            ax1.scatter(RawFreq_Water, 1000*RawIsoExp_Water[:,idx], s=60, marker=markers[LegOrder[i]], lw = 1.5, facecolors='w', edgecolors = colors[LegOrder[i]+1], zorder=3, label = r'$\mathrm{%s\/\/ps}$' % round(Time_Water[idx],1))


        else:
            ax1.plot(Freq_Water[:-2], 1000*IsoFit_Water[:-2,idx], color = colors[LegOrder[i]+1], lw=3, zorder=2)
#            ax1.errorbar(RawFreq_Water, 1000*RawIsoExp_Water[:,idx], yerr = 1000*RawIsoExpStDev_Water[:,idx], color = colors[i+1], fmt='.',lw=2)

            ax1.scatter(RawFreq_Water, 1000*RawIsoExp_Water[:,idx], s=60, marker=markers[LegOrder[i]], lw = 1.5, facecolors='w', edgecolors = colors[LegOrder[i]+1], zorder=3, label = r'$\mathrm{%s\/\/ps}$' % int(Time_Water[idx]))

    
    ax1.text(2411,-59,r'$\mathrm{[NaOH] = } \mathrm{0\/M}$', size = 22)        
        
    handles, labels = ax1.get_legend_handles_labels()  
#    ax1.set_xlim(2400,2590)
    ax1.yaxis.set_ticks([0,-20,-40,-60])
    ax1.set_ylim(-65,12)
#    ax1.legend(handles[::-1], labels[::-1],loc=4,fontsize=14,frameon=False,labelspacing=0.1,borderpad=0)
    
    ax1.legend(loc='upper center',ncol=3,bbox_to_anchor=(0.5, 1.355),fontsize=20,labelspacing=0.1,borderpad=0.19,scatterpoints = 2)
    
    ax1.set_ylabel(r'$\Delta \alpha_\mathrm{iso} \times 10^3$', size = 25)
#    ax1.set_xlabel(r'$\mathrm{Frequency}$ $\mathrm{[cm^{-1}]}$', size = 25)
    ax1.tick_params(labelsize=16)
    
    



    ax2.plot(RawFreq_NaOH2, 0.0*RawFreq_NaOH2, 'k--', lw=1, zorder=1)

    for i in range(len(n)):
        idx = (np.abs(Time_NaOH2 - n[i])).argmin()
        Rawidx = (np.abs(RawTime_NaOH2 - n[i])).argmin()
        if i < 2:
            ax2.plot(Freq_NaOH2, 1000*IsoFit_NaOH2[:,idx], color = colors[LegOrder[i]+1], lw=3, label = r'$\mathrm{%s\/\/ps}$' % round(Time_NaOH2[idx],1), zorder=2)
#            ax3.errorbar(RawFreq_NaOH5, 1000*RawIsoExp_NaOH5[:,Rawidx], yerr = 1000*RawIsoExpStDev_NaOH5[:,Rawidx], color = colors[i+1], fmt='.',lw=2)

            ax2.scatter(RawFreq_NaOH2, 1000*RawIsoExp_NaOH2[:,Rawidx], s=60, marker=markers[LegOrder[i]], lw = 1.5, facecolors='w', edgecolors = colors[LegOrder[i]+1], zorder=3)


        else:
            ax2.plot(Freq_NaOH2, 1000*IsoFit_NaOH2[:,idx], color = colors[LegOrder[i]+1], lw=3, label = r'$\mathrm{%s\/\/ps}$' % int(Time_NaOH2[idx]), zorder=2)
#            ax3.errorbar(RawFreq_NaOH5, 1000*RawIsoExp_NaOH5[:,Rawidx], yerr = 1000*RawIsoExpStDev_NaOH5[:,Rawidx], color = colors[i+1], fmt='.',lw=2)

            ax2.scatter(RawFreq_NaOH2, 1000*RawIsoExp_NaOH2[:,Rawidx], s=60, marker=markers[LegOrder[i]], lw = 1.5, facecolors='w', edgecolors = colors[LegOrder[i]+1], zorder=3)



    ax2.text(2411,-29.3,r'$\mathrm{[NaOH] = } \mathrm{2\/M}$', size = 22)     
#    ax3.text(2415,-10,r'$\times 3$', size = 30)        
        
    handles, labels = ax2.get_legend_handles_labels()  
    ax2.xaxis.set_ticks([2420,2460,2500,2540, 2580])
    ax2.set_xlim(2406,2595)
    ax2.yaxis.set_ticks([0,-10,-20,-30])
    ax2.set_ylim(-32,2)
#    ax3.legend(handles[::-1], labels[::-1],loc=4,fontsize=14,frameon=False,labelspacing=0.1,borderpad=0)
    ax2.set_ylabel(r'$\Delta \alpha_\mathrm{iso} \times 10^3$', size = 25)
    ax2.set_xlabel(r'$\mathrm{Frequency}$ $\mathrm{[cm^{-1}]}$', size = 25)
    ax2.tick_params(labelsize=16)










    ax3.plot(RawFreq_NaOH2, 0.0*RawFreq_NaOH5, 'k--', lw=1, zorder=1)

    for i in range(len(n)):
        idx = (np.abs(Time_NaOH5 - n[i])).argmin()
        Rawidx = (np.abs(RawTime_NaOH5 - n[i])).argmin()
        if i < 2:
            ax3.plot(Freq_NaOH5, 1000*IsoFit_NaOH5[:,idx], color = colors[LegOrder[i]+1], lw=3, label = r'$\mathrm{%s\/\/ps}$' % round(Time_NaOH5[idx],1), zorder=2)
#            ax3.errorbar(RawFreq_NaOH5, 1000*RawIsoExp_NaOH5[:,Rawidx], yerr = 1000*RawIsoExpStDev_NaOH5[:,Rawidx], color = colors[i+1], fmt='.',lw=2)

            ax3.scatter(RawFreq_NaOH5, 1000*RawIsoExp_NaOH5[:,Rawidx], s=60, marker=markers[LegOrder[i]], lw = 1.5, facecolors='w', edgecolors = colors[LegOrder[i]+1], zorder=3)


        else:
            ax3.plot(Freq_NaOH5, 1000*IsoFit_NaOH5[:,idx], color = colors[LegOrder[i]+1], lw=3, label = r'$\mathrm{%s\/\/ps}$' % int(Time_NaOH5[idx]), zorder=2)
#            ax3.errorbar(RawFreq_NaOH5, 1000*RawIsoExp_NaOH5[:,Rawidx], yerr = 1000*RawIsoExpStDev_NaOH5[:,Rawidx], color = colors[i+1], fmt='.',lw=2)

            ax3.scatter(RawFreq_NaOH5, 1000*RawIsoExp_NaOH5[:,Rawidx], s=60, marker=markers[LegOrder[i]], lw = 1.5, facecolors='w', edgecolors = colors[LegOrder[i]+1], zorder=3)



    ax3.text(2411,-17.3,r'$\mathrm{[NaOH] = } \mathrm{5\/M}$', size = 22)     
#    ax3.text(2415,-10,r'$\times 3$', size = 30)        
        
    handles, labels = ax3.get_legend_handles_labels()  
    ax3.xaxis.set_ticks([2420,2460,2500,2540, 2580])
    ax3.set_xlim(2406,2595)
    ax3.yaxis.set_ticks([0,-5,-10,-15,-20])
    ax3.set_ylim(-19,2)
#    ax3.legend(handles[::-1], labels[::-1],loc=4,fontsize=14,frameon=False,labelspacing=0.1,borderpad=0)
    ax3.set_ylabel(r'$\Delta \alpha_\mathrm{iso} \times 10^3$', size = 25)
    ax3.set_xlabel(r'$\mathrm{Frequency}$ $\mathrm{[cm^{-1}]}$', size = 25)
    ax3.tick_params(labelsize=16)





    pdf.savefig(bbox_inches='tight')
#    pdf.savefig()
    plt.show()
    plt.close()







#%%

################################################
############     DECOMPOSITION     #############
################################################

with PdfPages('Decomposition1.pdf') as pdf:
    fig, ax1 = plt.subplots(1, 1, figsize=(8,4))
#            fig.subplots_adjust(hspace=0.03)
    
#    ax1.add_patch(Rectangle((2411, 0.0),178,3.5,fill=True, color = 'b', alpha=0.2))
    
    ax1.plot(Freq_NaOH2, Signatures_2M[0]*100.0, '-', color = 'blue', lw=3, label = r'$\sigma_\mathrm{OD}$')
    ax1.plot(Freq_NaOH2, Signatures_2M[1]*100.0, '-', color = 'Forestgreen', lw=3, label = r'$\sigma_\mathrm{OH}$')
    ax1.plot(Freq_NaOH2, Signatures_2M[2]*100.0, '-', color = 'red', lw=3, label = r'$\sigma_\mathrm{therm.}$')
    
    linee = np.arange(2400,2600,1)
    ax1.plot(linee, linee*0.0, 'k--', lw=2)
      
    ax1.xaxis.set_ticks([2430,2475,2510,2545, 2580])
    ax1.set_xlim(2416,2595)
    
    ax1.set_ylim(-7.5,1.5)
#    ax1.yaxis.set_ticks([0,1,2,3])

    ax1.legend(loc=3, fontsize=23, frameon=False, labelspacing=0.25, borderpad=0)
    
#    ax1.text(2750,2,r'$\mathrm{3.0}\/\/\mathrm{M}$', size = 24)        
#    ax1.text(2750,1.18,r'$\mathrm{1.5}\/\/\mathrm{M}$', size = 24)
#    ax1.text(2750,0.38,r'$\mathrm{0.0}\/\/\mathrm{M}$', size = 24)

    ax1.set_ylabel(r'$\Delta \alpha_\mathrm{iso} \times 10^{-2}$', size = 25)
    ax1.set_xlabel(r'$\mathrm{Frequency}$ $\mathrm{[cm^{-1}]}$', size = 25)
    ax1.tick_params(labelsize=16)
   
    ax1.text(2537,0.5, r'$\mathrm{[NaOH] = } \mathrm{2\/M}$', fontsize=22)   
   
   
##    left, bottom, width, height = 
#    ax11 = fig.add_axes([0.17, 0.57, 0.4, 0.29])    
#    
#    ax11.plot(Freq1, Spectrum3MH, '-', color = 'tomato', lw=2.5)
#    ax11.plot(Freq1, Spectrum15MH, '-', color = 'Forestgreen', lw=2.5)
#    ax11.plot(Freq1, Spectrum0MH, 'k', lw=2.5)
#    ax11.yaxis.set_ticks(np.arange(0,4.1,1))
#    ax11.xaxis.set_ticks([1500,2500,3500])
#    ax11.set_xlim(1000,4000)
#    ax11.set_ylim(0.0,4.0)
#    ax11.tick_params(labelsize=14)
    
    pdf.savefig(bbox_inches='tight')
    plt.show()
    plt.close()



with PdfPages('Decomposition2.pdf') as pdf:
    fig, ax1 = plt.subplots(1, 1, figsize=(8,4))
#            fig.subplots_adjust(hspace=0.03)
        
    ax1.plot(Population_2M[0], Population_2M[1]*1.15, '-', color = 'blue', lw=3, label = r'$\sigma_\mathrm{OD}$',zorder=4)
#    ax1.errorbar(Population_2M[0], Population_2M[1], yerr=Population_2M_stdev[1], color='blue', fmt='o',lw=2.5, zorder=2)#,label='isotropic signal '+str(t[i]))

    ax1.plot(Population_2M[0], Population_2M[2], '-', color = 'Forestgreen', lw=3, label = r'$\sigma_\mathrm{OH}$',zorder=4)
#    ax1.errorbar(Population_2M[0], Population_2M[2], yerr=Population_2M_stdev[2], color='Forestgreen', fmt='o',lw=2.5, zorder=2)#,label='isotropic signal '+str(t[i]))
    
    ax1.plot(Population_2M[0], Population_2M[3], '-', color = 'red', lw=3, label = r'$\sigma_\mathrm{hgs}$',zorder=4)
#    ax1.errorbar(Population_2M[0], Population_2M[3], yerr=Population_2M_stdev[3], color='red', fmt='o',lw=2.5, zorder=2)#,label='isotropic signal '+str(t[i]))

    
    ax1.text(0.7,0.30, r'$N_\mathrm{OD} (t)$', fontsize = 22)

    ax1.text(0.49,0.06, r'$N_\mathrm{OH} (t)$', fontsize = 22)

    ax1.text(3.2,0.73, r'$N_\mathrm{therm} (t)$', fontsize = 22)

#    ax1.text(4.25,0.8, r'$N_\mathrm{hgs} (t)$', fontsize = 22)

    linee = np.arange(0,7,0.1)
    ax1.plot(linee, linee*0.0, 'k--', lw=2,zorder=1)

    ax1.text(4.7,-0.14, r'$\mathrm{[NaOH] = } \mathrm{2\/M}$', fontsize=22)   

    
#    linee = np.arange(2400,2600,1)
#    ax1.plot(linee, linee*0.0, 'k--', lw=2)
      
#    ax1.xaxis.set_ticks([2430,2475,2510,2545, 2580])
    ax1.set_xlim(0,7)
    
    ax1.set_ylim(-0.2,1.1)
    ax1.yaxis.set_ticks([0,0.2,0.4,0.6,0.8,1.0])

#    ax1.legend(loc=3, fontsize=23, frameon=False, labelspacing=0.25, borderpad=0)
    
#    ax1.text(2750,2,r'$\mathrm{3.0}\/\/\mathrm{M}$', size = 24)        
#    ax1.text(2750,1.18,r'$\mathrm{1.5}\/\/\mathrm{M}$', size = 24)
#    ax1.text(2750,0.38,r'$\mathrm{0.0}\/\/\mathrm{M}$', size = 24)

    ax1.set_ylabel(r'$\mathrm{Population\/\/dynamics}$', size = 25)
    ax1.set_xlabel(r'$\mathrm{Time\/\/delay\/\/[ps]}$ ', size = 25)
    ax1.tick_params(labelsize=16)
   
   
##    left, bottom, width, height = 
#    ax11 = fig.add_axes([0.17, 0.57, 0.4, 0.29])    
#    
#    ax11.plot(Freq1, Spectrum3MH, '-', color = 'tomato', lw=2.5)
#    ax11.plot(Freq1, Spectrum15MH, '-', color = 'Forestgreen', lw=2.5)
#    ax11.plot(Freq1, Spectrum0MH, 'k', lw=2.5)
#    ax11.yaxis.set_ticks(np.arange(0,4.1,1))
#    ax11.xaxis.set_ticks([1500,2500,3500])
#    ax11.set_xlim(1000,4000)
#    ax11.set_ylim(0.0,4.0)
#    ax11.tick_params(labelsize=14)
    
    pdf.savefig(bbox_inches='tight')
    plt.show()
    plt.close()


#%%

################################################
##########     FORSTER POPULATION     ##########
################################################




cd = '/Users/Cota/Documents/PhD/UltraFast/ForsterTransferOH/NaOH/'


scal = [1.82,2.1,2.6,2.95,3.1,5.4,6.24]
Conc = [1.0,1.5,2.0,2.5,3.0,4.0,5.0]




def Survival(t,Molar,kf,a0):
    Gra = np.exp(-(4.0*np.pi*0.0006022*Molar/3.0)*(np.power(a0,3.0)*(np.exp(-(kf*t)/np.power(a0,6.0))-1.0) + np.sqrt(kf)*np.sqrt(np.pi)*np.sqrt(t)*erf(np.sqrt(kf*t)/np.power(a0,3.0))))
#    Gra = (4.0*np.pi*Rho/3.0)*(np.sqrt(kf)*np.sqrt(np.pi)*np.sqrt(t)*erf(np.sqrt(kf*t)/np.power(a0,3.0)))
    return Gra


def residualForster(param,conc,Surv):
    
    T1, kf = param[0], param[1]
    Scalar = param[2:]
    
    Res = []    
    
    for i in range(len(conc)):
        t, popW = np.loadtxt(cd + 'Matrices/Population_NaOH_%smM.dat' % int(1000*conc[i]), unpack=True)[:2]
        
        diff = np.power(popW*Scalar[i] - np.exp(-t/T1)*Survival(t,conc[i],kf,1.5),2.0)
        
        Res.extend(diff)
        
    return Res



Par0 = np.array([1.8,300.0,1.82,2.1,2.6,2.95,3.1,5.4,6.24])    
    
Popt = least_squares(residualForster, Par0, args=(Conc, Survival), bounds=(0.0,np.inf),  method='trf', loss='linear', max_nfev = 20000)

print (Popt['x'])

ValOpt = Popt['x']

T1_opt, rf_opt = ValOpt[0], ValOpt[1]
scal_opt = ValOpt[2:]



PlotConc = [0,2,4,5,6]

with PdfPages('Survival_DelayCurves.pdf') as pdf:

    fig, ax1 = plt.subplots(1, 1, figsize=(7.52,5))

    colors = plt.cm.gnuplot(np.linspace(0,1,7))
    
    markers = ['o','p','s','>','v','d','^']
    
    j=0
    for i in PlotConc:
        t, popW, popLoc, popH = np.loadtxt(cd + 'Matrices/Population_NaOH_%smM.dat' % int(1000*Conc[i]), unpack=True)
        tz, popW_StDev, popLoc_StDev, popH_StDev = np.loadtxt(cd + 'Matrices/Population_StDev_NaOH_%smM.dat' % int(1000*Conc[i]), unpack=True)

        Timet = np.arange(0,10,0.05)

        ax1.plot(Timet, np.exp(-Timet/T1_opt)*Survival(Timet,Conc[i],rf_opt,1.5), lw=2.5, color = colors[j+1], zorder=2)# label = r'$\mathrm{%d\/\/M}$' % Conc[i], zorder=1)

#        ax1.plot(t, scal_opt[i]*popW, '.', lw = 2.0, color = colors[i])
        
#        ax1.errorbar(t, scal_opt[i]*popW, yerr=scal_opt[i]*popW_StDev, color=colors[i], fmt='.',lw=2.5, zorder=2)#,label='isotropic signal '+str(t[i]))
#        ax1.scatter(t, scal_opt[i]*popW, s=20, marker='o', lw = 2.5, facecolors=colors[i], edgecolors = colors[i], zorder=3)
        
        
        plt.errorbar(t, scal_opt[i]*popW, yerr=scal_opt[i]*popW_StDev, lw=2.5, color=colors[j+1],ecolor=colors[j+1], fmt=markers[i],label = r'$\mathrm{%d\/\/M}$' % Conc[i], mec=colors[j+1],ms=7,mfc='w',mew=2,zorder=3)#,label='isotropic signal '+str(t[i]))
        
        j+=1
        
    tt = np.arange(0,6,0.1)   
    ax1.plot(tt,np.exp(-tt/1.73),'k--',lw=1.5,zorder = 1)
    ax1.plot(tt,np.exp(-tt/1.20),'k--',lw=1.5,zorder = 1)

#    ax1.plot(tt,0.0*tt+1.0,'k-.', lw=0.5, zorder = 1)

#    ax1.arrow(3,0.27, 0.0, -0.21, head_width=0.05, head_length=0.01,fc='k', ec='k',lw=2)

    ax1.add_patch(Rectangle((0.17, 0.025),1.8,0.075,fill=True, color = 'w', ec="black",zorder=2))
    ax1.text(0.3,0.070,r'$\mathrm{Fit}\/\/\mathrm{parameters}$', size = 20,zorder=3)
    ax1.text(0.4,0.045,r'$\mathrm{T_1} = 1.78\/\/\mathrm{ps}$', size = 20,zorder=3)
    ax1.text(0.42,0.03,r'$\mathrm{R_F} = 2.98\/\/\mathrm{\AA}$', size = 20,zorder=3)        
          
#    ax1.text(2.8,0.033,r'$\mathrm{5}\/\/\mathrm{M}$', size = 20)
#    ax1.text(2.8,0.033,r'$\mathrm{5}\/\/\mathrm{M}$', size = 20)

    ax1.set_xlabel(r'$\mathrm{Time}$ $\mathrm{delay}$ $\mathrm{[ps]}$', size = 25)
    ax1.set_ylabel(r'$\mathrm{Population}$  $\mathrm{\mathbf{OD}} \cdots \mathrm{H_2O}$',size=25)

    
    ax1.set_xlim(0.0,5.0)
    ax1.set_ylim(0.02,1.2)
    
    ax1.set_yscale('log')
#    ax1.xaxis.set_major_formatter(LogFormatter(base=10,labelOnlyBase=False))
    ax1.set_yticks([1,0.1])
#    ax1.yaxis.set_minor_formatter(mticker.ScalarFormatter())    
    
#    ax1.get_yaxis().set_major_formatter(ScalarFormatter())
#    ax1.get_yaxis().set_minor_formatter(NullFormatter())    

#    ax1.yaxis.set_major_formatter(FuncFormatter(lambda y,pos: ('{{:.{:1d}f}}'.format(int(np.maximum(-np.log10(y),0)))).format(y)))
    
    ax1.tick_params(labelsize=16)

    ax1.legend(loc=1,fontsize=20,labelspacing=0.2,borderpad=0.5,frameon=False,numpoints=1,handletextpad=-0.2)
    
    
    pdf.savefig(bbox_inches='tight')
    plt.show()
    plt.close()




#%%
cd = '/Users/Cota/Documents/PhD/UltraFast/ForsterTransferOH/NaOH/'

with PdfPages('Survival_DelayCurves_4TOC.pdf') as pdf:

    fig, ax1 = plt.subplots(1, 1, figsize=(7.2,6.4))

    colors = plt.cm.gnuplot(np.linspace(0,1,7))
    
    markers = ['o','p','s','>','v','d','^']
    
    j=0
    for i in PlotConc:
        t, popW, popLoc, popH = np.loadtxt(cd + 'Matrices/Population_NaOH_%smM.dat' % int(1000*Conc[i]), unpack=True)
        tz, popW_StDev, popLoc_StDev, popH_StDev = np.loadtxt(cd + 'Matrices/Population_StDev_NaOH_%smM.dat' % int(1000*Conc[i]), unpack=True)

        Timet = np.arange(0.4,10,0.05)

        ax1.plot(Timet, np.exp(-Timet/T1_opt)*Survival(Timet,Conc[i],rf_opt,1.5), lw=2.5, color = colors[j+1], zorder=2)# label = r'$\mathrm{%d\/\/M}$' % Conc[i], zorder=1)

#        ax1.plot(t, scal_opt[i]*popW, '.', lw = 2.0, color = colors[i])
        
#        ax1.errorbar(t, scal_opt[i]*popW, yerr=scal_opt[i]*popW_StDev, color=colors[i], fmt='.',lw=2.5, zorder=2)#,label='isotropic signal '+str(t[i]))
#        ax1.scatter(t, scal_opt[i]*popW, s=20, marker='o', lw = 2.5, facecolors=colors[i], edgecolors = colors[i], zorder=3)
        
        
        plt.errorbar(t, scal_opt[i]*popW, yerr=scal_opt[i]*popW_StDev, lw=2.5, color=colors[j+1],ecolor=colors[j+1], fmt=markers[i],label = r'$\mathrm{%d\/\/M}$' % Conc[i], mec=colors[j+1],ms=7,mfc='w',mew=2,zorder=3)#,label='isotropic signal '+str(t[i]))
        
        j+=1
        
    tt = np.arange(0,6,0.1)   
#    ax1.plot(tt,np.exp(-tt/1.73),'k--',lw=1.5,zorder = 1)
#    ax1.plot(tt,np.exp(-tt/1.20),'k--',lw=1.5,zorder = 1)

#    ax1.plot(tt,0.0*tt+1.0,'k-.', lw=0.5, zorder = 1)

#    ax1.arrow(3,0.27, 0.0, -0.21, head_width=0.05, head_length=0.01,fc='k', ec='k',lw=2)

#    ax1.add_patch(Rectangle((0.17, 0.025),1.8,0.075,fill=True, color = 'w', ec="black",zorder=2))
#    ax1.text(0.3,0.070,r'$\mathrm{Fit}\/\/\mathrm{parameters}$', size = 20,zorder=3)
#    ax1.text(0.4,0.045,r'$\mathrm{T_1} = 1.78\/\/\mathrm{ps}$', size = 20,zorder=3)
#    ax1.text(0.42,0.03,r'$\mathrm{R_F} = 2.98\/\/\mathrm{\AA}$', size = 20,zorder=3)        

    ax1.text(0.42,0.185,r'$\mathrm{[OH^-]}$', size = 24,zorder=3)        
    ax1.hlines(0.157,0.3,1.43,colors='k')
          
#    ax1.text(2.8,0.033,r'$\mathrm{5}\/\/\mathrm{M}$', size = 20)
#    ax1.text(2.8,0.033,r'$\mathrm{5}\/\/\mathrm{M}$', size = 20)

    ax1.set_xlabel(r'$\mathrm{Time}$ $\mathrm{[ps]}$', size = 30)
    ax1.set_ylabel(r'$\mathrm{Excited}$ $\mathrm{OD}$  $\mathrm{vibrations} $',size=30)

    
    ax1.set_xlim(0.0,5.0)
    ax1.set_ylim(0.02,1.2)
    
    ax1.set_yscale('log')
#    ax1.xaxis.set_major_formatter(LogFormatter(base=10,labelOnlyBase=False))
    ax1.set_yticks([1,0.1])
#    ax1.yaxis.set_minor_formatter(mticker.ScalarFormatter())    
    
#    ax1.get_yaxis().set_major_formatter(ScalarFormatter())
#    ax1.get_yaxis().set_minor_formatter(NullFormatter())    

#    ax1.yaxis.set_major_formatter(FuncFormatter(lambda y,pos: ('{{:.{:1d}f}}'.format(int(np.maximum(-np.log10(y),0)))).format(y)))
    
    ax1.tick_params(labelsize=20)

    ax1.legend(loc=3,fontsize=24,labelspacing=0.4,borderpad=0.1,frameon=False,numpoints=1,handletextpad=0.2)
    
    
    pdf.savefig(bbox_inches='tight')
    plt.show()
    plt.close()


#%%


#ax2_T=fig.add_axes([0.775, 0.775, .125, .125], axisbg=None, xscale='log', xlim=(1,50),ylim=(0,80))
#ax2_T.xaxis.set_major_formatter(LogFormatter(base=10,labelOnlyBase=False))
#ax2_T.tick_params(axis='both', labelsize=6.5)
#ax2_T.set_xticks([1,5,10,20,40])
#ax2_T.set_yticks([0,20,40,60])





Time, ForsterD = np.loadtxt(cd + 'ForsterDistance.dat', delimiter = ',', unpack=True)

Time1, Time2 = Time[:443], Time[443:]
ForsterD11, ForsterD22 = ForsterD[:443], ForsterD[443:]

Time1 = Time1[35:]
Time2 = Time2[39:]

ForsterD1 = ForsterD11[35:]
ForsterD2 = ForsterD22[39:]


f = interp1d(Time1, ForsterD1)
Distance = np.arange(1.514,49.9,0.005)
DisForster1 = f(Distance)

#for i in range(len(Distance)):
#    Rate = sum(DisForster1[:i+1])/sum(DisForster1)
##    print(Rate)
#    if Rate<0.5:
#        ForD = Distance[i+1]
     
ForD1 = sum(DisForster1*Distance)/sum(DisForster1)    
     
print ('Mean hydroxide = ', ForD1)
    
    
    
f2 = interp1d(Time2, ForsterD2)
Distance = np.arange(1.514,49.9,0.005)
DisForster2 = f2(Distance)

#for i in range(len(Distance)):
#    Rate = sum(DisForster2[:i+1])/sum(DisForster2)
##    print(Rate)
#    if Rate<0.5:
#        ForD2 = Distance[i+1]
    

ForD2 = sum(DisForster2*Distance)/sum(DisForster2)    

print ('Mean proton (Rutger) = ', ForD2)
    

#%%


with PdfPages('ForsterDistribution.pdf') as pdf:

    fig, ax1 = plt.subplots(1, 1, figsize=(7.52,5))

    colors = plt.cm.gnuplot(np.linspace(0,1,8))


    ax1.axvline(x=ForD1, color = 'r', ls= '--')
    ax1.axvline(x=ForD2, color = 'b', ls= '--')
    

    ax1.plot(Time1, ForsterD1/sum(ForsterD1), lw=1.5, color = 'r', label = r'$\mathrm{OH^-}$')
    ax1.plot(Time2, ForsterD2/sum(ForsterD2), lw=1.5, color = 'b', label = r'$\mathrm{H^+}$')
#    ax1.plot(Time1[:47], ForsterD1[:47], lw=1.5, color = 'k')
    
#    ax1.plot(Distance, DisForster, 'k--',lw=1.0)
    
#    ax1.arrow(ForD, 0.0102, 0, -0.00931, head_width=0.2, head_length=0.0007,fc='k', ec='k',lw=2)
    
#    ax1.arrow(ForD2, 0.00302, 0, -0.00215, head_width=0.2, head_length=0.0007,fc='k', ec='k',lw=2)

#    ax1.axvline(x=2.6)
#    ax1.axvline(x=3.2)

    ax1.set_xlabel(r'$\mathrm{Distance}$ $\mathrm{[\AA]}$', size = 24)
    ax1.set_ylabel(r'$\mathrm{F\"orster}$ $\mathrm{transfer}$ $\mathrm{distribution}$',size=24)

    ax1.set_xlim(0.0,10.0)
    
#    ax1.yaxis.set_ticks([-1,1])

#    ax1.set_ylim(0.000,0.016)
    ax1.tick_params(labelsize=16)
    
    ax1.legend(loc=1,fontsize=20)
    
    pdf.savefig(bbox_inches='tight')
#    plt.show()
    plt.close()














cd = ('/Users/Cota/Documents/PhD/UltraFast/ForsterTransferOH/NaOD/')

Con = [1.0, 2.0, 3.0, 4.5, 6.0]

Scalar = [1.5,1.75,2.2,3.5,4.1]

with PdfPages('Survival_DelayCurves3400.pdf') as pdf:

    fig, ax1 = plt.subplots(1, 1, figsize=(7.52,5.075))

    markers = ['o','s','v','d','^']

    colors = plt.cm.gnuplot(np.linspace(0,1,7))
#    colors = colors[::-1]
    
    for i in range (5):
        t, popW, popIon, popH = np.loadtxt(cd + 'Matrices/Population_NaOD_%smM.dat' % int(1000*Con[i]), unpack=True)
        tz, popW_StDev, popIon_StDev, popH_StDev = np.loadtxt(cd + 'Matrices/Population_StDev_NaOD_%smM.dat' % int(1000*Con[i]), unpack=True)

        Timet = np.arange(0,4,0.02)

#        ax1.plot(Timet, np.exp(-Timet/T1_opt)*Survival(Timet,Conc[i],rf_opt,1.5), lw=1.5, color = colors[i], label = r'$\mathrm{%s\/\/M}$' % Conc[i])
#
#        ax1.plot(t, scal_opt[i]*popW, '.', lw = 2.0, color = colors[i])
        
        plt.errorbar(t, Scalar[i]*popW, yerr=Scalar[i]*popW_StDev, lw=2.5, color=colors[i+1],ecolor=colors[i+1], fmt=markers[i],label = r'$\mathrm{%s\/\/M}$' % Con[i], mec=colors[i+1],ms=7,mfc='w',mew=2,zorder=2)#,label='isotropic signal '+str(t[i]))
#        plt.scatter(t, Scalar[i]*popW, s=60, marker=markers[i], lw = 0.1, facecolors=colors[i], edgecolors = colors[i], zorder=3)#,label='isotropic signal '+str(t[i]))

#        ax1.scatter(t, scal_opt[i]*popW, s=20, marker='o', lw = 2.5, facecolors=colors[i], edgecolors = colors[i], zorder=3)

        
    tt=np.arange(0,4,0.02)
    ax1.plot(tt,np.exp(-tt/0.75),'k', lw=2,zorder=1)

    ax1.set_xlabel(r'$\mathrm{Time}$ $\mathrm{delay}$ $\mathrm{[ps]}$', size = 24)
    ax1.set_ylabel(r'$\mathrm{Population}$  $\mathrm{\mathbf{OH}} \cdots \mathrm{D_2O}$',size=24)

    ax1.set_yscale('log')
    ax1.set_xlim(0.0,3.0)
    ax1.set_ylim(0.02,1.2)
    ax1.tick_params(labelsize=16)

    ax1.legend(loc=3,labelspacing=0.2,borderpad=0.5,frameon=False,numpoints=1,handletextpad=-0.2,fontsize=20)
    
    ax11 = fig.add_axes([0.5, 0.617, 0.38, 0.25])    
    
#    ax11.plot(Freq1, Spectrum3MH, '-', color = 'tomato', lw=2.5)
#    ax11.plot(Freq1, Spectrum15MH, '-', color = 'Forestgreen', lw=2.5)
#    ax11.plot(Freq1, Spectrum0MH, 'k', lw=2.5)
    ax11.yaxis.set_ticks(np.arange(0,4.1,1))
    ax11.xaxis.set_ticks([1500,2500,3500])
    ax11.set_xlim(1000,4000)
    ax11.set_ylim(0.0,4.0)
    ax11.tick_params(labelsize=14)
    
    
    ax11.plot(Freq2, 1.25*Spectrum4MD, color = 'tomato', lw=2, label = r'$\mathrm{6\/M}$')
    ax11.plot(Freq2, 1.15*Spectrum3MD, color = 'Forestgreen', lw=2, label = r'$\mathrm{3\/M}$')
#    ax11.plot(Freq2, 1.22*Spectrum15MD, color = 'Forestgreen', lw=2, label = '1.5 M')
    ax11.plot(Freq2, Spectrum0MD, 'k', lw=2, label = r'$\mathrm{0\/M}$')
    
    ax11.legend(loc=1,fontsize=15,labelspacing=0.0,borderpad=0.0,frameon=False)
    
    ax11.set_ylabel(r'$\mathrm{Absorbance}$', size = 17)
    ax11.set_xlabel(r'$\mathrm{Frequency}$ $\mathrm{[cm^{-1}]}$', size = 17)
    
    
    pdf.savefig(bbox_inches='tight')
    plt.show()
    plt.close()





       