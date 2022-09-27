# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 11:52:18 2017

@author: cthompson
"""

from ramanTools.RamanSpectrum import *

def removewater(ftirspectrum):
    
    water = RamanSpectrum("C:/Users/cthompson/Dropbox/Data/WaterVapor.csv")
    
    if len(water[900:1700])!=len(ftirspectrum[900:1700]):
            print ('error')
            return None
    def func(A):return sum(abs(diff(ftirspectrum[1290:2100]-A*water[1290:2100])))
    A = scipy.optimize.minimize(func, 0.01).x[0]
    ftirspectrum[1290:2100]-= A*water[1290:2100]
    
#    def func(A):return sum(abs(diff(ftirspectrum[3340:4000]-A*water[3340:4000])))
#    B = scipy.optimize.minimize(func, A).x[0]
#    ftirspectrum[3340:4000]-= B*water[3340:4000]
    return ftirspectrum
    
def removeCO2(ftirspectrum):
    
    water = RamanSpectrum("C:/Users/cthompson/Dropbox/Data/WaterVapor.csv")
    
    if len(water[2310:2400])!=len(ftirspectrum[2310:2400]):
            print ('error')
            return None
    def func(A):return sum(abs(diff(ftirspectrum[2310:2400]-A*water[2310:2400])))
    A = scipy.optimize.minimize(func, 0.01).x[0]
    ftirspectrum[2310:2400]-= A*water[2310:2400]
    
    

    return ftirspectrum

def fittosum(spectrumtofit, spectratofitwith,rnge=None):
    
    
    xs =array(spectrumtofit.index[spectrumtofit.nearest(600):spectrumtofit.nearest(1800)])
    def func(x, A, B,m,b):return A*spectratofitwith[0][600:1800] + B*spectratofitwith[1][600:1800]+m*xs+b
    A = scipy.optimize.curve_fit(func,xs,spectrumtofit[600:1800], [1,1,0,0])
    return A#xs,func(xs,*A[0])
    
    
    
if __name__ == "__main__":
#    a =RamanSpectrum("T:/01 Active Programs/Anaerobic Membrane Bioreactor/Data/032217/CT146-60A 1 hr.csv")
#    water = RamanSpectrum("C:/Users/cthompson/Dropbox/Data/WaterVapor.csv")
#    water.plot()
#    #a+=0.16*water
#    a.plot(color='g')
#    removewater(a)
##    a+=water*0.5
#    
#    a.plot(color='r')
    pass

    