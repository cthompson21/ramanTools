# -*- coding: utf-8 -*-
"""
Created on Wed May  6 15:50:30 2015

@author: chris
"""

import sys
if 'C:/Users/cthompson/Dropbox/PyScripts/NT' not in sys.path:
    sys.path.append('C:/Users/cthompson/Dropbox/PyScripts/NT')
import pdb
import numpy as np
from numpy import *
import scipy.optimize
import pandas
from ramanTools.SPETools import File
from copy import deepcopy,copy
import scipy.optimize
import matplotlib.pyplot as plt
from collections import namedtuple
import inspect
import os
import numpy.polynomial.polynomial





class IRSpectrum(pandas.Series):

    def __init__(self,filename,avg = True):
        
      
        if type(filename) == pandas.Series:
           
            obj = pandas.Series.__init__(self,filename.values,array(filename.index))
            self.info = str()
            self.num_frames = 0
            self.accum_time = float()
            
            
        elif filename[-4:] == '.SPE' or filename[-4:] == '.spe':
            
                
                fid = File(filename)
                fid._load_size()
                x_offset = fid.read_at(3103, 1, float64)
                a = fid.read_at(3263, 6, float64)# wavelengths in nm
                wl = fid.read_at(3311,1, float64)  ## laser wavelength
                spectrum = fid.load_img()
                
                if avg== True:
                    
                    if fid._numframes>1:
                        
                        spectrum=sum(spectrum,axis=0)
                        spectrum=spectrum.flatten()
                    elif fid._numframes==1:
                        spectrum=spectrum.flatten()
                        

                a = list(a)
                a.reverse()
                a = polyeval(a,np.arange(len(spectrum))+x_offset)
                a = -(1E7/a-1E7/wl)
                
                
                obj = pandas.Series.__init__(self,spectrum,a)
                self.info = fid.get_info()
                self.accum_time=fid.get_accum_time()
                self.num_frames = fid._numframes
                self.name = filename
                fid.close()
        elif filename[-4:] == '.SSM':
            a = loadtxt(filename,delimiter ='  ', unpack = True,skiprows =2)
            obj = pandas.Series.__init__(self,a[1],a[0],dtype = float)#pandas.Series.__init__(self,a[1][0:cutoff],a[0][0:cutoff],dtype = float)
            
            
        elif filename[-4:] == '.csv':
            try: ### for UVVis files
                a= np.loadtxt(filename,
                          delimiter = ',',
                          usecols = (0,1),
                          unpack = True)
                
            except:
                try:
                    a= np.loadtxt(filename,
                              delimiter = '\t',
                              usecols = (0,1),
                                unpack = True)

                  
                except:
                    try: ### For FTIR spectra
                        a= np.genfromtxt(filename, delimiter = ',',skip_header = 19, skip_footer =35,unpack=True)
                       
    
                      
                    except:
                        obj = pandas.Series.__init__(self,array([0]),array([0]))
                        self.info = str()
                        self.accum_time=-1
                        self.num_frames = -1
                        self.name = os.path.basename(filename)
                        print ('error opening Raman Spectrum File', filename, ' error 1')
                        raise
                        return obj
            
            if any(abs(np.diff(a[0]))>100):   
                    cutoff = argmax(abs(diff(a[0])))
                    a[1][0:cutoff] = (a[1][0:cutoff]+a[1][cutoff:cutoff*2])/2
                    a = a[:,0:cutoff]
                    
            else:
                    cutoff = a[1].size
                        
            if a[0,-1]<a[0,0]:
                a=a[:,::-1]
            obj = pandas.Series.__init__(self,a[1],a[0],dtype = float)#pandas.Series.__init__(self,a[1][0:cutoff],a[0][0:cutoff],dtype = float)
            self.info = str()
            self.accum_time=-1
            self.num_frames =-1
            self.name = os.path.basename(filename)
        
        return obj
     
    def __array_finalize__(self,obj):
        if obj is None: return

        
    def reverse(self):
        return RamanSpectrum(pandas.Series(self.values[::-1],self.index[::-1]))
    def addoffset(self,offset):
        
        return RamanSpectrum(pandas.Series(self.values,array(self.index)+offset))
    def nearest(self, x):
        if x<self.index[0]:
            return 0
        elif x> self.index[-1]:
            return len(self.index)-1
        else:
            return np.argmax(np.diff(sign(self.index-x)))+1
    def __getitem__(self,key):
        
        try:
            if isinstance(key, slice):
                if key.stop!=None and key.start!=None:
                    new_key=slice(self.nearest(key.start),self.nearest(key.stop))
                elif key.stop!=None and key.start==None:
                    new_key=slice(None,self.nearest(key.stop))
                elif key.stop==None and key.start!=None:
                    new_key=slice(self.nearest(key.start),None)
                else:
                    new_key=slice(None,None)
                return self.values[new_key]
            elif isinstance(key,float) or isinstance(key,int):
                return self.get_values()[self.nearest(key)]
            elif isinstance(key,list) or isinstance(key,np.ndarray):
                return array([self[i] for i in key])
            else:
                
                return []
                
                
        except:
            raise
        return -1
        
    def findminbetween(self,rnge):
        return findminbetween(self,rnge)


        
    def copy(self):
        return self
    def smooth(self,window_len=3,window='flat'):
        self = smooth(self,window_len=window_len,window=window)
        return None
    def correctjump(self,xvalue):
        self=correctjump(self,xvalue)
        return None
    def autobaseline(self,rnge,order = 0,join=None,partial_fit_region=None,specialoption = None):
        self = autobaseline(self,rnge,order = order,join=join,partial_fit_region=partial_fit_region,specialoption =specialoption)
        return None

        
    def set_name(self,filename):
        self.name = filename
        return 0

    def smoothbaseline(self, rnge1,rnge2,ax=None,_plot= False):
        return smoothbaseline(self, rnge1,rnge2,_plot= _plot,ax=ax)

        
    def calc_noise(self,rnge,rnge_type = 'data'):
        return calc_noise(self,rnge,rnge_type = rnge_type)
        
    def calc_area(self,rnge,fill=False):
        
        return calc_area(self,rnge,fill=fill)
    def normalize(self):
        rnge =(self.index[0], self.index[-1])
        self = normalize(self,rnge)
        return None
        
        
def SGsmooth(x,y, width=11,order = 2):#data,rnge,rnge_type = 'data'):
    
    retval = np.ndarray(y.shape)

    for i in range(y.size):

        i_min = max(0,i-width/2)
        i_max = min( i+width/2-1, y.size)
        fit = np.polyfit(x[i_min:i_max+1],y[i_min:i_max+1],order)
        retval[i] = polyeval(fit,x[i])
    return retval


def polyeval(constants,x):
    n = 0 
    for i in range(len(constants)):
       
        n+= constants[i]*x**(len(constants)-i-1)
    return n
    
def add_RamanSpectra(_x,_y,mode = 'union'):
    if mode == 'union':
        xmax = min(max(array(_x.index)),max(array(_y.index)))
    
        xmin = max(min(array(_x.index)),min(array(_y.index)))
     
        _x2 = _x.truncate(xmin,xmax)
        _y2 = _y.truncate(xmin,xmax)
        _y2 = _y2.reindex(_x2.index,method = 'ffill')
    elif mode == 'join':
        _y2 = _y2.reindex(_x.index,method = 'ffill')
        
    
    output = _x2+_y2
    if any(output.isnull()):
        output = output.fillna(method = 'bfill')
    
    return RamanSpectrum(output)

def subtract_RamanSpectra(_x,_y):
    
    xmax = min(max(array(_x.index)),max(array(_y.index)))
    xmin = max(min(array(_x.index)),min(array(_y.index)))
   
    _x2 = _x.truncate(xmin,xmax)
    _y2 = _y.truncate(xmin,xmax)
    _y2 = _y2.reindex(_x2.index,method = 'ffill')
    
    output = _x2-_y2
    if any(output.isnull()):
        output = output.fillna(method = 'bfill')
    
    return RamanSpectrum(output)
    
    

    
    
def FourierFilter(input_array,width =900,demo = False):
    

    
    r = fft.fft(input_array)
    if demo == True:
        plot(1/input_array.index,r.real)
        plot(r.imag)
        vlines(1/(r.size/2-width),0,100)
        vlines(1/(r.size/2+width),0,100)

    r[r.size/2-width:r.size/2+width]=0
    s = fft.ifft(r).real
    output = pandas.Series(s,array(input_array.index))
    if demo==True:
        pass
    
    

    return RamanSpectrum(output)
    
    
def calc_area(spectrum,rnge,fill = False):
    
    """calculates the area between data and a straight line (a linear baseline)"""
    start = spectrum.nearest(rnge[0])#argmin(abs(array(spectrum.index)-rnge[0]))
    end = spectrum.nearest(rnge[1])#argmin(abs(array(spectrum.index)-rnge[1]))+1
    
    xs = array(spectrum.index[start:end])
    ys= spectrum[rnge[0]:rnge[1]]#.values[start:end]   

    try:
        slope =(ys[-1]-ys[0])/(xs[-1]-xs[0])
        baseline = slope*(xs-xs[0])+ys[0]
        
    except:
        print( 'error in area calculation.')
        return 0
    
    
   
    if fill:
        plt.fill_between(xs,ys,baseline)
        

    return sum((ys-baseline)[1:]*np.diff(xs))
    
def findminbetween(spectrum,rnge):
     
    """returns frequqncy of miniumum intensity within rnge"""
    start = spectrum.nearest(rnge[0])#argmin(abs(array(spectrum.index)-rnge[0]))
    end = spectrum.nearest(rnge[1])#argmin(abs(array(spectrum.index)-rnge[1]))+1
    
    xs = array(spectrum.index[start:end])
    ys= spectrum[rnge[0]:rnge[1]]#.values[start:end]   

    
    return xs[argmin(ys)]
    
def correctjump(spectrum,xvalue):
    x= spectrum.nearest(xvalue)
    b=(spectrum.values[x+1]-spectrum.values[x])#*ones(spectrum.values[x+1:].shape)
    spectrum.values[x+1:]-=b
    return spectrum
    
def autobaseline(spectrum,rnge,order = 0,join = None,partial_fit_region=None,specialoption = None):
    
   
    start = spectrum.nearest(rnge[0])
    end = spectrum.nearest(rnge[-1])
    
 
    if start == end:
        print ('error.  same start and end points')
        return spectrum
    xs = array(spectrum.index[start:end])
    ys= spectrum.values[start:end]
    if specialoption=='points':
        xfits = array(list(rnge))
        yfits = spectrum[xfits]
        r = np.polyfit(xfits,yfits,order)
        spectrum.values[start:end] = spectrum.values[start:end]-polyeval(r,xs)
        spectrum[:]-=min(spectrum)
        if join==None:
            pass
        elif join == 'start':
            if start==0:
                pass
            else:
                offset = spectrum.iloc[start]-spectrum.iloc[start-1]
                spectrum.values[start:] -= offset
        elif join == 'end':
            if end == spectrum.size-1:
                print ('end point too far')
            else:
                offset = spectrum.iloc[end]-spectrum.iloc[end-1]
                print (spectrum.index[end],spectrum.index[end-1])
                spectrum.values[end:] -= offset
        else:
            print( 'wrong join option in autobaseline:', join)
        return spectrum
    elif partial_fit_region!=None:
        lo_start = spectrum.nearest(partial_fit_region[0])
        lo_end = spectrum.nearest(partial_fit_region[1])
#        xfits = np.delete(xs,range(lo_start,lo_end))
#        yfits = np.delete(ys,range(lo_start,lo_end))
        xfits = array(spectrum.index[lo_start:lo_end])
        yfits= spectrum.values[lo_start:lo_end]
        
    else:
        xfits = xs
        yfits = ys

    if order == 0:
        slope =(yfits[-1]-yfits[0])/(xfits[-1]-xfits[0])
        b = yfits[0]-slope*xfits[0]
        baseline = slope*(xs)+b
        spectrum.values[start:end] = spectrum.values[start:end]-baseline
    else:

        r = np.polynomial.polynomial.polyfit(xfits,yfits,order)
        
        spectrum.values[start:end] = spectrum.values[start:end]-np.polynomial.polynomial.polyval(xs,r)
    if join == None:
        pass
    elif join=='start':
        if start==0:
            pass
        else:
            offset = spectrum.iloc[start]-spectrum.iloc[start-1]
            spectrum.values[start:] -= offset
    elif join == 'end':
        if end == spectrum.size-1:
            print ('end point too far')
        else:
            offset = spectrum.iloc[end]-spectrum.iloc[end-1]
            print( spectrum.index[end],spectrum.index[end-1])
            spectrum.values[end:] -= offset
    else :
        print ('wrong join option in autobaseline:', join)
    spectrum[:]-=min(spectrum.values)
    return spectrum
    


def smooth(spectrum,window_len=5,window='flat'):
    
    spectrum =spectrum.copy()
    if spectrum.ndim != 1:
        raise ValueError ("smooth only accepts 1 dimension arrays.")

    if spectrum.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")


    if window_len<3:
        return spectrum
 
    if not window in ['flat','SG']:
        raise ValueError( "Window is on of 'flat', 'SG")


    s=np.r_[spectrum.values[window_len-1:0:-1],spectrum,spectrum.values[-1:-window_len:-1]]
   
    if window == 'SG':
        
        order = 2
        y=spectrum.values
        x=np.array(spectrum.index)
        retval = np.ndarray(y.shape)
        for i in range(y.size):
            
            i_min = max(0,i-window_len/2)
            i_max = min( i+window_len/2-1, y.size)
            fit = np.polynomial.polynomial.polyfit(x[i_min:i_max+1],y[i_min:i_max+1],min(order,i_max-i_min))
            retval[i] = polyeval(fit,x[i])
        

        spectrum[:]=retval[:]
        return spectrum
        
            
    elif window == 'flat': #moving average
        w=np.ones(window_len,'d')
        spectrum.values[:] =np.convolve(w/w.sum(),s,mode='valid')[(window_len-1)/2:-(window_len-1)/2]
    else:
        w=eval(window+'(window_len)')
        spectrum.values[:] =np.convolve(w/w.sum(),s,mode='valid')[(window_len-1)/2:-(window_len-1)/2]
   
    
    return spectrum 
    
def calc_noise(spectrum,rnge,rnge_type = 'data'):
    
    start = spectrum.nearest(rnge[0])
    end = spectrum.nearest(rnge[-1])
    
   
    xs = array(spectrum.index[start:end])
    ys= spectrum.values[start:end]
  
    r = np.polynomial.polynomial.polyfit(xs,ys,2)
    average = np.polynomial.polynomial.polyval(r,xs)
    
    standarddev = np.sqrt(np.mean((ys-average)**2) - np.mean(ys-average)**2)
    
    return standarddev

def removespikes(spectrum,thresholds=[10,5]):
   
    
    spectrum = spectrum.copy()
    start = 0
    end = len(spectrum)-1
    for s in range(start,end,64):
        e = min(s+63,end)
        for threshold in thresholds:
            try:
                noise = calc_noise(spectrum, (spectrum.index[s],spectrum.index[e]))
            except:
               
               continue
            
            change = np.append(0,np.diff(spectrum.values[s:e]))
            spikes_x = np.where(abs(change)>noise*threshold)[0]+s
            
            for i in spikes_x:
                try:
                    spectrum.values[i] = np.mean(spectrum.values[np.array([i-4,i-3,i-2,i+2,i+3,i+4])])
                except:
                    pass
    return spectrum

def normalize(spectrum,rnge):
    start = spectrum.nearest(rnge[0])
    end = spectrum.nearest(rnge[1])
    spectrum = spectrum.copy()
    spectrum.values[:]-=min(spectrum.values[start:end])
    spectrum.values[:]/=max(spectrum.values[start:end])
    return spectrum



    
    
def ics(orca_output_file, normalize = False,color='k',labelpeaks = True):
    fileopen = open(orca_output_file,'rb')
    f=fileopen.readlines()
    fileopen.close()
    for l in f:
        if 'RAMAN SPECTRUM' in l:
            start =  f.index(l)+5
    table=list()
    for z in f[start:]:
        i = z.split(' ')
        while '' in i: 
            i.remove('')
       
        if i[0][0].isdigit():
           
            table.append([float(i[1]),float(i[2])])
            if labelpeaks:   
                plt.gca().annotate(i[0], (float(i[1]), float(i[2])+0.2), color=color,fontsize = 8,horizontalalignment='center') 
        elif not i[0][0].isdigit():
            break
    table = transpose(table)
    if normalize == True:
            table[1]/=max(table[1])
    return plt.vlines(table[0],0,table[1],linewidth = 2,color=color)
    
def irs(orca_output_file, normalize = False,color='k',labelpeaks = True):

    fileopen = open(orca_output_file,'rb')
    f=fileopen.readlines()
    fileopen.close()
    for l in f:
        if 'IR SPECTRUM' in l:
            start =  f.index(l)+5
    table=list()
    for z in f[start:]:
        i = z.split(' ')
        while '' in i: 
            i.remove('')
       
        if i[0][0].isdigit():
           
            table.append([float(i[1])-0.5,float(i[2])])
            if labelpeaks:   
                plt.gca().annotate(i[0], (float(i[1]), float(i[2])+0.2), color=color,fontsize = 8,horizontalalignment='center') 
        elif not i[0][0].isdigit():
            
            break
    table = transpose(table)
    if normalize == True:
            table[1]/=max(table[1])
    return plt.vlines(table[0],0,table[1],linewidth = 2,color=color)
    

    

def etchegoin_analysis(folder):
    
    os.chdir(folder)
    data = zeros((0,1024))
    frequency = zeros(data.shape)
    l = os.listdir(os.curdir)
    l.sort()
    for x in l:
        if 'notes' not in x and '.txt' in x:
            a = RamanSpectrum(x)
            a = removespikes(a)
            data = append(data,array([a.values]), axis = 0)
            frequency = append(frequency, array([a.index]),axis = 0)
    print (data.shape[0], 'spectra averaged')
    #plot(frequency[:,512], 's')
   
    ##########now the data is in the proper form.
    
    tup = copy(data)
   
    pad = zeros((tup.shape[0],512))
    tup = append(pad,tup,axis = 1)
    tup = append(tup,pad,axis = 1)   #### pad with zeroes
    Snoise = array([mean(tup,axis = 0)]*tup.shape[0])

    tup2 = copy(tup)
    for i in range(tup.shape[0]):
        
        z = mean(tup[i])/mean(Snoise[0])
        
        tup[i]-=z*Snoise[0]
        tup[i] = roll(tup[i],i)  ###### correct the pixel offset
        tup2[i] = roll(tup2[i],i)  # rolling tup2 but not subtracting noise

    tup_av = mean(tup[:,512:1536],axis = 0)
    tup2_av = mean(tup2[:,512:1536],axis = 0)
    tup2_av-=8800
    
    tup_av = SGsmooth(arange(1024),tup_av, width = 11, order = 3)

    return RamanSpectrum(pandas.Series(tup_av,frequency[0]))
    
def fitspectrum(spectrum, rnge, functiontype, guess, function=None):
    numpeaks = int((len(guess)-2)/3)
    print(numpeaks)
    if functiontype == 'Custom':
        pass
    elif functiontype=='xGaussian':
        numpeaks = int((len(guess)-2)/3)
        def function(x,*guess):
            
            y = guess[-2]*x/1000+guess[-1]
            for i in range(numpeaks):
                y+= guess[i]*exp(-(x-guess[i+numpeaks])**2/guess[i+2*numpeaks])
            return y
    elif functiontype=='xGaussianNoBase':
        numpeaks = (len(guess))/3
        def function(x,*guess):
            
            y = np.zeros(x.shape)
            for i in range(numpeaks):
                y+= guess[i]*exp(-(x-guess[i+numpeaks])**2/guess[i+2*numpeaks])
            return y
    elif functiontype=='xGaussianExpBase':
        numpeaks = (len(guess)-3)/3
        def function(x,*guess):
            with np.errstate(divide='ignore'):
                y = exp(guess[-2]*(x-guess[-1]))+guess[-3]
                
                for i in range(numpeaks):
                    y+= guess[i]*exp(-(x-guess[i+numpeaks])**2/guess[i+2*numpeaks])
            return y
    elif functiontype=='xGaussianFix':
        numpeaks = (len(guess))/3
        x0list = guess[numpeaks:2*numpeaks]
        def function(x,*guess):
            
            y = np.zeros(x.shape)
            for i in range(numpeaks):
                y+= guess[i]*exp(-(x-x0list[i])**2/guess[i+2*numpeaks])
                
            return y

    elif functiontype=='xVoigt':
        numpeaks = (len(guess))/4   
        def function(x,*guess):
            y = np.zeros(x.shape)
            for i in range(numpeaks):
                A  = guess[i] 
                x0 = guess[i+numpeaks]
                G = guess[i+2*numpeaks]
                g= guess[i+3*numpeaks]
                y = np.zeros(x.shape)
                y+= A*exp(-(x-x0)**2/G)*g/(pi*((x-x0)**2+g**2))
            return y
    else: 
        print ("function not understood")
        def function(x,m,b):return m*x/1000+b
    start = spectrum.nearest(rnge[0])
    end = spectrum.nearest(rnge[1])
    if start>=end:
        print ("error with range of fit")
        return -1
    x = array(spectrum.index[start:end])
    y = spectrum.values[start:end]
    y_in = array(spectrum.values[start:end])
    try:
        result = scipy.optimize.curve_fit(function,x,y,guess,sigma=True)
    except RuntimeError:
        result = [np.ndarray((len(guess),))*np.nan,[0],[0]]
        print ('problem with fitting')
        
    
    fitresult = namedtuple('fitresult', 'params x y_in y peaks areas report')
    peaks = list()
    areas = list()
    listofvariables = str()
   
    listofvariables = functiontype
    if functiontype == "Custom":
        listofvariables = strnum(inspect.getargspec(function)[0][1:])
    elif 'Gaussian' in functiontype:
        order = np.argsort(result[0][numpeaks:2*numpeaks])
        
        result[0][:numpeaks] = result[0][order] 
        result[0][numpeaks:2*numpeaks] = result[0][numpeaks+order] 
        result[0][2*numpeaks:3*numpeaks] = result[0][2*numpeaks+order] 
        for i in range(numpeaks):
            A = result[0][i]
            x0 = result[0][i+numpeaks]
            G = result[0][i+2*numpeaks]
            
            
            if 'NoBase' in functiontype:
                y= A*exp(-(x-x0)**2/G)
            elif 'Fix' in functiontype:
                A = result[0][i]
                G = result[0][i+numpeaks]
                y= A*exp(-(x-x0list[i])**2/G)+result[0][-2]/1000*x+result[0][-1]
            elif 'Exp' in functiontype:
                y= A*exp(-(x-x0)**2/G)+exp(result[0][-2]*(x-result[0][-1]))+result[0][-3]
            else:
                y= A*exp(-(x-x0)**2/G)+result[0][-2]/1000*x+result[0][-1]
            ar1 = A*np.sqrt(pi*G)
            areas.append(ar1)
            peaks.append(y)
    elif 'Lorentzian' in functiontype:
        for i in range(numpeaks):
            if 'NoBase' in functiontype:
                y= result[0][i]**2/((x-result[0][i+numpeaks])**2+result[0][i+2*numpeaks]**2)
            else:
                y= result[0][i]**2/((x-result[0][i+numpeaks])**2+result[0][i+2*numpeaks]**2)+result[0][-2]/1000*x+result[0][-1]
            peaks.append(y)
            areas.append(result[0][i]**2*pi/result[0][i+2*numpeaks])
    elif 'Voigt' in functiontype:
        for i in range(numpeaks):
            A = result[0][i]
            x0 = result[0][i+numpeaks]
            G = result[0][i+2*numpeaks]
            g = result[0][i+3*numpeaks]
            
            if 'NoBase' in functiontype:
                y= A*exp(-(x-x0)**2/G)*g/(pi*((x-x0)**2+g**2))
            else:
                y= A*exp(-(x-x0)**2/G)*g/(pi*((x-x0)**2+g**2))+result[0][-2]/1000*x+result[0][-1]
            peaks.append(y)
            ar1 = A*np.sqrt(pi*G)
            areas.append(ar1)
    
    report=str(spectrum.name)+','+listofvariables+','+strnum(result[0])+'\n'
    pt1 = fitresult(result,x,y_in,function(x,*result[0]),peaks,areas,report)
    
    return pt1
    





###################
##create wavelength dependent correction for SPID RamanMicroscope.  Then remove these artifacts from spectrum.

    
def quickoffset(ax,rnge=None,offset = None,autolim=True):
    offset_val = 0
    if offset is None:
        for i in range(1,len(ax.lines)):
            
            if rnge==None:
                start = 0
                end = -1
            else:
                start = np.argmin(abs(rnge[0]-ax.lines[i].get_xdata()))
                end = np.argmin(abs(rnge[1]-ax.lines[i].get_xdata()))
            
            y = ax.lines[i].get_ydata()[start:end]
            ax.lines[i].set_ydata(ax.lines[i].get_ydata()-min(y))
            y = ax.lines[i].get_ydata()[start:end]
            offset_val=max(offset_val,max(y))
    else:
        offset_val = offset
    for i in range(len(ax.lines)):
        ax.lines[i].set_ydata(ax.lines[i].get_ydata()+(len(ax.lines)-i-1)*offset_val)
    if autolim:

        ax.set_xlim(rnge[0],rnge[1])
        ax.autoscale(enable=True)
        ax.relim()
        ax.autoscale_view(scalex=False)
        plt.draw()
        
    return 0
    
def smoothbaseline(spectrum, rnge1,rnge2,_plot= False,ax=None):
    start = spectrum.nearest(rnge1[0])
    end = spectrum.nearest(rnge1[-1])
    
 
    if start == end:
        print ('error.  same start and end points')
        return spectrum
    xs = array(spectrum.index[start:end])
    
    ys= spectrum.values[start:end]
    x1 = xs[int(len(xs)/2)]
    y1 = ys[int(len(ys)/2)]
    [k1,intercept1] = list(np.polynomial.polynomial.polyfit(xs,ys,1))
    
    start = spectrum.nearest(rnge2[0])
    end = spectrum.nearest(rnge2[-1])
    
 
    if start == end:
        print ('error.  same start and end points')
        return spectrum
    xs = array(spectrum.index[start:end])
    ys= spectrum.values[start:end]
    x2 = xs[int(len(xs)/2)]
    y2 = ys[int(len(ys)/2)]
    [k2,intercept2] = list(np.polynomial.polynomial.polyfit(xs,ys,1))
    
    
#    a3 = ((k1+k2)*(x1-x2)-2*y1+2*y2)/(x1-x2)**3
#    a2 = (-k1*(x1-x2)*(x1+2*x2)+k2*(-2*x1**2+x1*x2+x2**2)+3*(x1+x2)*(y1-y2))/(x1-x2)**3
#    a1 = (k2*x1*(x1-x2)*(x1+2*x2)-x2*(k1*(-2*x1**2+x1*x2+x2**2)+6*x1*(y1-y2)))/(x1-x2)**3
#    a0 = (x2*(x1*(-x1+x2)*(k2*x1+k1*x2)-x2*(-3*x1+x2)*y1)+x1**2*(x1-3*x2)*y2)/(x1-x2)**3    
#    backup
    
    a3 = ((k1+k2)*(x1-x2)-2*y1+2*y2)/(x1-x2)**3
    a2 = (-k1*(x1-x2)*(x1+2*x2)+k2*(-2*x1**2+x1*x2+x2**2)+3*(x1+x2)*(y1-y2))/(x1-x2)**3
    a1 = (k2*x1*(x1-x2)*(x1+2*x2)-x2*(k1*(-2*x1**2+x1*x2+x2**2)+6*x1*(y1-y2)))/(x1-x2)**3
    a0 = (x2*(x1*(-x1+x2)*(k2*x1+k1*x2)-x2*(-3*x1+x2)*y1)+x1**2*(x1-3*x2)*y2)/(x1-x2)**3
    r = array([a3,a2,a1,a0])
    
    start = spectrum.nearest(x1)
    end = spectrum.nearest(x2)
    xs = array(spectrum.index[start:end]) 
    
    
    spectrum.values[start:end]=spectrum.values[start:end]-np.polynomial.polynomial.polyval(r,xs)
    if _plot:
        
        if ax==None:
            plt.plot(xs,np.polynomial.polynomial.polyval(r,xs))
            plt.plot([x1],[y1],'s')
        else:
            ax.plot(xs,np.polynomial.polynomial.polyval(r,xs),linewidth=1)
           
            
        
    
    return spectrum
    
    
    
def inflection_point(spectrum, rnge):
    start = spectrum.nearest(rnge[0])
    end = spectrum.nearest(rnge[1])
    xs = array(spectrum.index[start:end])
    ys=spectrum.values[start:end]
    p = numpy.polynomial.polynomial.polyfit(xs,ys,5)
    second_derivative = numpy.polynomial.polynomial.polyder(p,m=2)
    return xs[argmin(numpy.polynomial.polynomial.polyval(xs,second_derivative)**2)]
    
    

    
def save_parameters_from_fit(filename, fitresult):
    paramarray = copy(fitresult.params[0])
    numpeaks = (paramarray.size-2)/3
    np.savetxt(filename, fitresult.params[0], delimiter=0)
       
    return 0

def strnum(arr):
    out = str(arr[0])
    for i in arr[1:]:
        out+=','
        out+=str(i)
    return out