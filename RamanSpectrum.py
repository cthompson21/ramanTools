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
from copy import copy, deepcopy
import scipy.optimize
import matplotlib.pyplot as plt
from collections import namedtuple
import inspect
import os
import numpy.polynomial.polynomial
import struct

import matplotlib


"""
### An example of how to use the fitting function
a = RamanSpectrum(i, filetype='FTIR GVD')
a.autobaseline((530,675,1367,1480,1830,1920), join='start',specialoption='points', order = 4)
a.autobaseline((1920,2100,2500,3900,4000), join='start',specialoption='points', order = 3)
a-=min(a[850:1400])
a.plot()

z=fitspectrum(a, (950, 1300), 'xGaussian',guess )
print(z.printable_report)
z.plot(label=i, color = 'k', linestyle='--')#_peaks()
z.plot_peaks(linewidth=0.5,color = 'k')
savefig(i[:-4]+'.png',dpi=180)   
        
"""

class RamanSpectrum(pandas.Series):
    
    def __init__(self,fname_or_other,avg = True, filetype='unknown'):
        
        
        if filetype == 'Raman GVD':
            
            a=np.loadtxt(fname_or_other,unpack=True, delimiter = '\t', skiprows=0, usecols=(0,1))
            pandas.Series.__init__(self,a[1],a[0],dtype = float)
            self.name = fname_or_other
            
            self.filetype = filetype
            self.comment=""
        elif filetype == 'ATR-FTIR':  
            a=np.loadtxt(fname_or_other,unpack=True, delimiter = ',', usecols=(0,1))
            a=a[:,::-1]  ### reverse values 
            pandas.Series.__init__(self,a[1],a[0],dtype = float)
            self.name = fname_or_other
            self.filetype = filetype
            self.comment=""
        elif filetype == 'FTIR GVD':
            if fname_or_other[-4:].lower() == '.csv':
                a=np.loadtxt(fname_or_other,unpack=True, delimiter = ',', skiprows=2, usecols=(0,1))
                a=a[:,::-1]  ### reverse values 
                a-=np.min(a)
                pandas.Series.__init__(self,a[1],a[0],dtype = float)
                self.name = fname_or_other
                self.filetype = filetype
                self.comment=""
                self._name=self.name
            elif fname_or_other[-4:].lower() == '.asc':
                a=np.loadtxt(fname_or_other,unpack=True, delimiter = '\t', skiprows=26, usecols=(0,1))
                a=a[:,::-1]  ### reverse values 
                a-=np.min(a)
                pandas.Series.__init__(self,a[1],a[0],dtype = float)
                self.name = fname_or_other
                self.filetype = 'FTIR ASC'
                self.comment=""
                self._name=self.name
            
        elif filetype == 'ATR-FTIR GVD':
            if fname_or_other[-4:].lower() == '.csv':
                a=np.loadtxt(fname_or_other,unpack=True, delimiter = ',', skiprows=0, usecols=(0,1))
                
                a-=np.min(a)
                pandas.Series.__init__(self,a[1],a[0],dtype = float)
                self.name = fname_or_other
                self.filetype = filetype
                self.comment=""
                self._name=self.name
                
        elif filetype == 'JCAMP':  
            a=np.array(readJCAMP(fname_or_other))
            a=a[:,::-1]  ### reverse values 
            pandas.Series.__init__(self,a[1],a[0],dtype = float)
            self.name = fname_or_other
            self.filetype = filetype
            self.comment=""
        elif filetype == 'FTIR ASC':  
            a=np.loadtxt(fname_or_other,unpack=True, delimiter = '\t', skiprows=26, usecols=(0,1))
            a=a[:,::-1]  ### reverse values 
            a-=np.min(a)
            pandas.Series.__init__(self,a[1],a[0],dtype = float)
            self.name = fname_or_other
            self.filetype = 'FTIR ASC'
            self.comment=""
       
        elif filetype=='pandas.Series':# == pandas.Series:
            pandas.Series.__init__(self,fname_or_other.values,array(fname_or_other.index))
            
            self.name = ""
            self.comment=""
            self.info = ""
            self.num_frames = 0
            self.accum_time = 0
            self.filetype ='manual input'
            
        elif fname_or_other[-4:] == '.jdx':
            a=np.array(readJCAMP(fname_or_other))
            a=a[:,::-1]  ### reverse values 
            pandas.Series.__init__(self,a[1],a[0],dtype = float)
            self.name = fname_or_other
            self.filetype = filetype
            self.comment=""
            
        elif fname_or_other[-4:] == '.SPE' or fname_or_other[-4:] == '.spe':
            
            
            fid = File(fname_or_other)
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
            a.sort(ascending=True)
            a = polyeval(a,np.arange(len(spectrum))+x_offset)
            a = -(1E7/a-1E7/wl)
            pandas.Series.__init__(self,spectrum,a)
            
            self.comment=""
            self.info = fid.get_info()
            self.accum_time=fid.get_accum_time()
            self.num_frames = fid._numframes
            self.name = fname_or_other
            self.filetype ='SPE Raman'
            fid.close()

        elif fname_or_other[-4:] == '.SSM':
            a = loadtxt(fname_or_other,delimiter ='  ', unpack = True,skiprows =2)
            pandas.Series.__init__(self,a[1],a[0],dtype = float)#pandas.Series.__init__(self,a[1][0:cutoff],a[0][0:cutoff],dtype = float)

            self.filetype ='SSM UVVis'
            self.comment=""
            self.info = ""
            self.accum_time=0
            self.num_frames = 0
            self.name = fname_or_other
            
        elif fname_or_other[-4:] == '.txt' or fname_or_other[-4:] == '.TXT':
            try:   
                a=np.loadtxt(fname_or_other,unpack=True, delimiter = ',', skiprows=0, usecols=(0,1))
            except:
                print('error in opening text file')
                raise
            pandas.Series.__init__(self,a[1],a[0],dtype = float)
            self.filetype = 'txt'
            self.comment=""
            self.info = ""
            self.accum_time=0
            self.num_frames = 0
            self.name = fname_or_other

        elif fname_or_other[-4:]=='.xlsx':
            a=np.loadtxt(fname_or_other,unpack=True, delimiter = ',', skiprows=2, usecols=(0,1))
            pandas.Series.__init__(self,a[1],a[0],dtype = float)
            self.name = fname_or_other
            self.filetype = filetype
            self.comment=""
            
            
        elif fname_or_other[-3:]=='.ch':
            with open(fname_or_other,'rb') as f:
                binary = f.read()
                f.close()
            z=array(struct.unpack('1559h',binary[1025:4143]))
            to_remove= np.intersect1d(where(diff(z)>3000)[0]+1,where(diff(z)<-3000)[0])
            z = numpy.delete(z,to_remove)
            chromatogram = array(list(sum(z[:i]) for i in range(len(z))), dtype=float)
            chromatogram*=2/4096
            time_len = len(chromatogram)*0.16666/26#(to_remove[1]-to_remove[0])#(len(to_remove)+1)/6 ## in minutes
            
            time = linspace(-0.042,time_len-0.042,len(chromatogram))

            pandas.Series.__init__(self,chromatogram,time,dtype = float)
            self.name = struct.unpack('7s',binary[25:32])[0].decode('utf-8').replace('\x00','')#str(binary[25:32])#str(struct.unpack('7c',binary[25:32]))
            self.wavelength=struct.unpack('11s',binary[604:615])[0].decode('utf-8').replace('\x00','')# str(struct.unpack('7s',binary[604:611]))
#        self._name = self.name
        return 
#     
#    def __array_finalize__(self,obj):
#        if obj is None: return
            
    def apply(self,func,):
        self.values[:]=  func(self.values[:])
   

            
    def newadd(self,other,mode = 'union'):
        
        if mode == 'union':
            xmax = min([np.max(array(self.index)),np.max(array((other.index)))])
        
            xmin = max([min(array(self.index)),np.min(array(other.index))])
         
            _x2 = self.truncate(xmin,xmax)
            _y2 = other.truncate(xmin,xmax)
            _y2 = _y2.reindex(_x2.index,method = 'ffill')
        elif mode == 'join':
            _y2 = _y2.reindex(self.index,method = 'ffill')

        output = _x2+_y2
        if any(output.isnull()):
            output = output.fillna(method = 'bfill')
        
        return RamanSpectrum(output)

        
    def reverse(self):
        """don't use this anymore -- use 'sort'"""
        return RamanSpectrum(pandas.Series(self.values[::-1],self.index[::-1]))
    def addoffset(self,offset):
        self.loc[:]=self.loc[:]+offset
        return self#RamanSpectrum(pandas.Series(self.values,array(self.index)+offset))
    def nearest(self, x):
        if x in array(self.index):
            return np.argmin(abs(array(self.index)-x))
        if x<self.index[0]:
            return 0
        elif x>self.index[-1]:
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
#                return self.get_values()[self.nearest(key)]
                return self.values[self.nearest(key)]
            elif isinstance(key,list):
                return array([self[i] for i in key])
            elif isinstance(key,np.ndarray):
                return array([self[i] for i in key.astype(list)])
            else:
                
                return []
                
                
        except:
            raise
        return -1
        
    def findminbetween(self,rnge):
        return findminbetween(self,rnge)
#
#    def copy(self):
#        return copy(self,deep=deep)

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
    def derivative(self,order):
        
        w = 11
        a = np.array([nan,]*w)
        for x in range(w,len(self.values)-w+1):
            xs = np.array(self.index)[x-w:x+w]
            ys = self.values[x-w:x+w]
            r = np.polynomial.polynomial.polyfit(xs,ys, order)
            a=np.append(a,2*r[2])
        a = np.append(a,np.array([nan,]*(w-1)))
        self.values[:] = a
        return 

class ramansubplot(matplotlib.axes.Axes):

    def __init__(self,n):
        obj = matplotlib.axes._subplots.AxesSubplot.__init__(self,n)
        self.set_xlabel('Raman Shift')
        self.set_ylabel('Intensity')
        return obj
    
    
    
    
    
    
    
    
        
def SGsmooth(x,y, width=11,order = 2):#data,rnge,rnge_type = 'data'):
    
    retval = np.ndarray(y.shape)

    for i in range(y.size):
        i_min = int(max(0,i-width/2))
        i_max = int(min( i+width/2-1, y.size))
        fit = np.polyfit(x[i_min:i_max+1],y[i_min:i_max+1],order)
        retval[i] = polyeval(fit,x[i])
    return retval

def secondderivative(x,y, width=11,order = 2):#data,rnge,rnge_type = 'data'):
    
    retval = np.ndarray(y.shape)

    for i in range(y.size):
        i_min = int(max(0,i-width/2))
        i_max = int(min( i+width/2-1, y.size))
        fit = np.polyfit(x[i_min:i_max+1],y[i_min:i_max+1],order)
        retval[i] = fit[0]
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
    if len(xs)==0 or len(ys)==0:
        print('empty range data')
    if specialoption=='points':
        #pdb.set_trace()
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
                spectrum.iloc[start:] -= offset
        elif join == 'end':
            if end == spectrum.size-1:
                print ('end point too far')
            else:
                offset = spectrum.iloc[end]-spectrum.iloc[end-1]
                spectrum.iloc[end:]-= offset
        else:
            print( 'wrong join option in autobaseline:', join)
        return spectrum
    elif partial_fit_region!=None:
        lo_start = spectrum.nearest(partial_fit_region[0])
        lo_end = spectrum.nearest(partial_fit_region[1])

        xfits = array(spectrum.index[start:end])
        yfits= spectrum.values[start:end]
        r = np.polyfit(xfits,yfits,order)
        spectrum.values[lo_start:lo_end] = spectrum.values[lo_start:lo_end]-polyeval(r,array(spectrum.index[lo_start:lo_end]))
        spectrum.values[lo_start:lo_end]-=np.min(spectrum.values[start:end])
        
    else:
        xfits = xs
        yfits = ys

    if order == 0:
        slope =(yfits[-1]-yfits[0])/(xfits[-1]-xfits[0])
        b = yfits[0]-slope*xfits[0]
        baseline = slope*(xs)+b
        spectrum.iloc[start:end] = spectrum.iloc[start:end]-baseline
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
            spectrum.iloc[start:] -= offset
    elif join == 'end':
        if end == spectrum.size-1:
            print ('end point too far')
        else:
            offset = spectrum.iloc[end]-spectrum.iloc[end-1]
            spectrum.iloc[end:]-= offset
    else :
        print ('wrong join option in autobaseline:', join)
    spectrum[:]-=min(spectrum.values)
    return spectrum
    


#def smooth(spectrum,window_len=5,window='flat'):
#    
#    spectrum =spectrum.copy()
#    if spectrum.ndim != 1:
#        raise ValueError ("smooth only accepts 1 dimension arrays.")
#
#    if spectrum.size < window_len:
#        raise ValueError("Input vector needs to be bigger than window size.")
#
#
#    if window_len<3:
#        return spectrum
# 
#    if not window in ['flat','SG']:
#        raise ValueError( "Window is on of 'flat', 'SG")
#
#
#    s=np.r_[spectrum.values[window_len-1:0:-1],spectrum,spectrum.values[-1:-window_len:-1]]
#   
#    if window == 'SG':
#        
#        order = 2
#        y=spectrum.values
#        x=np.array(spectrum.index)
#        retval = np.ndarray(y.shape)
#        for i in range(y.size):
#            
#            i_min = max(0,i-window_len/2)
#            i_max = min( i+window_len/2-1, y.size)
#            fit = np.polynomial.polynomial.polyfit(x[i_min:i_max+1],y[i_min:i_max+1],min(order,i_max-i_min))
#            retval[i] = polyeval(fit,x[i])
#        
#
#        spectrum[:]=retval[:]
#        return spectrum
#        
#            
#    elif window == 'flat': #moving average
#        w=np.ones(window_len,'d')
#        spectrum.values[:] =np.convolve(w/w.sum(),s,mode='valid')[(window_len-1)/2:-(window_len-1)/2]
#    else:
#        w=eval(window+'(window_len)')
#        spectrum.values[:] =np.convolve(w/w.sum(),s,mode='valid')[(window_len-1)/2:-(window_len-1)/2]
#   
#    
#    return spectrum 

def smooth(spectrum,window_len=5,window='flat'):
    if spectrum.ndim != 1:
        raise ValueError ("smooth only accepts 1 dimension arrays.")

    if spectrum.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")
    if window_len<3:
        return spectrum
    if not window in ['flat','SG']:
        raise ValueError( "Window is on of 'flat', 'SG")
    if window == 'SG':
        order = 2
        if order>window_len-1:
            print('under-parametrized smoothing.  Window_len must be greater than order.')
        y=spectrum[:]
        x=np.array(spectrum.index)
        rval = np.ndarray(y.shape)
        for i in range(y.size):
            
            i_min = int(max(0,i-(window_len-1)/2))
            i_max = int(min(i+(window_len-1)/2, y.size))
            fit = np.polyfit(x[i_min:i_max+1],spectrum.iloc[i_min:i_max+1],order)
            rval[i] = np.polyval(fit,x[i])
        spectrum[:]=rval[:]
        return 
     
    elif window == 'flat': #moving average
        w=np.ones(window_len,'d')
        
        s=np.r_[spectrum.iloc[window_len-1:0:-1],spectrum.values[:],spectrum.values[-1:-window_len:-1]]        
        spectrum[:]=np.convolve(w/w.sum(),s,mode='valid')[int((window_len-1)/2):int(-(window_len-1)/2)]
    else:
        w=eval(window+'(window_len)')
        s=np.r_[spectrum.values[window_len-1:0:-1],spectrum.values[:],spectrum.values[-1:-window_len:-1]]        
        spectrum[:]=np.convolve(w/w.sum(),s,mode='valid')[int((window_len-1)/2):int(-(window_len-1)/2)]
    return
    
def calc_noise(spectrum,rnge,rnge_type = 'data'):
    
    start = spectrum.nearest(rnge[0])
    end = spectrum.nearest(rnge[-1])
    
   
    xs = array(spectrum.index[start:end])
    ys= spectrum.values[start:end]
  
    r = np.polynomial.polynomial.polyfit(xs,ys,2)
    average = np.polynomial.polynomial.polyval(xs, r)
    
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



    
    
def ics(orca_output_file, normalize = False,color='k',labelpeaks = True,x_factor=1):
    fileopen = open(orca_output_file,'rt')
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
                plt.gca().annotate(i[0], (float(i[1])*x_factor, float(i[2])+0.2), color=color,fontsize = 8,horizontalalignment='center') 
        elif not i[0][0].isdigit():
            break
    table = transpose(table)
    if normalize == True:
            table[1]/=max(table[1])
    return plt.vlines(table[0]*x_factor,0,table[1],linewidth = 2,color=color)
    
def irs(orca_output_file, normalize = False,color='k',labelpeaks = True, x_factor=1):

    fileopen = open(orca_output_file,'rt')
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
                plt.gca().annotate(i[0], (float(i[1])*x_factor, float(i[2])+0.2), color=color,fontsize = 8,horizontalalignment='center') 
        elif not i[0][0].isdigit():
            
            break
    table = transpose(table)
    if normalize == True:
            table[1]/=max(table[1])
    return plt.vlines(table[0]*x_factor,0,table[1],linewidth = 2,color=color)
    

    


class fitspectrum(pandas.Series):
    
    def __init__(self,spectrum, rnge, functiontype, guess, function=None):
        
        
        start = spectrum.nearest(rnge[0])
        end = spectrum.nearest(rnge[1])
        if start>=end:
            print ("error with range of fit")
            return -1
        x = array(spectrum.index[start:end])
        y = array(spectrum.values[start:end])
        y_in = array(spectrum.values[start:end])
        
        pandas.Series.__init__(self,np.zeros(x.shape),x)
        self.bounds=(-inf,inf)
        self.guess=guess
        self.x0list = []
        self.functiontype = functiontype
        
        if functiontype == 'Custom':
            pass
        elif functiontype=='xGaussian':
            if len(guess)%3 !=2:
                print('bad guess detected')
                return 
            numpeaks = int((len(guess)-2)/3)
            self.set_default_bounds(numpeaks)
            def function(x,*guess):
                """guess in the format A1, A2,...,x1, x2,..., G1, G2..., m,b"""
                y = guess[-2]*x/1000+guess[-1]
                for i in range(numpeaks):
                    y+= guess[i]*exp(-(x-guess[i+numpeaks])**2/guess[i+2*numpeaks])
                
                return y
            
        elif functiontype=='xGaussianPoly':
            if len(guess)%3 !=0:
                print('bad guess detected', len(guess), len(guess)%3)
                return 
            numpeaks = int((len(guess)-3)/3)
            def function(x,*guess):
                """guess in the format A1, A2,...,x1, x2,..., G1, G2..., m,b"""
                y = guess[-3]*x**2/1e6 + guess[-2]*x/1e3 + guess[-1]
                for i in range(numpeaks):
                    y+= guess[i]*exp(-(x-guess[i+numpeaks])**2/guess[i+2*numpeaks])
                return y
            
        elif functiontype=='xGaussianNoBase':
            if len(guess)%3 !=0:
                print('bad guess detected')
                return 
            numpeaks = int(len(guess)/3)
            self.set_default_bounds(numpeaks)
            def function(x,*guess):
                
                y = np.zeros(x.shape)
                for i in range(numpeaks):
                    y+= guess[i]*exp(-(x-guess[i+numpeaks])**2/guess[i+2*numpeaks])
                
                    
                return y
        elif functiontype=='xGaussianExpBase':
            if len(guess)%3 !=3:
                print('bad guess detected')
                return 
            numpeaks = (len(guess)-3)/3
            def function(x,*guess):
                with np.errstate(divide='ignore'):
                    y = exp(guess[-2]*(x-guess[-1]))+guess[-3]
                    
                    for i in range(numpeaks):
                        y+= guess[i]*exp(-(x-guess[i+numpeaks])**2/guess[i+2*numpeaks])
                return y
        elif functiontype=='xGaussianFix':
            if len(guess)%3 !=0:
                print('bad guess detected')
                return 
            numpeaks = int((len(guess))/3)
            self.set_default_bounds(numpeaks)
            self.x0list = list(guess[slice(numpeaks,2*numpeaks)])
            guess = guess[:numpeaks]+guess[2*numpeaks:]
            
            def function(x,*guess):
                y = np.zeros(x.shape)
                for i in range(numpeaks):
                    y+= guess[i]*exp(-(x-self.x0list[i])**2/guess[i+numpeaks])   
                return y
            

        elif functiontype=='xVoigt':
            if len(guess)%4 !=0:
                print('bad guess detected')
                return 
            numpeaks = int((len(guess))/4)   
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
        elif functiontype =='xLorentzian':
            if len(guess)%3 !=2:
                print('bad guess detected')
                return 
            numpeaks = int((len(guess))/3)   
            def function(x,*guess):
                
                y=guess[-2]/1000*x+guess[-1]
                for i in range(numpeaks):
                    A  = guess[i] 
                    x0 = guess[i+numpeaks]
                    G = guess[i+2*numpeaks]
                    y+= A**2/((x-x0)**2+G**2)
                
                return y
        elif functiontype=='xGaussianLimBase':
            if len(guess)%3 !=1:
                print('bad guess detected')
                return 
            numpeaks = int(len(guess)/3)
            lowerBounds = np.zeros(numpeaks*3)
            upperBounds = np.ones(numpeaks*3)*1e6
            for i in range(numpeaks):
                lowerBounds[i+numpeaks] = guess[i+numpeaks]-guess[-1]
                upperBounds[i+numpeaks] = guess[i+numpeaks]+guess[-1]
            self.bounds = (lowerBounds, upperBounds)
            guess = guess[:-1]
            def function(x,*guess):
                """guess in the format A1, A2,...,x1, x2,..., G1, G2..., bound"""
                y = np.zeros(x.shape)
                for i in range(numpeaks):
                    y+= guess[i]*exp(-(x-guess[i+numpeaks])**2/guess[i+2*numpeaks])
                    
                return y
        elif functiontype=='xGaussianLimBase_SiArtifact':
            if len(guess)%3 !=2:
                print('WARNING: bad guess detected')
                return 
            numpeaks = int(len(guess)/3)
            lowerBounds = np.zeros(numpeaks*3+1)
            upperBounds = np.ones(numpeaks*3+1)*1e6
            for i in range(numpeaks):
                lowerBounds[i+numpeaks] = guess[i+numpeaks]-guess[-2]
                upperBounds[i+numpeaks] = guess[i+numpeaks]+guess[-2]
            lowerBounds[-1] = -1
            self.bounds = (lowerBounds, upperBounds)
            guess = guess[:-1]
            def function(x,*guess):
                """guess in the format A1, A2,...,x1, x2,..., G1, G2..., bound"""
                y = np.zeros(x.shape)+guess[-1]*exp(-(x-1107)**2/360)
                for i in range(numpeaks):
                    y+= guess[i]*exp(-(x-guess[i+numpeaks])**2/guess[i+2*numpeaks])
                    
                return y

            
        elif functiontype=='OxideExilisDeconv_SiArtifact':
            if len(guess)%3 !=1:
                print('bad guess detected:continuing')
                
            numpeaks = int(len(guess)/3)
            lowerBounds = [0,0,0,0,0,0,0,0,
                 960,1003,1075,1060,1164,1180,1260,1271,
                 50,500,4000,2000,500,500,20,20,0]
            upperBounds =  np.array([2,2,2,2,2,2,2,2,
                 960,1003,1075,1060,1164,1280,1260,1271,
                 150,1500,6000,4000,1500,1500,80,80,1])
            for i in range(numpeaks):
                lowerBounds[i+numpeaks] = guess[i+numpeaks]-10
                upperBounds[i+numpeaks] = guess[i+numpeaks]+10
            lowerBounds[-1] = -1
            self.bounds = (lowerBounds, upperBounds)
            
            def function(x,*guess):
                """guess in the format A1, A2,...,x1, x2,..., G1, G2..., bound"""
                y = np.zeros(x.shape)+guess[-1]*exp(-(x-1107)**2/360)
                for i in range(numpeaks):
                    y+= guess[i]*exp(-(x-guess[i+numpeaks])**2/guess[i+2*numpeaks])
                    
                return y
        else: 
            print ("function not understood")
            def function(x,m,b):return m*x/1000+b
      
        listofvariables = str()
       
        listofvariables = functiontype
        
        
        
        try:
            self.numpeaks = int(numpeaks)
        except:
            raise
        self.function=function
        self.printable_report ='                '+spectrum.name+'            \n'+"{0:15s}{1:15s}{2:15s}{3:15s}{4:15s}\n".format('A', 'x0','G','FWHM','Area')
     
        
        self.params = []
        self.y_in=np.array([])
        self.peaks = []
        self.areas = []
        self.report = ''
        
        
        try:
            self.result = scipy.optimize.curve_fit(function,x,y,guess,sigma=None, bounds = self.bounds)
            self.values[:] = self.function(array(self.index), *self.result[0])
            
        except RuntimeError:
            print('runtime error')
            self.result=np.array([guess])
            self.result[:self.numpeaks]=0
            self.valid=-1
            print(self.result,self.numpeaks)
        try:    
            self.compute_peaks(x)
        except:
            raise
            
                

        
                
        
        self.params = self.result[0]
        self.y_in=y_in
     
        self.report = str(spectrum.name)+','+listofvariables+','+strnum(self.result[0])+'\n'
        
            
        return 
    def compute_peaks(self,x,):
        
        if self.functiontype == "Custom":
            listofvariables = strnum(inspect.getargspec(function)[0][1:])
            
        elif self.functiontype == 'xGaussianFix':
          
            for i in range(self.numpeaks):
                A = self.result[0][i]
                x0 = self.x0list[i]
                G = self.result[0][i+self.numpeaks]
                y= A*np.exp(-(x-x0)**2/G)
                ar1 = A*np.sqrt(pi*G)
                self.areas.append(ar1)
                self.peaks.append(y)
            
#            newres= list(self.result[0][0:int(len(self.result[0])/2)])+self.x0list+list(self.result[0][int(len(self.result[0])/2):])
            newres= list(self.result[0][0:self.numpeaks])+self.x0list+list(self.result[0][self.numpeaks:])
            self.result=[list(newres)]

                
                
        elif 'Gaussian' in self.functiontype:
            order = np.argsort(self.result[0][self.numpeaks:2*self.numpeaks])
            self.result[0][:self.numpeaks] = self.result[0][order] 
            self.result[0][self.numpeaks:2*self.numpeaks] = self.result[0][self.numpeaks+order] 
            self.result[0][2*self.numpeaks:3*self.numpeaks] = self.result[0][2*self.numpeaks+order] 
            for i in range(self.numpeaks):
                A = self.result[0][i]
                x0 = self.result[0][i+self.numpeaks]
                G = self.result[0][i+2*self.numpeaks]

                if 'NoBase' in self.functiontype:
                    y= A*exp(-(x-x0)**2/G)

                elif 'Exp' in self.functiontype:
                    y= A*exp(-(x-x0)**2/G)+exp(result[0][-2]*(x-self.result[0][-1]))+self.result[0][-3]
                elif 'Poly' in self.functiontype:
                    y= A*exp(-(x-x0)**2/G)+self.result[0][-3]*x**2/1e6 + self.result[0][-2]*x/1e3 + self.result[0][-1]
                elif 'LimBase' in self.functiontype:
                    y= A*exp(-(x-x0)**2/G)
                else:
                    y= A*exp(-(x-x0)**2/G)+self.result[0][-2]/1000*x+self.result[0][-1]
                ar1 = A*np.sqrt(pi*G)
                self.areas.append(ar1)
                self.peaks.append(y)
                
        elif 'OxideExilisDeconv_SiArtifact' in self.functiontype:
            order = np.argsort(self.result[0][self.numpeaks:2*self.numpeaks])
            self.result[0][:self.numpeaks] = self.result[0][order] 
            self.result[0][self.numpeaks:2*self.numpeaks] = self.result[0][self.numpeaks+order] 
            self.result[0][2*self.numpeaks:3*self.numpeaks] = self.result[0][2*self.numpeaks+order] 
            for i in range(self.numpeaks):
                A = self.result[0][i]
                x0 = self.result[0][i+self.numpeaks]
                G = self.result[0][i+2*self.numpeaks]

                
                y= A*exp(-(x-x0)**2/G)
                ar1 = A*np.sqrt(pi*G)
                self.areas.append(ar1)
                self.peaks.append(y)
            
        elif 'Lorentzian' in self.functiontype:
            for i in range(self.numpeaks):
                if 'NoBase' in self.functiontype:
                    y= result[0][i]**2/((x-result[0][i+numpeaks])**2+result[0][i+2*numpeaks]**2)
                else:
                    y= self.result[0][i]**2/((x-self.result[0][i+self.numpeaks])**2+self.result[0][i+2*self.numpeaks]**2)+self.result[0][-2]/1000*x+self.result[0][-1]
                self.peaks.append(y)
                self.areas.append(self.result[0][i]**2*pi/result[0][i+2*self.numpeaks])
        elif 'Voigt' in self.functiontype:
            y=np.zeros(x.shape)
            for i in range(numpeaks):
                A = self.result[0][i]
                x0 = self.result[0][i+self.numpeaks]
                G = self.result[0][i+2*self.numpeaks]
                g = self.result[0][i+3*self.numpeaks]
                
                
                y+= A*exp(-(x-x0)**2/G)*g/(pi*((x-x0)**2+g**2))
                self.peaks.append(y)
                ar1 = A*np.sqrt(pi*G)
                self.areas.append(ar1)
                
        return 
    def plot_peaks(self, offset=0, **kwargs):
        for i in self.peaks:
            plt.plot(array(self.index), i+offset, **kwargs)
    def set_default_bounds(self,numpeaks):
        if self.functiontype == 'xGaussianNoBase':
            self.bounds =tuple((array([0,]*numpeaks+[0,]*numpeaks+[0,]*numpeaks),array([1e6]*numpeaks+[1e6]*numpeaks+[1e6]*numpeaks)))
        elif self.functiontype == 'xGaussian':
            self.bounds =tuple((array([0,]*numpeaks+[0,]*numpeaks+[0,]*numpeaks+[-100,-1000]),array([1e6]*numpeaks+[1e6]*numpeaks+[1e6]*numpeaks+[100,1000])))
        elif self.functiontype == 'xGaussianFix':
            self.bounds =tuple((array([0,]*numpeaks+[0,]*numpeaks),array([1e6]*numpeaks+[1e6]*numpeaks)))
        return
    def plot_guess(self):
        plt.plot(np.array(self.index), self.function(np.array(self.index), *self.guess))
        
    def print_report(self):
        self.printable_report ='                '+str(self.name)+'            \n'+"{0:15s}{1:15s}{2:15s}{3:15s}{4:15s}{5:15s}\n".format('A', 'x0','G','FWHM','Area',"Norm. Area")
        
        if 'Gaussian' in self.functiontype or 'conv' in self.functiontype:
            for i in range(self.numpeaks):
                try:
                    self.printable_report += "{0:0.3f}{1:13.1f}{2:15.3f}{3:15.3f}{4:15.3f}{5:15.3f}\n".format(*self.params[slice(i,self.numpeaks*2+i+1, self.numpeaks)], 
                                              2*sqrt(np.log(2)*self.params[self.numpeaks*2+i]), 
                                              self.areas[i], 
                                              self.areas[i]*100/np.max(self.areas))
                except:
                    self.printable_report += "error in peak report generation"
        self.printable_report += "-----------END REPORT---------------"
        print(self.printable_report)
        return 
    



    
def quickoffset(ax,rnge=None,offset = None,autolim=True, lines=False):
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
        if lines:
            plt.hlines((len(ax.lines)-i-1)*offset_val, ax.lines[i].get_xdata()[0],ax.lines[i].get_xdata()[-1],linewidth=0.5)
    if autolim:
        if rnge==None:
            rnge=(np.min(list(np.min(i.get_xdata()) for i in ax.lines)),np.max(list(np.max(i.get_xdata()) for i in ax.lines)))
        ax.set_xlim(rnge[0],rnge[1])
        ax.autoscale(enable=True, axis='y')
#        ax.relim()
#        ax.autoscale_view(scalex=False)
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

def readJCAMP(filename):
    startdata = 0
    deltax=1
    with open(filename, 'rt') as f:
        w=f.readlines()
        f.close()
    for i in range(len(w)):
        if "##DELTAX" in w[i]:
            deltax = float(w[i][9:15])
            continue
        if "##YUNITS" in w[i]:
            yunits = w[i][9:]
            
            continue
        if '##XYDATA' in w[i]:
            startdata = i
            break
    yfactor=0.001
       
        
    x=np.array([])
    y=np.array([])

    j0 = float(w[startdata+1].replace('#','').split(' ')[0])
    j1 = float(w[startdata+2].replace('#','').split(' ')[0])
    
    xdirection = np.sign(j1-j0)
    print(j0,j1,xdirection)
    for i in w[startdata:]:
        if 'END' in i:
            break
        i=i.replace('#','')
        j=i.split(' ')
        
        for z in arange(len(j)-1):
            try:
                x = np.append(x, float(j[0])+z*deltax*xdirection)
                y=np.append(y, float(j[1+z]))
            except:
                raise
            
    
    y[y<0.00001]=0.00001
    if xdirection==1:
        x=x[::-1]
        y=y[::-1]
    if 'TRANSMITTANCE' in yunits:
        
        y=-log10(y)
        y-=np.min(y) 
    
    return x,y

def readKelli(filename,direction = 1,yfactor=1):
    startdata = 0
    with open(filename, 'rt') as f:
        w=f.readlines()
        f.close()
    for i in range(len(w)):

#        if "##DELTAX" in w[i] or "##XFACTOR" in w[i]:
#            deltax = float(w[i][9:15])
#            continue
        
        if "##YUNITS" in w[i]:
            yunits = w[i][9:]
            
            continue
        if '##XYDATA' in w[i]:
            startdata = i
            break
    yfactor=0.001
    deltax=1
       
        
    x=np.array([])
    y=np.array([])

    j0 = float(w[startdata+1].replace('#','').split(' ')[0])
    j1 = float(w[startdata+2].replace('#','').split(' ')[0])
    
    xdirection = np.sign(j1-j0)
    xdirection=1
    print(j0,j1,xdirection)
    for i in w[startdata+1:]:
        if 'END' in i:
            break
        i=i.replace('#','')
        i=i.replace('\n', '')
        j=i.split(' ')
        print(j)
        
#        for z in arange(len(j)-1):
#            try:
#                x = np.append(x, float(j[0])+z*deltax*xdirection)
#                y=np.append(y, float(j[1+z]))
#            except:
#                raise
        if direction ==1:
            x=np.append(x, np.arange(float(j[0]),float(j[0])+len(j)-1,1.0))
            y=np.append(y, list((float(b) for b in j[1:])))
        elif direction==-1:
         
            x=np.append(x, np.linspace(float(j[0]),float(j[0])+40,11)[1:])
            y=np.append(y, list((float(b) for b in j[1:])))
    
    
    y*=yfactor
    if 'TRANSMITTANCE' in yunits:
        
        y=-log10(y)
        y-=np.min(y) 
    print(x[-1])
    return x,y


def interpolate():
    figure()
    x, y = readKelli("U:/temp/116-14-3-IR.jdx")
    np.savetxt('U:/temp/116-14-3-IR.csv', np.transpose([x, y]), delimiter = ',')
    x,y=readKelli("U:/temp/116-14-3-IR (1).jdx",direction=-1, yfactor = 0.000035919)

    newy = np.array([])
    for i in range(2,len(x)-2):
        a = np.polyfit(x[i-2:i+2],y[i-2:i+2], 2)
        xs=np.arange(x[i]-2,x[i]+2)
        newy=np.append(newy,np.polyval(a,xs))
        print(len(newy))
    newxs = np.arange(x[2],x[-2],1)  
    plot(newxs,newy,'k')
    np.savetxt('U:/temp/116-14-3-IR (1).csv', np.transpose([newxs, newy]), delimiter = ',')
    
    
    
    