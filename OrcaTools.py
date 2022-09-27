# -*- coding: utf-8 -*-
"""
Created on Sat Dec  5 12:06:45 2015

@author: chris
"""

import pdb
import numpy as np

from numpy import array,float64, pi,exp, transpose,sign
import scipy.optimize
import pandas
#from ramanTools.SPETools import File
from copy import deepcopy,copy
import scipy.optimize
import matplotlib.pyplot as plt
from collections import namedtuple
import inspect



def import_IR_spectrum(orca_output_file, normalize = False,color='k',labelpeaks = True, sticks = True, freq_factor = 1, width=0):

    fileopen = open(orca_output_file,'r')
    f=fileopen.readlines()
    fileopen.close()
    mode_numbers = []
    for l in f:
        try:
            if 'IR SPECTRUM' in l:
                
                start =  f.index(l)+5
        except:
            print('error')
    table=list()
    for z in f[start:]:
        i = z.split(' ')
        while '' in i: 
            i.remove('')
        
       
        if i[0][0].isdigit():# if i[0][0].isdigit():
            print(i[0])
            
            
            mode_numbers.append(int(i[0].replace(':','')))
            table.append([float(i[1]),float(i[2])])
#            if labelpeaks:   
#                plt.gca().annotate(i[0], (float(i[1]), float(i[2])+0.2), color=color,fontsize = 8,horizontalalignment='center') 
        elif not i[0][0].isdigit():
            
            break
    table = transpose(table)
    if normalize == True:
            table[1]/=max(table[1])
    print(table)
    table[0]*=freq_factor
    
    if labelpeaks:
        
        for i in np.arange(table.shape[1]):
            
            plt.gca().annotate(str(i+min(mode_numbers)), (float(table[0][i]), float(table[1][i])*1.1), color=color,fontsize = 8,horizontalalignment='center')
    if sticks==True:
        ret = plt.vlines(table[0],0,table[1],linewidth = 2,color=color)
    if width>0:
        x=np.arange(0,4000,0.1)
        y=np.zeros(x.shape)
        for i in range(table.shape[1]):
            
            x0 = table[0][i]
            A = table[1][i]
            y+= A*exp(-(x-x0)**2/width**2)
        ret= plt.plot(x,y,color=color)
    return ret
            
        
    
#def import_TDDFT_spectrum(orca_output_file,axistoploton=None,broaden=False, normalize = False,color='k',labelpeaks = True):
#
#    fileopen = open(orca_output_file,'rb')
#    f=fileopen.readlines()
#    fileopen.close()
#    for l in f:
#        if 'ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS' in l:
#            start =  f.index(l)+5
#    table=list()
#    for z in f[start:]:
#        i = z.split(' ')
#        while '' in i: 
#            i.remove('')
#       
#        if i[0][0].isdigit():
#           
#            table.append([float(i[2]),float(i[4])])
#            if labelpeaks:   
#                plt.gca().annotate(i[0], (float(i[2]), float(i[4])+0.2), color=color,fontsize = 8,horizontalalignment='center') 
#        elif not i[0][0].isdigit():
#            
#            break
#    table = transpose(table)
#
#    if normalize == True:
#            table[1]/=max(table[1])
#    if broaden:
#        gamma=3
#        freqs = table[0]
#        amps = table[1]
#        x = arange(freqs[-1]-10,freqs[0]+100,0.1)
#        
#        y = zeros(x.shape)
#        
#        for i in range(amps.size):
#        
#            y+=    (1/pi)*amps[i]*(0.5*gamma)/((x-freqs[i])**2+(0.5*gamma)**2)
#        
#        if axistoploton==None:
#            axistoploton=gca()
#        return axistoploton.plot(x,y)
#        
#    return plt.vlines(table[0],0,table[1],linewidth = 2,color=color)
#    
#def separateNWgeo(filename):
#
#    fileopen = open(filename,'rb')
#    f=fileopen.readlines()
#    fileopen.close()
#    namefile = 0
#    indexlist = []
#    
#    
#    
#    for l in range(len(f)):
#        if 'Geometry "geometry" -> "geometry"' in f[l]:
#            indexlist.append(l) 
#            print (indexlist)
#    with open('/home/chris/Desktop/geos/animation.v000.xyz','wb') as animfile:
#        animfile.write('')
#        animfile.close()
#    tablelist = []
#    lengthoflist = 0
#    for start in indexlist:
#            
#            start =  start+7#f.index(l)+7
#            table=str()
#            
#            for z in f[start:]:
#                
#                i = z.split(' ')
#                while '' in i: 
#                    i.remove('')
#               
#                if i[0][0].isdigit():
#                    lengthoflist = i[0]
#                    table+=(str(i[1])+ ' ' + str(i[3]) + ' ' + str(i[4])+' '+ str(i[5][:-1]) +'\n')# ' 0 0 0\n')
#                    
#                elif not i[0][0].isdigit():
#                    xyzfiletext = str(lengthoflist)+'\nnocomment\n'+table
#                    tablelist.append(table)
#                    with open('/home/chris/Desktop/geos/xyzfile'+str(namefile)+'.xyz', 'wb') as writefile:
#                        writefile.write(xyzfiletext)
#                        writefile.close()
#                        namefile+=1
#                        writefile.close
#                    break
#    with open('/home/chris/Desktop/geos/animation.v058.xyz','wb') as animfile:
#        
#        for i in range(len(tablelist)):
#           
#            table = tablelist[i]
#            
#            animfile.write(str(lengthoflist)+'\n')
#            
#            if i==0:
#                animfile.write('* (null), Energy   -1000.0000000 ' +str(i)+'\n' )
#                
#            else:    
#                animfile.write('* Frequency(58) 1172.5400 ' +str(i)+'\n' )
#            
#            animfile.write(table+'\n')
#            
##        for i in range(len(tablelist)-1,-1,-1):
##            print i
##            table = tablelist[i]
##            
##            animfile.write(str(lengthoflist)+'\n')
##            animfile.write('* Frequency(58) 1172.5400 ' +str(i)+'\n' )
##            animfile.write(table+'\n')
##            
#            
#            
#        animfile.close()
#                            
#                   
# 
#                
#    indexlist = []            
#    for l in range(len(f)):
#        if 'Total DFT energy =' in f[l]:
#            indexlist.append(l)
#    energies=array([])
#    for start in indexlist:
#
#            z= f[start]
#            #print z
#                
#            i = z.split(' ')
#            while '' in i: 
#                i.remove('')
#                
#            energies = np.append(energies, float(i[4]))
#    plot(energies)
#    ylabel('energy')
#    xlabel('step')
#          
#    
#    return 0
#    
#def NWchemTDDFT(filename,energyunits = 'ev',triplets=False,lines=True,spectrum=True,axistoploton=None,peakstoplot=None):
#    
#    if axistoploton==None:
#        axistoploton=plt.gca()
#    if triplets:
#        pass
#    
#    
#    fileopen = open(filename,'rb')
#    f=fileopen.readlines()
#    fileopen.close()
#    namefile = 0
#    indexlist = []
#    
#    start = 0
#    for i in f:
#        if 'NWChem TDDFT Module' in i:
#            start = f.index(i)
#    if start == 0: return array([-1])
#            
#    singlet_energies = array([])
#    triplet_energies = array([])
#    singlet_osc_strengths = array([])
#    triplet_osc_strengths = array([])
#   
#    for index in range(start, len(f)):
#        
#        i =f[index]
#        if 'Excited state energy' in i:
#            print( 'found maker', i, index)
#            start = index
#            break
#        elif 'Root' in i:
#            z = i.split(' ')
#            energy = float(z[-3])
#            singlet_energies = np.append(singlet_energies, [energy])
#            oscillator = float(f[f.index(i)+5].split(' ')[-1])
#            singlet_osc_strengths = np.append(singlet_osc_strengths,oscillator)
#          
#    for index in range(start,len(f)):
#        if 'Davidson' in f[index]:
#            start = index
#  
#    if triplets:
#        for index in range(start+1, len(f)):
#            
#            i =f[index]
#            
#            if 'Excited state energy' in i:
#              
#                break
#            elif 'Root' in i:
#                
#                z = i.split(' ')
#                
#                energy = float(z[-3])
#                
#                triplet_energies = np.append(triplet_energies, [energy])
#                oscillator = 1
#                triplet_osc_strengths = append(triplet_osc_strengths,oscillator)
#    
#
#    if True:
#        if peakstoplot==None:
#            peakstoplot=len(singlet_energies)
#        gamma=0.005
#       
#        freqs = singlet_energies[:peakstoplot]
#        amps = singlet_osc_strengths[:peakstoplot]
#        x = np.arange(1,5,0.0001)
#        
#        y = np.zeros(x.shape)
#        
#        for i in range(amps.size):
#            print (amps[i], freqs[i])
#            y+=    (1/pi)*amps[i]*(0.5*gamma)/((x-freqs[i])**2+(0.5*gamma)**2)
#            if lines:
#                plt.vlines(freqs[i], 0, amps[i]*10,colors = 'k')
#        if triplets:
#            
#            for i in range(len(triplet_energies)):plt.vlines(triplet_energies[i], 0, 1,colors = 'r')
#
#        
#        if energyunits=='nm':
#            if spectrum:
#                axistoploton.plot(1240/x,y)
#            return 1240/freqs,amps
#        elif energyunits=='ev' or  energyunits=='eV':
#            if spectrum:
#                axistoploton.plot(x,y)
#            return freqs,amps
#    
#    return freqs
#    
#def NWChemDOS(filename,axis = plt,_plot=False):
#    
#              
#    fileopen = open(filename,'rb')
#    f=fileopen.readlines()
#    fileopen.close()
#    namefile = 0
#    indexlist = []
#    
#    
#    for i in range(len(f)-1, 0,-1):
#       
#        if 'DFT Final Molecular Orbital Analysis' in f[i]:
#            start = i
#            break
#            
#    singlet_energies = np.array([])
#    occs = array([])
#    singlet_osc_strengths = np.array([])
#    triplet_osc_strengths = np.array([])
#    
#
#    for index in range(start, len(f)):
#        i =f[index]
#        if 'center of mass' in i:
#            
#            break
#        elif 'Vector' in i:
#            z = i.split('E=')
#            occupation = float((z[0].split('Occ=')[1][0:4]))
#            
#            x = float(z[-1].replace('D','E'))
#           
#            if x == '':
#                continue
#            
#            orbital_energy = x#float(z[-1].split(' ')[0].replace('D','E'))
#
#            singlet_energies = np.append(singlet_energies, [orbital_energy])
#            occs = np.append(occs, occupation)
#    energy = 0
#    for l in range(0,len(f)):
#        if 'Total DFT energy =' in f[l]:
#            z= f[l]
#                
#            i = z.split(' ')
#            while '' in i: 
#                i.remove('')
#            
#            energy=float(i[4])
#            print( 'energy', energy)
#                
#    print ('the energy of the groundstate is', energy)
#            
#    
#    energies = singlet_energies*27.21138386  
#    
#    occupied_energies = energies[occs==2.0]  
#    unoccupied_energies = energies[occs==0.0]
#    homo = np.max(occupied_energies)
#    lumo = np.min(unoccupied_energies)
#    gap = lumo-homo
#    hist_bins = np.arange(min(energies),np.max(energies),0.2)
#    print ("HOMO: number",len(occupied_energies),'energy:', np.max(occupied_energies))
#    print ("LUMO: number",len(occupied_energies)+1,'energy:', np.min(unoccupied_energies))
#    print ("Gap:", gap)
#    if _plot:
#        if axis == None:
#            axis = gca()
#        realDOS = histogram(occupied_energies,range = (min(occupied_energies)-2,max(occupied_energies)+2),bins=round((max(occupied_energies)+2-min(occupied_energies)-2)/0.2))
#        virtualDOS = histogram(unoccupied_energies,range = (min(unoccupied_energies)-2,max(unoccupied_energies)+2),bins=round((2+max(unoccupied_energies)-min(unoccupied_energies)-2)/0.2))
#        print (realDOS)
#        axis.fill_between(realDOS[1][1:],realDOS[0],color= 'b')
#        axis.fill_between(virtualDOS[1][1:],virtualDOS[0],color = 'r')
#        axis.legend(['occupied', 'unoccupied'])
#    
#    data = (energy,gap,homo,lumo)#np.array([(energy,gap,homo,lumo)], dtype=[ ('energy', float),('gap', float), ('homo', float),('lumo',float)])
#    return data
#        
#        
#    
#def makeaspectrum():
#    a = loadtxt('/home/chris/Desktop/TDDFT_.csv', unpack=True,delimiter = ',')
#    
#  
#    x = arange(0.01,10,0.01)
#    y=zeros(x.shape)
#    gamma=0.1
#
#    for i in range(a.shape[1]):
#        
#            y+=    (1/pi)*a[1,i]*(0.5*gamma)/((x-a[0,i])**2+(0.5*gamma)**2)
#        
#    
#    return (x,y)
#    
#def blacktored(ax):
#    dist = len(ax.lines)
#    c=0
#    for i in ax.lines:
#        c+=1.0/dist
#        
#        i.set_color((c,0,0))
#    return 0
#    
#def saveaxisexcel(filename,axis,header = ''):
#
#    lengths = list((len(x.get_xdata()) for x in axis.lines))
#    
#    maxlen = max(lengths)
#    
#    tosave= np.ndarray((0,maxlen))
#    
#    for l in range(len(lengths)):
#        padnumber = maxlen-lengths[l]
#        
#        addtosave = np.append(array(axis.lines[l].get_data()),np.zeros((2,padnumber)),axis=1)
#        
#        tosave=  np.append(tosave,addtosave ,axis = 0)
#    np.savetxt(filename,np.transpose(tosave), delimiter=',',header = header )
#    return tosave
#        
#        
#
#    