# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 10:38:30 2015

@author: chris
"""

#
#
#global root 
#    
#root = Tk()
#FittingWindow(root)
#
#root.mainloop()
#  
#  
#import sys


import ramanTools.RamanTools, ramanTools.RamanSpectrum
from tkinter import *
import os
import matplotlib.pyplot
matplotlib.pyplot.ioff()
os.chdir('C:/Users/cthompson')
import sys
sys.path.append('C:/Users/cthompson/Dropbox/PyScripts')

root = Tk()
root.withdraw()
ramanTools.RamanTools.FittingWindow(root)
root.mainloop()