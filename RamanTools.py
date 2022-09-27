# -*- coding: utf-8 -*-
"""
Created on Thu Oct 16 11:37:39 2014

@author: chris
<<<<<<< HEAD

Raman Tools contains classes for GUI of spectral processing
=======
RamanTools contains the classes and objects for GUI of data processing
>>>>>>> newramanToolbranchFeb19
"""



from numpy import *
import scipy.optimize
import pandas
from tkinter import *
import tkinter as Tkinter
import tkinter.filedialog as tkFileDialog
import tkinter.filedialog, tkinter.simpledialog, tkinter.messagebox
from matplotlib.pyplot import *
from matplotlib.axes import *
import ramanTools.SFG_Notebook 

from matplotlib.figure import Figure 
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg,NavigationToolbar2Tk
import os

from copy import copy

import pdb
from matplotlib.widgets import SpanSelector
from collections import deque
from matplotlib.colors import ColorConverter

from ramanTools.RamanSpectrum import *
from matplotlib.colors import to_hex

#colordict= ColorConverter()

#list_of_display_windows = []

def RGBtohex(rgbtuple): return '#%02x%02x%02x'%  (int(rgbtuple[0]*255),int(rgbtuple[1]*255),int(rgbtuple[2]*255))

def modify_plot(ax):
    root = Tk()
    root.withdraw()
    DisplayWindow(root,ax=ax)
    root.mainloop()
    
class DisplayWindow(tkinter.Toplevel):
        
        color_list = deque(['g','r','c','m','y','k','r','b'])
        current_checker = None
        list_of_display_windows = list()

        def __init__(self, rootref,ax = None):
            
                tkinter.Toplevel.__init__(self)

                self.title(string = "Viewer")
               
                self.checker_list = []
                self.checker_frame=Frame(master=self, width =500)

                self.menubar = Menu(master=self)

                self.filemenu = Menu(self.menubar, tearoff=0)
                self.filemenu.add_command(label="New window", command = DisplayWindow)
                self.filemenu.add_command(label="Open", command = self.DisplaySpectrum)
                self.filemenu.add_command(label="Save",command = self.SaveSpectrum)
                self.filemenu.add_command(label="SaveFig",command = self.SaveFigure)
                self.filemenu.add_command(label="Exit", command = self.quitproc)
                self.menubar.add_cascade(label="File", menu=self.filemenu)
                # create more pulldown menus
                self.editmenu = Menu(self.menubar, tearoff=0)
                self.editmenu.add_command(label="Smooth", command = self.SmoothSpectrum)
                self.editmenu.add_command(label="Baseline", command = self.Baseline)
                self.editmenu.add_command(label="Copy")
                self.editmenu.add_command(label="Paste")
                self.menubar.add_cascade(label="Edit", menu=self.editmenu)

                self.config(menu= self.menubar)
                self.fig = figure()
                self.ax=self.fig.add_subplot(111)
                
                
                
                self.canvas = FigureCanvasTkAgg(self.fig,master = self)
                self.canvas.draw()
                self.toolbar = NavigationToolbar2Tk(self.canvas,self)
                self.canvas._tkcanvas.pack(side=TOP,expand =True)
                
                
                
                
                self.checker_frame.pack(side=BOTTOM)
                self.toolbar.pack()
                
                
                self.toolbar.update()

                self.protocol("WM_DELETE_WINDOW", self.quitproc)
                #-------------------------------------------------------------
                
                self.rootref = rootref

#                
                self.list_of_display_windows.append(self)
           
                return None
        def popup(self,event):
            self.left_click_menu.post(event.x_root, event.y_root)
        def update_legend(self):
            
            if self.legend_var.get() == 0:
                self.Plot.a.get_legend().set_visible(False)
            else:
                l = list()
                for entry in self.checker_list:
                    if True:#entry.get_visible():
                        l.append(entry.commentbox.get())
                        self.Plot.a.legend(l)
                   
          
            self.Plot.canvas.draw()
            return 0
        def showinfo(self):
            tkMessageBox.showinfo('Spectrum info',self.current_checker.spectrum.info)
            return 0
            
        def getcolor(self): 
            self.color_list.rotate(1)
            
            return self.color_list[0]
                
        def DisplaySpectrum(self, str_name=None, spectrum=None):
            import re
            options = {}
            options['defaultextension'] = '.csv'
            options['filetypes'] = [('all files', '.*')]
            options['initialdir'] = 'C:\sfg\data'
            options['multiple'] = True
            options['title'] = 'Open Spectrum...'
            options['parent'] = self.master

            if type(spectrum) == RamanSpectrum:
                a=[spectrum.index, spectrum.values]
            elif type(str_name)== str:# != None: 
                a=np.genfromtxt(str_name,delimiter = ' ',unpack=True)
            else:
                str_name_list = tkFileDialog.askopenfilenames(**options)
                if str_name_list == None:
                    return
                elif str_name_list == '':
                    return
                elif type(str_name_list) == str:
                    a=np.genfromtxt(str_name,delimiter = '\t',unpack=True)

                elif type(str_name_list) == tuple or type(str_name_list) == list:
                    for str_name in str_name_list:
                        try:
                            a=RamanSpectrum(str_name, filetype='FTIR ASC')#np.genfromtxt(str_name,delimiter = '\t',unpack=True)
                            
                        except:
                            raise
                            filetype = filetypewindow()
                            
                            ## if filetype in [,]:
                            return 
                        r=Channel(array(a.index),a.values,None)
                        self.ax.add_line(r)
                        r.set_color(self.getcolor())
                        self.checker_list.append(checker(self.checker_frame,r, spectrum_name = os.path.basename(str_name)))
                    self.ax.autoscale_view()
                    self.fig.canvas.draw()
                    return
                        
                else:
                    print('I dont know what to do with this', str_name_list)
 
            r=Channel(a[0],a[1],None)
            self.ax.add_line(r)
            r.set_color(self.getcolor())
            self.checker_list.append(checker(self.checker_frame,r, spectrum_name = os.path.basename(str_name)))
            self.ax.autoscale_view()
            self.fig.canvas.draw()
            return

            
        def SaveSpectrum(self):
                
            file_opt = options = {}
            options['defaultextension'] = '.txt'
            options['filetypes'] = [('all files', '.*')]
            
            options['title'] = 'Open Spectrum...'
            options['initialfile'] = os.path.basename(self.current_checker.name)
            options['parent'] = self.master

            if self.current_checker  == None:
                    return 0
            str_filename = tkFileDialog.asksaveasfilename(**options)
           
            if str_filename == '':
                    return 0
            else:
                    data = transpose([self.current_checker.get_xdata(),self.current_checker.get_ydata()])
                    
                    savetxt(str_filename,data,delimiter = ',',comments = '#'+str(self.current_checker.spectrum.info))#SaveSpectrum(self.current_checker.channel.spec_array,str_filename)
            os.chdir(os.path.dirname(str_filename))
            return 0
        def ViewNotebook(self):
                pass
                return 0 
        def zeroall(self):
                for c in ax.lines:
                    c.set_ydata(c.get_ydata()-min(c.get_ydata()))
                self.fig.canvas.draw()
                return 0
                
                        
        def SmoothSpectrum(self):

                if self.current_checker == None:
                    return
                self.current_checker.spectrum.smooth()

                self.current_checker.channel.set_ydata((self.current_checker.spectrum))
                self.fig.canvas.draw()
                return
            
        def Baseline(self):
            self.baseline = np.ndarray((0,2))
            self.connection_id = self.fig.canvas.mpl_connect('button_press_event', lambda event: self.onclick(event))
            self.connection_id2 = self.fig.canvas.mpl_connect('key_press_event', lambda event: self.retclick(event))
            self.baseline_dots = self.ax.plot([],[],'s')[-1]
        def onclick(self, event):
            
            if event.dblclick:
                if self.baseline.shape[0]==0:
                    ### if this is the first point, add it to the baseline
                    self.baseline = np.append(self.baseline, np.array([[event.xdata,self.current_checker.spectrum[event.xdata]]]), axis=0)
                    self.baseline=self.baseline[self.baseline[:,0].argsort()]#.sort(axis = 0)
                    self.baseline_dots.set_data(self.baseline[:,0], self.baseline[:,1])
                    self.fig.canvas.draw()
                    return
                
                z = array([[event.xdata,event.ydata]]*self.baseline.shape[0])
                b = (self.baseline-z)/z
                if np.any(np.sum(b**2, axis = 1)<0.05):
                    ## if double click close to already existing point, delete that point.
                    to_remove = np.argmin(np.sum(b**2, axis = 1))
                    self.baseline=np.delete(self.baseline, to_remove, axis = 0)
                    
                else:
                    ## if the indicated point is new, add the point to baseline
                    self.baseline = np.append(self.baseline, np.array([[event.xdata,self.current_checker.spectrum[event.xdata]]]), axis=0)
                    self.baseline=self.baseline[self.baseline[:,0].argsort()]
                    
                self.baseline_dots.set_data(self.baseline[:,0], self.baseline[:,1])
                self.fig.canvas.draw()
            return
                
        def retclick(self,event):
            
            if event.key == 'a': 
                if self.baseline.shape[0]<2:
                    pass
                else:
                    self.current_checker.spectrum.autobaseline(self.baseline[:,0], specialoption='points', order = self.baseline.shape[0]-2)
                    self.current_checker.channel.set_ydata((self.current_checker.spectrum))
                self.fig.canvas.mpl_disconnect(self.connection_id2)
                self.fig.canvas.mpl_disconnect(self.connection_id)
                self.ax.lines.remove(self.baseline_dots)
                self.fig.canvas.draw()
                

        def open_next_spectrum_in_folder(self):
                pass
                
                return 
        def open_next_spectrum(self,event):
                
                pass
                return 
                
                
        def normalizeall(self):
            for check in self.checker_list:
                data = check.get_ydata()
                data[:]-=min(data)
                data/=max(data)
                check.set_ydata(data)
            self.Plot.a.relim()
            self.Plot.a.set_ylim(-0.5,1.5)
            self.Plot.a.autoscale_view(tight = False)
            self.Plot.canvas.draw()
            return 0
                
        def start_calc_noise(self):
              
                self.span = SpanSelector(self.Plot.a, self.calc_noise, 'horizontal')
                self.span.connect_event('pick_event',self.calc_noise)
                gcf().canvas.mpl_connect('button_press_event',self.disconnect)
                return 0
        def start_calc_area(self):
               
                self.span = SpanSelector(self.Plot.a, self.calc_area, 'horizontal')
                self.span.connect_event('pick_event',self.calc_area)
                gcf().canvas.mpl_connect('button_press_event',self.disconnect)
                
                return 0
                
        def calc_noise(self,start,end):
                try:
                    print ("STD =", calc_noise(pandas.Series(self.current_checker.channel.get_ydata(),self.current_checker.channel.get_xdata()),(start,end)))
                 
                except:
                    print ('error')
                return 0
        def calc_area(self, start,end):
                try:
                    print ("Area =", calc_area(pandas.Series(self.current_checker.channel.get_ydata(),self.current_checker.channel.get_xdata()),(start,end)) )
                except:
                    print( 'error')
                return 0

        def removenoise(self):
            pass
            return 

        def RemoveSpectrum(self,channel):
            self.checker_list.remove(channel.checker)
            channel.checker.destroy()
            self.ax.lines.remove(channel)
            self.ax.autoscale_view()
            self.fig.canvas.draw()
            return 
                

        def FFT(self):
                return 0
                
        def ShowChannel(self,channel):
                channel.set_visible(True)
                self.fig.canvas.draw()
                return 0

        def HideChannel(self,channel):
                channel.set_visible(False)
                self.fig.canvas.draw()
                return 0
  
                        
        def SaveFigure(self):
                self.fig.savefig('figure.png')
                return 0

        
        def quitproc(self):
                self.destroy()
                self.rootref.destroy()
                self.rootref.quit()

                return

class filetypewindow(tkinter.Toplevel):
    
    def __init__(self, variable_to_change):
        tkinter.Toplevel.__init__(self)
        self.listbox = Listbox(self)
        
        
        self.listbox.insert(END, "a list entry")

        for item in ["Raman GVD", "FTIR GVD",]:
            self.listbox.insert(END, item)
            
        self.selectbutton = tk.Button(self, text = 'OK', command = self.ok)
        
        self.listbox.grid(row=0,column = 0)
        self.selectbutton.grid(row=3,column=3)
            
        return 
    
    def ok(self):
        
        self.destroy()

    

class checker(tkinter.Frame):
         
        def __init__(self,master,channel, spectrum_name=""):
                self.master=master
                tkinter.Frame.__init__(self, master  = master)
                self.window = self.master.master
                self.xoffset = 0
                
                self.pack()
                                   
                self.channel = channel
                self.channel.checker = self
                self.spectrum = RamanSpectrum(pandas.Series(channel.get_ydata(),channel.get_xdata()))   
                self.spectrum.set_name(spectrum_name)
                self.visible_var = IntVar()
                self.visible_var.set(1)
                
                self.spectrumname_var = StringVar()
                self.spectrumname_var.set(spectrum_name)
                self.box = Checkbutton(self,
                                       variable=self.visible_var,
                                       command=self.cb)
                self.box.grid(row=0,column = 5)

                self.DeleteButton = tkinter.Button(self,
                                        text = "Delete",
                                        command = lambda: self.window.RemoveSpectrum(self.channel),
                                        width=5,
                                        height=1)
                self.DeleteButton.grid(row=0,column = 7)

                self.NameLabel =  Label(self,width=30, height = 1, text = self.spectrum.name)
                self.NameLabel.bind("<Button-1>",self.set_current)
                self.NameLabel.grid(row=0,column = 1)
                
                self.w = Canvas(master = self,width = 20,height=20)
                self.w.grid(row = 0, column=10)#self.w.pack(side = LEFT,padx=10)
                
               
                fillcolor = to_hex(self.channel.get_color())#RGBtohex(colordict.colors[])

                self.w.create_rectangle(0,0,100,100,  fill=fillcolor)
                self.w.bind("<Button-3>", self.popup)

                
                self.var_channel_type = StringVar()
                self.var_channel_type.set('SFG')
            
                
                return
        
                 
        def set_current(self,event):
                
                if self.window.current_checker == None:
                        pass
                else:
                        self.window.current_checker.NameLabel.config(relief = FLAT)
                self.window.current_checker = self
                self.NameLabel.config(relief = SUNKEN)
                return 0 
                
                
        def cb(self):
                
                if self.visible_var.get() == 0:
                        self.window.HideChannel(self.channel)
                        
                else:
                        self.window.ShowChannel(self.channel)
                
                return 0
       
        def SetChannelType(self,channel_type):
                
                return
        def set_xoffset(self,offset):
            self.channel.set_xdata(self.channel.get_xdata()-offset)
            self.xoffset = offset
            
            return

        def popup(self,event):
            return
        
        
        
class Channel(matplotlib.lines.Line2D):
    
    def __init__(self, xdata, ydata,checker):
        matplotlib.lines.Line2D.__init__(self, xdata,ydata)
        self.checker = checker
        
        return 


       
class FittingWindow(object):
    j=np.complex(0,1)

    def __init__(self, master,ramanspectrum =None):
        global guessplot,dataplot
        
        self.master  = master
        self.textframe= Frame(master = self.master)
        self.plotframe = Frame(master = self.master)
        self.buttonframe = Frame(master = self.master)
        self.frame = Frame(master = self.master)
        
        self.scroll = Tkinter.Scrollbar(self.textframe)
        self.scroll.grid(row=0,column=1)
        self.t = Tkinter.Text(self.textframe,yscrollcommand=self.scroll.set,width=30)
        self.scroll.config(command=self.t.yview)
        
        ##############################
        self.menubar = Menu(self.master)
        
        self.filemenu = Menu(self.menubar, tearoff=0)
        self.filemenu.add_command(label="New window", command = lambda: FittingWindow(Toplevel()))
        self.filemenu.add_command(label="Open", command = self.open_data)
        self.filemenu.add_command(label="Save",command = None)
        self.filemenu.add_command(label="SaveFig",command = None)
        self.filemenu.add_command(label="ViewNotebook",command = None)
        self.filemenu.add_command(label="Exit", command = self.quitproc)
        self.menubar.add_cascade(label="File", menu=self.filemenu)
        
        self.master.config(menu= self.menubar)
        ##############################
        
        self.fit_button = Tkinter.Button(master = self.buttonframe, command = self.fit_and_draw, width  = 5, text = 'FitNow')
        self.fit_button.grid(row = 3, column = 0)
        self.open_button =  Tkinter.Button(master = self.buttonframe, command = self.open_data, width  = 5, text = 'open')
        self.open_button.grid(row = 4, column = 0)

        
        self.smoothbutton = Tkinter.Button(master = self.buttonframe, command = self.smoothdata, width  = 5, text = 'Smooth')
        self.smoothbutton.grid(row = 5, column = 0)

        self.function_list=[
                              'OneGaussian',
                              'TwoGaussian',
                              'ThreeGaussian',
                              'FourGaussian',
                              'FiveGaussian',
                              '5coherent']
        
        self.var_func = StringVar()
        self.var_func.set(self.function_list[0])

        self.MotorMenu = OptionMenu(self.buttonframe,self.var_func, *self.function_list, command = self.init_function)

        self.normalization_constant = 1
        self.scale_list = []
        self.var_list = []
        
        self.fig = Figure(figsize=(5,3))
        self.ax1  = self.fig.add_subplot(111)
        dataplot = self.ax1.plot(arange(2800,3605,5),zeros((161,)),animated = True)[0]#(ax = self.ax1).lines[-1]
        guessplot = self.ax1.plot(arange(2800,3605,5),zeros((161,)),animated = True)[0]
        
        
        self.canvas = FigureCanvasTkAgg(self.fig,master = self.plotframe)

        
        self.toolbar = NavigationToolbar2TkAgg(self.canvas, self.plotframe)
        self.toolbar.update()

        self.reference = numpy.ndarray((161,))
        self.reference[:] = 1
        self.reference_name  = ''
        self.guess = list()

        self.name = ['']
        if ramanspectrum is None:
            self.a = pandas.Series(zeros(161),arange(2800,3605,5))
        else:
            self.a = copy(ramanspectrum)
            self.t.insert(END,'data multiplied by'+str(1/max(self.a)))
            self.a[:]/=max(self.a.values)
         
        self.startfreq_text = Entry(self.buttonframe,width = 5)
        self.startfreq = min(array(self.a.index))
        self.startfreq_text.insert(END,str(self.startfreq))
        self.startfreq_text.bind("<Return>", self.update_limits)

        self.endfreq_text = Entry(self.buttonframe,width = 5)
        self.endfreq = max(array(self.a.index))
        self.endfreq_text.insert(END,str(self.endfreq))
        self.endfreq_text.bind("<Return>", self.update_limits)
        
        
        if self.init_function('OneLorentzian') == -1:
            self.t.insert(END, 'error initializing fit function')

        self.textframe.grid(row=0,column=0,columnspan=1)
        self.plotframe.grid(row = 0, column = 3, columnspan = 11,sticky = 'e')
        self.frame.grid(row = 1,column = 1, columnspan = 11, rowspan = 3,sticky = 'se')
        self.buttonframe.grid( row =1,column=0,columnspan=1)

        
        
        
        
        ######Items in plot frame
        self.canvas._tkcanvas.pack(side=TOP, fill = BOTH,expand =True)    
        self.toolbar.pack(side= BOTTOM, fill = BOTH,expand =True)#,expand = True)#fill=BOTH)#, expand=1)
        
        ### Items in button frame
        self.endfreq_text.grid(row = 1 , column =1,sticky = W)
        self.startfreq_text.grid(row = 0 , column =1,sticky = W)
        self.MotorMenu.grid(row = 3,column =1,sticky = W)
        
        self.ax1.cla()
        self.canvas.show()
        self.background = self.canvas.copy_from_bbox(self.ax1.bbox)
        
        #### items in text frame
        self.t.grid(row=0,column=0)

        self.ax1.set_xlim(min(self.a.index),max(self.a.index))
        self.ax1.set_ylim(min(self.a.values),1)

        self.raw = self.a.copy   
        self.startidx = 0
        self.endidx = -1
        try:
            self.update_limits(0)
        except:
            pass
        return None

    def quitproc(self):
        
            self.destroy()
            return 0
        
        
    def init_function(self,functiontype):
        global function
        import inspect
        j=np.complex(0,1)
        
        if functiontype == None:
            functiontype = self.var_func.get()
        
        self.functiontype = functiontype
        
      
                    
        if functiontype == 'OneLorentzian':
            def function(x,A1,w1,G1,m,b): return m*x/1000+b + A1**2/((x-w1)**2+G1**2)
        elif functiontype == 'TwoLorentzian':
            def function(x,A1,A2,w1,w2,G1,G2,m,b): return m*x/1000+b + A1**2/((x-w1)**2+G1**2) +A2**2/((x-w2)**2+G2**2)
    
        elif functiontype == 'ThreeLorentzian':
            def function(x,A1,A2,A3,w1,w2,w3,G1,G2,G3,m,b): return m*x/1000+b + A1**2/((x-w1)**2+G1**2) +A2**2/((x-w2)**2+G2**2) +A3**2/((x-w3)**2+G3**2)
        elif functiontype == 'FourLorentzian':
            def function(x,A1,A2,A3,A4,w1,w2,w3,w4,G1,G2,G3,G4,m,b): return m*x/1000+b + A1**2/((x-w1)**2+G1**2) +A2**2/((x-w2)**2+G2**2) +A3**2/((x-w3)**2+G3**2) +A4**2/((x-w4)**2+G4**2)
             
        elif functiontype == 'FiveLorentzian':
            def function(x,A1,A2,A3,A4,A5,w1,w2,w3,w4,w5,G1,G2,G3,G4,G5,m,b): return m*x/1000+b + A1**2/((x-w1)**2+G1**2) +A2**2/((x-w2)**2+G2**2) +A3**2/((x-w3)**2+G3**2) +A4**2/((x-w4)**2+G4**2)+A5**2/((x-w5)**2+G5**2)
        elif functiontype == 'OneGaussian':
            def function(x,A1,w1,G1,m,b): return m*x/1000+b + A1*exp(-(x-w1)**2/G1)
        elif functiontype == 'TwoGaussian':
            def function(x,A1,A2,w1,w2,G1,G2,m,b): return A1*exp(-(x-w1)**2/G1) +A2*exp(-(x-w2)**2/G2) +m*x/1000+b 
    
        elif functiontype == 'ThreeGaussian':
            def function(x,A1,A2,A3,w1,w2,w3,G1,G2,G3,m,b): return m*x/1000+b + A1*exp(-(x-w1)**2/G1) +A2*exp(-(x-w2)**2/G2) +A3*exp(-(x-w3)**2/G3) 
        elif functiontype == 'FourGaussian':
            def function(x,A1,A2,A3,A4,w1,w2,w3,w4,G1,G2,G3,G4,m,b): return m*x/1000+b +A1*exp(-(x-w1)**2/G1) +A2*exp(-(x-w2)**2/G2)  +A3*exp(-(x-w3)**2/G3) +A4*exp(-(x-w4)**2/G4) 
             
        elif functiontype == 'FiveGaussian':
            def function(x,A1,A2,A3,A4,A5,w1,w2,w3,w4,w5,G1,G2,G3,G4,G5,m,b): return m*x/1000+b + (A1/G1*numpy.sqrt(2*pi))*exp(-(x-w1)**2/(2*G1**2)) +(A2/G2*numpy.sqrt(2*pi))*exp(-(x-w2)**2/(2*G2**2))  +(A3/G3*numpy.sqrt(2*pi))*exp(-(x-w3)**2/(2*G3**2))  +(A4/G4*numpy.sqrt(2*pi))*exp(-(x-w4)**2/(2*G4**2)) +(A5/G5*numpy.sqrt(2*pi))*exp(-(x-w5)**2/(2*G5**2)) 
        
        elif functiontype == '5coherent':
            def function(x,A1,A2,A3,A4,A5,w1,w2,w3,w4,w5,G1,G2,G3,G4,G5,m,b): 
                return np.real(m*x/1000+b + A1/((x-w1)+j*G1)+
                A2/((x-w2)+j*G2) +
                A3/((x-w3)+j*G3) +
                A4/((x-w4)+j*G4)+
                A5/((x-w5)+j*G5))**2
        else:
            tkMessageBox.showerror('Not regconized function')
            def function(x,m,b): return m*x/1000+b

#==============================================================================
        self.w_name = list() 
        self.w_name = inspect.getargspec(function).args[1:]
        
        self.guess= []
       
        for w in self.w_name:
            if 'A' in w:
                self.guess.append(0.1)
            elif 'w' in w:
                self.guess.append(1000)
            elif 'G' in w:
                self.guess.append(10)
            else:
                self.guess.append(0)
#================================================================================
        
        for i in range(len(self.scale_list),len(self.w_name)):
                self.scale_list.append(Scale(master  = self.frame,command = self.guessdraw))
                
        for i in range(len(self.w_name),len(self.scale_list)):
                self.scale_list[i].grid_forget()
           
        for i in range(len(self.w_name)):
                self.scale_list[i].config(label = self.w_name[i])
                self.scale_list[i].grid(row = int(i/8), column = i%8)
                
                if self.w_name[i][0] == 'A':
                     self.scale_list[i].config(from_ =0, to =  5,resolution = 0.1)
                     
                     self.scale_list[i].set(1)
                elif self.w_name[i][0] == 'w':
                    self.scale_list[i].config(from_ = self.startfreq-20, to =  self.endfreq+20,resolution = 1)
                    self.scale_list[i].set(2800)
                elif self.w_name[i][0] == 'G':
                    if 'coherent' in functiontype:
                        self.scale_list[i].config(from_ = -20, to =  20)
                        self.scale_list[i].set(7)
                    else:
                        self.scale_list[i].config(from_ = 1, to =  50)
                        self.scale_list[i].set(7)
                elif self.w_name[i][0] == 'm':
                    self.scale_list[i].config(from_ = -5, to =  5, resolution = 0.05)
                    self.scale_list[i].set(0)
                elif self.w_name[i][0] == 'b':
                    self.scale_list[i].config(from_ = -5, to =  5, resolution = 0.05)
                    self.scale_list[i].set(0)
                else:
                    self.t.insert(END, "Variable not recognized:", self.w_name[i] )
                self.scale_list[i].bind('<B1-Motion>',self.guessdraw)
        return 0
        
        
    def smoothdata(self):
        self.a.smoooth()
        self.guessdraw(0)
        
        
        
    def update_limits(self,extra):
        #functiontype = self.var_func.get()
        
        self.startfreq = max(float(self.startfreq_text.get()),min(array(self.a.index)))
        self.endfreq = min(float(self.endfreq_text.get()),max(array(self.a.index)))

        self.startidx = self.a.nearest(self.startfreq)
        self.endidx =  self.a.nearest(self.endfreq)        

        for i in range(len(self.guess)):
            self.guess[i] = self.scale_list[i].get()

        self.ax1.cla()
        self.canvas.show()
        self.background = self.canvas.copy_from_bbox(self.ax1.bbox)

        dataplot = self.ax1.plot(array(self.a.index[self.startidx:self.endidx]),self.a.values[self.startidx:self.endidx],animated=True)[0]#(ax = self.ax1).lines[-1]
        x= array(self.a.index[self.startidx:self.endidx])
        guessplot = self.ax1.plot(x,zeros(x.size),animated=True)[0]
        
        
        #self.ax1.set_xlim(min(self.a.index),max(self.a.index))
        #self.ax1.set_ylim(min(self.a.values),1)
        
        self.canvas.show()
        self.guessdraw(0)

        return 0
            
    def open_data(self): 
        global dataplot,guessplot
        self.name = self.DisplaySpectrum()
       
        if self.name == -1:
            return 0
        try:
            self.a = RamanSpectrum(self.name)
            self.normalization_constant= 1/max(self.a)
            self.startfreq = min(array(self.a.index))
            self.startfreq_text.insert(END,str(self.startfreq))
            self.endfreq = max(array(self.a.index))
            self.endfreq_text.insert(END,str(self.endfreq))
            
        except:
            self.t.insert(END, 'error opening spectrum')
            return -1

        self.raw = self.a.copy
        self.ax1.cla()
        self.canvas.show()
        self.background = self.canvas.copy_from_bbox(self.ax1.bbox)
        

        dataplot = self.ax1.plot(array(self.a.index),self.a.values,animated=True)[0]#(ax = self.ax1).lines[-1]
        x= array(self.a.index)
        guessplot = self.ax1.plot(x,zeros(self.a.size),animated=True)[0]
        self.ax1.set_xlim(min(self.a.index),max(self.a.index))
        self.ax1.set_ylim(min(self.a.values),1)
        
        self.canvas.show()
        self.update_limits(0)
        
        return 0
    
    def open_ref(self):

        pass
        
        return 0 
        
    def open_bkg(self):

        pass
        return 0 
  
#    def recalc_data(self):
#        if self.raw.shape == self.reference.shape:
#            self.a= SG_Smooth(self.raw,width = 11,order = 3)/self.reference/self.normalization_constant
#            self.ax1.cla()
#            plot(self.a[0],self.raw/self.reference)
#            plot(self.a[0],self.a[1])
#        else:
#            self.t.insert(END, "reference invalid")
#            self.reference_name = "No IR Reference"
#            
#            self.ax1.cla()
#            self.ax1.plot(self.a[0],self.a[1])
#            self.canvas.draw()
#        
#        return 0
        

    def guessdraw(self,ex):
        global function,guessplot,dataplot

        x = array(self.a.index[self.startidx:self.endidx+1])

        for i in range(len(self.guess)):
            self.guess[i] = self.scale_list[i].get()
 
        self.canvas.restore_region(self.background)
        self.ax1.lines[-1].set_xdata(x)
        self.ax1.lines[-1].set_ydata(function(x,*self.guess))
         
        for l in self.ax1.lines:
            self.ax1.draw_artist(l)
          
        self.canvas.blit(self.ax1.bbox)

        return 0
    def fit_and_draw(self):
        
        global function
        if type(self.a) == pandas.core.series.Series:
            self.a = RamanSpectrum(self.a)
       
        result = fitspectrum(self.a, (self.startfreq,self.endfreq),'Custom',self.guess,function = function)
        if result ==-1:
            self.t.insert(END, "\nError with Fitting ")
            return
        fittingparameters = result[0]
        z = list(fittingparameters[0])   
     
        for i in range(len(self.w_name)):
            self.scale_list[i].set(z[i])
        xfit = result[1]
        yfit = result[2]

        self.canvas.restore_region(self.background)
        self.ax1.lines[-1].set_xdata(xfit)
        self.ax1.lines[-1].set_ydata(yfit)
        
        
        for l in self.ax1.lines:

            self.ax1.draw_artist(l)
        for i in range(len(self.w_name)):
            #self.scale_list[i].set(result[0][i])
            if "w" in self.w_name[i]:
                self.ax1.axvline(x = z[i], ymin = 0, ymax = 1 ,color = 'r')
           
        self.canvas.blit(self.ax1.bbox)

        self.t.insert(END, "  \nResult Found for " + str(self.name))
        self.t.insert(END, " \nReferenced to IR spectrum " +self.reference_name)
        self.t.insert(END, "  \nNormalized by constant "+str(self.normalization_constant))
        for i in range(len(self.w_name)):
            self.t.insert(END, ' \n'+self.w_name[i]+': '+str(z[i]))

        return 0
    def DisplaySpectrum(self):
        import re
        global name_list, str_name_list
        file_opt = options = {}
        options['defaultextension'] = '.spe'
        options['filetypes'] = [('all files', '.*'),('SPE files', '.SPE'), ('csv files','.csv'),('text files','.txt')]
        options['title'] = 'Open Spectrum...'
        options['initialdir'] = '/home/chris/Documents/DataWeiss'

        str_name_list = tkFileDialog.askopenfilename(**options)
        if str_name_list == '':
            return -1
        return str_name_list                        
     



if __name__=='__main__':
    root = Tk()
    root.withdraw()
    DisplayWindow(root,ax=None)
    root.mainloop()