#!/usr/bin/ipython
# -*- coding: utf-8 -*-
"""
Created on Fri May 23 13:51:48 2014

@author: mkk
"""

import pylab
import numpy , sys
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection

from pylab import rcParams
rcParams['figure.figsize'] = 15, 3

font = {'family' : 'Sans',
        'weight' : 'normal',
        'size'   : 10}

matplotlib.rc('font', **font)

def index_min_max(data,col):
    i_min = 0
    i_max = 0
    range_iter = 0
    dsize = numpy.size(data[:,0])
    for i in range(dsize):
        if ( ( int(data[i,col]) > -100 ) and ( range_iter == 0 ) ):            
            i_min = i
            range_iter = 2
        if ( ( int(data[i,col]) < -99 ) and ( range_iter == 2 ) ):            
            i_max = i
            range_iter = 3    
    return i_min,i_max


def plot_rel(filename,polar_name):
    data           = numpy.loadtxt(filename)
    data_polar     = numpy.loadtxt(filename[0:len(filename)-4]+polar_name)
    mink = 0.0
    maxk = 0.0
    range_iter = 0
    dsize = numpy.size(data[:,0])
    for i in range(dsize):
        if ( ( int(data[i,2]) != 0 ) and ( range_iter == 0 ) ):            
            mink = data[i,0]
            range_iter = 2
        if ( ( int(data[i,2]) == 0 ) and ( range_iter == 2 ) ):            
            maxk = data[i,0]
            range_iter = 3
            
  
    max_modow = int(data[dsize/2,2]);
    print "max modow=",max_modow
    
    #fig1 = plt.figure()
    for i in range(max_modow):
        
        i_min,i_max = index_min_max(data_polar,i+3)
        x = data[i_min:i_max,0]        
        y = data[i_min:i_max,i+3]        
        z = data_polar[i_min:i_max,i+3]    
        points   = np.array([x, y]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        lc = LineCollection(segments, cmap=plt.get_cmap('jet'),
            norm=plt.Normalize(-1, 1))
        lc.set_array(z)
        lc.set_linewidth(2)        
        plt.gca().add_collection(lc)
        
        
    i_min,i_max = index_min_max(data_polar,3)    
    x = data[i_min:i_max,0]
    y = data[i_min:i_max,1]
    plt.plot(x,y,label="$E_F$",lw=2,c='r')
    
    plt.legend()
    plt.xlabel("k [1/nm]")
    #plt.ylabel("E [meV]")
    plt.xlim(mink,maxk)  


params = sys.argv
if( numpy.size(params) < 2 ):
    print "Zla liczba argumentow:"
    print "./relacja.py nazwa_pliku.txt"
else:
    print "Parametry wejsciowe wejsciowy:", str(sys.argv)
    print "Plik wejsciowy:", params[1]
    #print "Liczba modow :", params[2]
    Z = [[0,0],[0,0]]        
    levels = pylab.linspace(-1,1,10)
    CS3 = plt.contourf(Z, levels, cmap=plt.get_cmap('jet'))
    plt.clf()      
    
    filename = params[1]
    plt.subplot(1,3,1)    
    
    plt.title("Spin z")
    plot_rel(filename,"_polar_z.txt")
    plt.colorbar(CS3)

    plt.subplot(1,3,3)
    plt.title("Spin x")
    plot_rel(filename,"_polar_x.txt")
    plt.colorbar(CS3)

    plt.subplot(1,3,2)
    plt.title("Spin y")
    plot_rel(filename,"_polar_y.txt")
    
    plt.colorbar(CS3)

    
#plt.ylim(ymax=data[0,1]*1.5) 
#pylab.title(filename+" kolumna: "+params[2])
#pylab.savefig(filename + "-col=" + params[2] + '.png')

plt.show();
