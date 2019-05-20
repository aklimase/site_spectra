
#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  2 15:47:34 2017

@author: Alexis Klimasewski

input: reads in all of the site output files from secondo


"""

import numpy as np
import glob
import matplotlib.pyplot as plt
import os.path as path
import matplotlib.lines as mlines
import matplotlib as mpl
plt.style.use("classic")
mpl.rcParams['font.size'] = 30

from scipy.optimize import curve_fit

#freq bin number
f1 = 28 #28
f2 = 70

secondo_dir = 'Andrews_inversion_constrained'

working_dir = '/Users/aklimase/Desktop/USGS/project/'
site_dir =  working_dir + secondo_dir

site_list = glob.glob(site_dir + '/[!2]*.out')
site_list = sorted(site_list)

outfilename = working_dir + 'tstar_site.out'
spec_levels_file =  working_dir + 'spec_levels.out'

shiftind = -1#-1

def decay_func(f, A0, k):
    return np.log(A0) + (-np.pi*f*k)


def calc_site_tstar(f1,f2):
    I = len(site_list)
    J = 75
    residual = np.zeros((len(site_list), J))
    
    site = np.zeros((I, J))

    print 'number of sites: ', I
    tstar = []
    A = []
    L2norm = []
    tstar_std = []
    A_std = []

    cmap = plt.get_cmap('hsv')
    colors = cmap(np.linspace(0, 1.0, I))
    
    ###################################
    avg_l = []
    avg_m = []
    avg_h = []
    site = []
    
    avg_std_l = []
    avg_std_m = []
    avg_std_h = []

    #for each site, calculate the 1-6, 6-14, and 14-36 Hz spectral levels
    for j in range(I):
        data = np.genfromtxt(site_list[j], comments = '#', dtype = float)
        freq = data.T[0]
        Amplitude = (data.T[1]) #data from log to lin
        std = data.T[2]
#        print std

        l = []
        m = []
        h = []
        l_std = []
        m_std = []
        h_std = []
    
        for i in range(len(freq)):
            ########changed to 35#############
            if 14.0<freq[i]<=35.0:
                h.append(Amplitude[i])
                h_std.append(std[i])
            if 6.0<freq[i]<=14.0:
                m.append(Amplitude[i])
                m_std.append(std[i])
            if 1.0<freq[i]<=6.0:
                l.append(Amplitude[i])
                l_std.append(std[i])
        site.append(path.basename(site_list[j]).split('.')[0])
        avg_l.append(np.mean(l))
        avg_m.append(np.mean(m))
        avg_h.append(np.mean(h))
        
        avg_std_l.append(np.sqrt(sum([x**2 for x in l_std])))
        avg_std_m.append(np.sqrt(sum([x**2 for x in m_std])))
        avg_std_h.append(np.sqrt(sum([x**2 for x in h_std])))
        
        
    outfile = open(spec_levels_file, 'w')
    out = (np.array([site, avg_l,avg_std_l,  avg_m, avg_std_m, avg_h, avg_std_h]).T)
    outfile.write('#site \t 1-6Hz(m) \t l_std(m) \t  6-14Hz(m) \t m_std(m) \t 14-35Hz(m) \t h_std(m) \n')
    np.savetxt(outfile, out, fmt='%s', delimiter='\t')
    outfile.close()

    #for each site
    for i in range(I):
        data = np.genfromtxt(site_list[i], comments = '#', dtype = float)
        freq = data.T[0][f1:f2]
        Amplitude = (data.T[1][f1:f2]) #data from log to lin
        
        residual[i:,] = (data.T[1])#residual in matrix, lin
        
        G = np.zeros((len(freq), 2))
        d = np.zeros((len(freq), 1))
        std = (data.T[2][f1:f2])
        
        l = site_list[i].split('/')[-1].split('.')[0]
        
        plt.plot(data.T[0][f1:f2], data.T[1][f1:f2], color = colors[i], linestyle = '-', lw = 2, alpha= 1, label = l)#y = linear A
        for j in range(len(freq)):
            d[j][0] = np.log(Amplitude[j])#lin to ln
            G[j][0] = -1*np.pi*(freq[j])
            G[j][1] = 1
        
        G_inv = np.linalg.pinv(G, rcond=1e-10)
        m = np.dot(G_inv,d)
        t = float(m[0])
        A0 = np.exp(float(m[1]))
        tstar.append(float(t))
        A.append(A0)
        covd = np.diag(((std)/(Amplitude))**2.)
        covm = np.dot((np.dot(G_inv, covd)), G_inv.T)
        var = (covm.diagonal())
        tstar_std.append(np.sqrt(np.abs(var[0])))#error in kappa
        A_std.append(np.sqrt(np.abs(A0*var[1])))



    outfile = open(outfilename, 'w')
    names = [path.basename(site_list[i]).split('.')[0] for i in range(len(site_list))]
    out = (np.array([names, tstar, tstar_std, A, A_std])).T


    outfile.write('#site \t tstar(s) \t tstar_std \t A0(m) \t A0_std \n')
    np.savetxt(outfile, out, fmt='%s', delimiter='\t')
    outfile.close()

    return freq[0], freq[-1], np.mean(L2norm)


calc_site_tstar(f1,f2)