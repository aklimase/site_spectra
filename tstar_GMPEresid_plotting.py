#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 11 17:05:12 2018

@author: temp
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats import pearsonr
import statsmodels
import statsmodels.stats.power
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import MaxNLocator
import matplotlib.ticker as mtick
import matplotlib.lines as mlines


mpl.rcParams['legend.numpoints'] = 1


mpl.rcParams['figure.subplot.hspace'] = 0.3
mpl.rcParams['figure.subplot.left'] = 0.05
mpl.rcParams['figure.subplot.right'] = 0.98
mpl.rcParams['figure.subplot.top'] = 0.95
mpl.rcParams['figure.subplot.bottom'] = 0.05

#mpl.rcParams['xtick.labelsize'] = 10
#mpl.rcParams['ytick.labelsize'] = 10
mpl.rcParams['xtick.major.size'] = 7
mpl.rcParams['ytick.major.size'] = 7
mpl.rcParams['xtick.major.width'] = 2
mpl.rcParams['ytick.major.width'] = 2
mpl.rcParams['xtick.minor.size'] = 7
mpl.rcParams['ytick.minor.size'] = 7
mpl.rcParams['xtick.minor.width'] = 2
mpl.rcParams['ytick.minor.width'] = 2

plt.style.use("classic")



#import glob
#from mpl_defaults import plot_defaults
#plot_defaults()

#use site name as a key to match variables
#top_dir = '/Volumes/USGS_Data/project/'
top_dir = '/Users/aklimase/Desktop/USGS/project/'
out_dir = 'secondo_rebin3_constrained_2013_11_30_11_36_35'
#out_dir = 'secondo_unconstrained'
#out_dir = 'secondo_rebin3_constrained_2010_04_10_22_03_58'
#out_dir = 'secondo_rebin3_constrained_2015_05_06_16_11_01'



#read in tstar files
#tstar_data = np.genfromtxt(top_dir + 'tstar/site/' + out_dir + '/tstar_site.out', delimiter='\t', names = True, dtype = None)#, encoding = None)
tstar_data = np.genfromtxt(top_dir + 'Andrews_inversion/error_analysis/tstar_site_2013_11_30_11_36_35.out', delimiter='\t', names = True, dtype = None)#, encoding = None)

site_t = [tstar_data[i]['site'] for i in range(len(tstar_data))]
t = [tstar_data[i]['tstars'] for i in range(len(tstar_data))]
terr = [tstar_data[i]['tstar_std'] for i in range(len(tstar_data))]
tstar = zip(site_t,t, terr)

#read in vs30, only for our sites
vs30 = np.genfromtxt(top_dir + 'event_boxes/site/v3anza2013_vs30.txt', delimiter='\t', names = True, dtype = None)#, encoding = None)
s = []
vs = []
method = []
for i in range (len(vs30)):
    if vs30[i]['Sta'] in site_t:
#    if (vs30[i]['Sta'] in site_t) and (vs30[i]['Sta'] != 'ERR') and (vs30[i]['Sta'] != 'SOL'):
        s.append(vs30[i]['Sta'])
        vs.append(vs30[i]['Vs30'])
        method.append(vs30[i]['Method'])
vs30 = zip(s,vs, method)
lnvs30 = zip(s, np.log(vs), method)
        
        
#read in site terms, only for our sites
site_terms = np.genfromtxt(top_dir + 'event_boxes/site/site_terms.txt', delimiter='\t', names = True, dtype = None)#, encoding = None)
#site_terms = np.genfromtxt(top_dir + 'Andrews_inversion/error_analysis/tstar_site_2013_11_30_11_36_35.out', delimiter='\t', names = True, dtype = None)#, encoding = None)

s = []
t = []
for i in range (len(site_terms)):
#    if (site_terms[i]['Site'].split(':')[0] in site_t) and (site_terms[i]['Site'].split(':')[0] != 'ERR') and (site_terms[i]['Site'].split(':')[0] != 'SOL'):
    if (site_terms[i]['Site'].split(':')[0] in site_t):
        s.append(site_terms[i]['Site'].split(':')[0])
        t.append(site_terms[i]['Term'])
site_terms = zip(s,t)


#read in spectral levels
#spec_levels = np.genfromtxt(top_dir + 'tstar/site/' + out_dir + '/spec_levels.out', delimiter='\t', names = True, dtype = None)#, encoding = None)
spec_levels = np.genfromtxt(top_dir + 'Andrews_inversion/error_analysis/spec_levels_2013_11_30_11_36_35.out', delimiter='\t', names = True, dtype = None)#, encoding = None)

s = [spec_levels[i]['site'] for i in range(len(spec_levels))]
l = [spec_levels[i]['16Hzm'] for i in range(len(spec_levels))]
m = [spec_levels[i]['614Hzm'] for i in range(len(spec_levels))]
h = [spec_levels[i]['1435Hzm'] for i in range(len(spec_levels))]
l_levels = zip(s,l)
m_levels = zip(s,m)
h_levels = zip(s,h)


#Kilb kappa
Kilb = np.genfromtxt(top_dir + 'event_boxes/site/Kilb_2012_table5.txt', delimiter = None, names = True, dtype = None)#, encoding = None)
s = [Kilb[2*i]['Station'] for i in range(len(Kilb)/2)]
kAH = [np.mean([Kilb[2*i]['kAHms'], Kilb[2*i+1]['kAHms']])/1000. for i in range(len(Kilb)/2)]
ksfix = [np.mean([Kilb[2*i]['kSfixms'], Kilb[2*i+1]['kSfixms']])/1000. for i in range(len(Kilb)/2)]
Kilb_kAH = zip(s,kAH)
Kilb_ksfix  = zip(s,ksfix)

#Anderson kappa
Anderson = np.genfromtxt(top_dir + 'event_boxes/site/Anderson_1991_skappa.txt', delimiter='\t', names = True, dtype = None)#, encoding = None)
s = [Anderson[i]['stn'] for i in range(len(Anderson))]
k = [Anderson[i]['s_wave_kappams']/1000. for i in range(len(Anderson))]
And_k = zip(s,k)

#Joe's data
Joe = np.genfromtxt(top_dir + 'event_boxes/site/Joe_anza_tstar.txt', delimiter='\t', names = True, dtype = None)#, encoding = None)
site = [Joe[i]['stn'] for i in range(len(Joe))]
t = [Joe[i]['tstar'] for i in range(len(Joe))]
l = [Joe[i]['avglcm']/100. for i in range(len(Joe))]
m = [Joe[i]['avglcm']/100. for i in range(len(Joe))]
Joe_t = zip(site, t)
Joe_l = zip(site, l)
Joe_m = zip(site, m)




def main():
#    plot_3(x1 = l_levels, x2 = m_levels, x3 = h_levels, y = tstar ,vs30 = vs30, color = False, x1label = '1-6Hz amplitude (m)', x2label = '6-14Hz amplitude (m)', x3label = '14-35Hz amplitude (m)', ylabel = r'$\kappa_0$ (s)', outfilename = 'kappa_vs_spec_amp')
#    plot_3(x1 = l_levels, x2 = m_levels, x3 = h_levels, y = site_terms ,vs30 = vs30, color = True, x1label = '1-6Hz amplitude (m)', x2label = '6-14Hz amplitude (m)', x3label = '14-35Hz amplitude (m)', ylabel = 'ln(GMPE site residual)', outfilename = 'GMPE_term_vs_spec_amp_rev')
#    plot_3(x1 = And_k, x2 = Kilb_kAH, x3 = Kilb_ksfix, y = tstar ,vs30 = vs30, color = False, x1label = 'Anderson kappa (s)', x2label = 'Kilb kappa AH (s)', x3label = 'Kilb kappa sfix (s)', ylabel = r'$\kappa_0$ (s), this study', outfilename = 'kappa_comp_sameaxes_rev')
#    plot_3(x1 = l_levels, x2 = m_levels, x3 = h_levels, y = vs30 ,vs30 = vs30, color = False, x1label = '1-6Hz amplitude (m)', x2label = '6-14Hz amplitude (m)', x3label = '14-35Hz amplitude (m)', ylabel = 'Vs30 (m/s)', outfilename = 'vs30_vs_spec_amp_rev')
    plot_1(x = tstar, y = site_terms, yerr = , vs30 = vs30, color = True, xlabel = r'$\kappa_0$ (s)', ylabel = 'ln(GMPE site residual)', outfilename = 'GMPE_term_vs_kappa')
#    plot_1(x = tstar, y = vs30, vs30 = vs30, color = False, xlabel = r'$\kappa_0$ (s)', ylabel = 'Vs30 (m/s)', outfilename = 'lnvs30_vs_kappa')
#    plot_1(x = vs30, y = site_terms, vs30 = vs30, color = False, xlabel = 'Vs30 (m/s)', ylabel = 'ln(GMPE site residual)', outfilename = 'GMPE_term_vs_lnvs30_2R')
#    plot_1(x = tstar, y = Joe_t, vs30 = vs30, color = False, xlabel = r'$\kappa_0$ (s)', ylabel = r'$\kappa$ (s) from Joe', outfilename = r'kappa_vs_Joe_tstar')
#    plot_1(x = l_levels, y = Joe_l, vs30 = vs30, color = False, xlabel = '1-6Hz amplitude (m)', ylabel = '1-6Hz amplitude (m) from Joe', outfilename = '1-6Hzlevel')
#    plot_1(x = m_levels, y = Joe_m, vs30 = vs30, color = False, xlabel = '6-14Hz amplitude (m)', ylabel = '6-14Hz amplitude (m) from Joe', outfilename = '6-14Hzlevel')
#
    plot_1(x = tstar, y = site_terms, vs30 = vs30, color = True, xlabel = r'$\kappa_0$ (s)', ylabel = 'ln(GMPE site residual)', outfilename = 'GMPE_term_vs_kappa_No_SOL_ERR')
    plot_3(x1 = l_levels, x2 = m_levels, x3 = h_levels, y = site_terms ,vs30 = vs30, color = True, x1label = '1-6Hz amplitude (m)', x2label = '6-14Hz amplitude (m)', x3label = '14-35Hz amplitude (m)', ylabel = 'ln(GMPE site residual)', outfilename = 'GMPE_term_vs_spec_amp_No_SOL_ERR')

    
#function for a 1x3 plot 
#inputs two lists of variables (x and y) and T/F to color by vs30, 
#generates 3 subplots comparing either spectral levels or diff t*
def plot_3(x1, x2, x3, y, vs30, color, x1label, x2label, x3label, ylabel, outfilename):
#inputs two lists of variables (x and y) and T/F to color by vs30
    mpl.rcParams['font.size'] = 20
    y_list = zip(*y)[1]
    abc = ['(a)', '(b)', '(c)']

#    fig.subplots_adjust(wspace = 0)
#    plt.tight_layout()
    if color == True:
        fig = plt.figure(figsize = (20,7))
        gs = GridSpec(1,4,width_ratios=[1, 1, 1, 0.4])
    else:
        fig = plt.figure(figsize = (18,7))
        gs = GridSpec(1,3)

    fig.text(0.02, 0.5, ylabel, va='center', rotation='vertical')
    
    cmap = mpl.cm.viridis
    norm = mpl.colors.Normalize(vmin=230, vmax=800)
    s_m = mpl.cm.ScalarMappable(cmap = cmap, norm=norm)
    s_m.set_array([])

    
    #first match x 1-3 and y with site name
    for x in [x1,x2,x3]:
        site = []
        x_var = []
        y_var = []
        site_vs30 = []
        for j in range(len(y)):
            for k in range(len(x)):
                for m in range(len(vs30)):
                    if y[j][0] == x[k][0] and y[j][0] == vs30[m][0]:
                        site.append(y[j][0])
                        x_var.append(x[k][1])
                        y_var.append(y[j][1])
                        site_vs30.append(vs30[m][1])

        ind = [x1,x2,x3].index(x)
        ax = plt.subplot(gs.new_subplotspec((0, ind), colspan=1))
        ax.tick_params(axis='both', which='major')
        if ind == 0:
            ax.yaxis.tick_left()
        elif ind == 2:
            ax.yaxis.tick_right()
        else:
            ax.set_yticklabels([''])
        ax.grid(linestyle='--', linewidth=0.5, color = 'gray', which = 'both')
        ax.xaxis.set_major_locator(MaxNLocator(prune='both', nbins=5))
        
        plt.xlim(min(x_var)-0.5*np.abs(min(x_var)), max(x_var)+0.2*max(x_var))
        ###################################################
#        plt.xlim(0,0.08)
#        plt.ylim(0,0.08)
#        plt.xticks([0,0.02,0.04,0.06])
#        plt.yticks([0,0.02,0.04,0.06,0.08])
        
        if y == vs30:
            plt.ylim(min(y_list)-0.4*np.abs(min(y_list)), max(y_list)+0.3*max(y_list))
        elif y == tstar:
            plt.ylim(min(y_list)-0.6*np.abs(min(y_list)), max(y_list)+0.1*max(y_list))
        else:
            plt.ylim(min(y_list)-0.4*np.abs(min(y_list)), max(y_list)+0.6*max(y_list))
#            
        if x1 == And_k:
            plt.plot([0,1,2], [0,1,2], color = 'black', linewidth = 1.5)
#        
#        if y == vs30:
#            plt.yscale('log')
#            ax.set_yticklabels([''])
#            if ind == 0:
#                plt.gca().yaxis.set_major_formatter(mtick.FormatStrFormatter('%.0f')) 
#                plt.gca().yaxis.set_minor_formatter(mtick.FormatStrFormatter('%.0f')) 
#                ax.yaxis.tick_left()
#            if ind == 2:
#                plt.gca().yaxis.set_major_formatter(mtick.FormatStrFormatter('%.0f')) 
#                plt.gca().yaxis.set_minor_formatter(mtick.FormatStrFormatter('%.0f'))
#                ax.yaxis.tick_right()
#                 


        if color == True:
            plt.scatter(x_var, y_var, marker = 'o', s = 60, c = s_m.to_rgba(site_vs30), edgecolors = s_m.to_rgba(site_vs30))
        else:
            if x == vs30:
                for i in range(len(x)):
                    if x[i][2] == 'proxy':
                        plt.scatter(x_var[i], y_var[i], marker = 'o', s = 60, c = 'blue', label = 'proxy')
                    else:
                        plt.scatter(x_var[i], y_var[i], marker = 'x', s = 80, c = 'red', label = 'MASW')
                blue_line = mlines.Line2D([1], [1], color='blue', marker='o', markersize=10, label='proxy', ls = '')
                red_line = mlines.Line2D([1], [1], color='red', marker='x', markersize=10, label='MASW', ls = '')
                plt.legend(handles = [blue_line, red_line], loc = 'lower left', numpoints = 1, fontsize = 18)
            elif y ==vs30:
                for i in range(len(y)):
                    if y[i][2] == 'proxy':
                        plt.scatter(x_var[i], y_var[i], marker = 'o', s = 60, c = 'blue', label = 'proxy')
                    else:
                        plt.scatter(x_var[i], y_var[i], marker = 'x', s = 80, c = 'red', label = 'MASW')
                blue_line = mlines.Line2D([1], [1], color='blue', marker='o', markersize=10, label='proxy', ls = '')
                red_line = mlines.Line2D([1], [1], color='red', marker='x', markersize=10, label='MASW', ls = '')
                plt.legend(handles = [blue_line, red_line], loc = 'lower left', numpoints = 1, fontsize = 18)
            else:
                plt.scatter(x_var, y_var, marker = 'o', s = 60, c = 'blue')
        
        
        for i in range(len(x_var)):#default 2,5
            an = plt.annotate(site[i], xy = (x_var[i], y_var[i]), xytext = (2, 2), textcoords = 'offset points')
            an.draggable()


            
        plt.xlabel([x1label, x2label, x3label][ind])
        
        if y == vs30:
            rval, pval = pearsonr(x_var, np.log(y_var))
            power = statsmodels.stats.power.tt_solve_power(effect_size = rval, nobs = len(x_var), alpha = 0.05)
        else:
            rval, pval = pearsonr(x_var, y_var)
            power = statsmodels.stats.power.tt_solve_power(effect_size = rval, nobs = len(x_var), alpha = 0.05)
        
        label1 = 'Pearson R : ' + "{:.4f}".format(rval)
        label3 = 'power : ' + "{:.4f}".format(power)
        label2 = 'pvalue: ' + "{:.4f}".format(pval)    
                
        ax.annotate(label1 + '\n' + label2 + '\n' + label3, xy=(0.45, 0.04), xycoords='axes fraction', bbox=dict(facecolor='white', edgecolor='black', boxstyle='square,pad=0.2'))
#        ax.annotate(label2, xy=(0.02, 0.07), xycoords='axes fraction')
#        ax.annotate(label3, xy=(0.02, 0.01), xycoords='axes fraction')
        
        ax.annotate(abc[ind], xy = (0.00, 1.02), xycoords='axes fraction', fontsize = 24, weight = 'bold')
        
    
    if color == True:
    ##add color bar to right
        fig.subplots_adjust(right=0.85)
        cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.75])
        cbar = plt.colorbar(s_m, cax=cbar_ax)
        cbar.set_label(ur"Vs30 (m/s)", fontsize = 22)#"$azimuth$ (\u00b0)"
        cbar.ax.tick_params(labelsize = 18)
    
    
    plt.tight_layout()
    fig.subplots_adjust(left = 0.08, top = 0.94, wspace = 0)
    plt.show()
    plt.savefig(top_dir + '/Andrews_inversion/error_analysis/site/' + out_dir + outfilename + '.png')
    
def plot_1(x, y, yerr, vs30, color, xlabel, ylabel, outfilename):
    mpl.rcParams['font.size'] = 22
    mpl.rcParams['xtick.major.size'] = 7
    mpl.rcParams['ytick.major.size'] = 7
    mpl.rcParams['xtick.major.width'] = 2
    mpl.rcParams['ytick.major.width'] = 2
    mpl.rcParams['xtick.minor.size'] = 7
    mpl.rcParams['ytick.minor.size'] = 7
    mpl.rcParams['xtick.minor.width'] = 2
    mpl.rcParams['ytick.minor.width'] = 2


    
    if color == True:
        fig = plt.figure(figsize = (16,12))
    else: 
        fig = plt.figure(figsize = (13,12))
    plt.tight_layout()
    plt.grid(linestyle='--', linewidth=0.5, color = 'gray',  which = 'both')
    cmap = mpl.cm.viridis
    norm = mpl.colors.Normalize(vmin=230, vmax=800)
    s_m = mpl.cm.ScalarMappable(cmap = cmap, norm=norm)
    s_m.set_array([])


    site = []
    x_var = []
    y_var = []
    site_vs30 = []
    for j in range(len(y)):
        for k in range(len(x)):
            for m in range(len(vs30)):
                if y[j][0] == x[k][0] and y[j][0] == vs30[m][0]:
                    site.append(y[j][0])
                    x_var.append(x[k][1])
                    y_var.append(y[j][1])
                    site_vs30.append(vs30[m][1])
    plt.tick_params(axis='both', which='major')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    
    if x == vs30:
        plt.xscale('log')
        plt.gca().xaxis.set_major_formatter(mtick.FormatStrFormatter('%.0f')) 
        plt.gca().xaxis.set_minor_formatter(mtick.FormatStrFormatter('%.0f')) 
    if y == vs30:
        plt.yscale('log')
        plt.gca().yaxis.set_major_formatter(mtick.FormatStrFormatter('%.0f')) 
        plt.gca().yaxis.set_minor_formatter(mtick.FormatStrFormatter('%.0f')) 
        
    plt.xlim(min(x_var)-0.05*np.abs(min(x_var)), max(x_var)+0.05*max(x_var))
    if color == True:
        plt.scatter(x_var, y_var, marker = 'o', s = 60, c = s_m.to_rgba(site_vs30), edgecolors = s_m.to_rgba(site_vs30))
    else:
        if x == vs30:
            for i in range(len(x)):
                if x[i][2] == 'proxy':
                    plt.scatter(x_var[i], y_var[i], marker = 'o', s = 60, c = 'blue', label = 'proxy')
                else:
                    plt.scatter(x_var[i], y_var[i], marker = 'x', s = 80, c = 'red', label = 'MASW')
            blue_line = mlines.Line2D([1], [1], color='blue', marker='o', markersize=10, label='proxy', ls = '')
            red_line = mlines.Line2D([1], [1], color='red', marker='x', markersize=10, label='MASW', ls = '')
            plt.legend(handles = [blue_line, red_line], loc = 'lower left', numpoints = 1)
        elif y ==vs30:
            for i in range(len(y)):
                if y[i][2] == 'proxy':
                    plt.scatter(x_var[i], y_var[i], marker = 'o', s = 60, c = 'blue', label = 'proxy')
                else:
                    plt.scatter(x_var[i], y_var[i], marker = 'x', s = 80, c = 'red', label = 'MASW')
            blue_line = mlines.Line2D([1], [1], color='blue', marker='o', markersize=10, label='proxy', ls = '')
            red_line = mlines.Line2D([1], [1], color='red', marker='x', markersize=10, label='MASW', ls = '')
            plt.legend(handles = [blue_line, red_line], loc = 'lower left', numpoints = 1)

        else:
            plt.scatter(x_var, y_var, marker = 'o', s = 60, c = 'blue')
   
    for i in range(len(x_var)):#default 2,5
        if site[i] == 'KNW':
            an = plt.annotate(site[i], xy = (x_var[i], y_var[i]), xytext = (2, -25), textcoords = 'offset points')
            an.draggable()
        else:
            an = plt.annotate(site[i], xy = (x_var[i], y_var[i]), xytext = (2, 2), textcoords = 'offset points')
            an.draggable()
        
        if y == vs30:
            rval, pval = pearsonr(x_var, np.log(y_var))
            power = statsmodels.stats.power.tt_solve_power(effect_size = rval, nobs = len(x_var), alpha = 0.05)
        elif x == vs30:
            rval, pval = pearsonr(np.log(x_var), y_var)
            power = statsmodels.stats.power.tt_solve_power(effect_size = rval, nobs = len(x_var), alpha = 0.05)
        else:
            rval, pval = pearsonr(x_var, y_var)
            power = statsmodels.stats.power.tt_solve_power(effect_size = rval, nobs = len(x_var), alpha = 0.05)


    label1 = 'Pearson R : ' + "{:.4f}".format(rval)
    label3 = 'power : ' + "{:.4f}".format(power)
    label2 = 'pvalue: ' + "{:.4f}".format(pval)
    
#    plt.annotate(label1 + '\n' + label2 + '\n' + label3, xy=(0.72, 0.02), xycoords='axes fraction', bbox=dict(facecolor='white', edgecolor='black', boxstyle='square,pad=0.2'))
    
##################################
    #masw
    label1a = 'Pearson R : ' + "{:.4f}".format(-0.9025)
    label3a = 'power : ' + "{:.4f}".format(0.2465)
    label2a = 'pvalue: ' + "{:.4f}".format(0.0975)
    #proxy
    label1b = 'Pearson R : ' + "{:.4f}".format(0.0842)
    label3b = 'power : ' + "{:.4f}".format(0.0582)
    label2b = 'pvalue: ' + "{:.4f}".format(0.7947)
    
    plt.annotate(r'proxy' + '\n' +label1b + '\n' + label2b + '\n' + label3b, xy=(0.72, 0.02), xycoords='axes fraction', bbox=dict(facecolor='white', edgecolor='black', boxstyle='square,pad=0.3'))
    plt.annotate(r'proxy', xy=(0.717, 0.1300), xycoords='axes fraction', bbox=dict(facecolor='white', edgecolor='black', boxstyle='square,pad=0.15'))
    plt.annotate('MASW' +'\n' +label1a + '\n' + label2a + '\n' + label3a, xy=(0.42, 0.02), xycoords='axes fraction', bbox=dict(facecolor='white', edgecolor='black', boxstyle='square,pad=0.3'))
    plt.annotate('MASW', xy=(0.417, 0.1300), xycoords='axes fraction', bbox=dict(facecolor='white', edgecolor='black', boxstyle='square,pad=0.15'))

    if color == True:
    ##add color bar to right
#        cbar_ax = fig.add_axes([0.98, 0.1, 0.02, 0.8])
        cbar = plt.colorbar(s_m)
        cbar.set_label(ur"Vs30 (m/s)")#"$azimuth$ (\u00b0)"
        cbar.ax.tick_params(labelsize = 20)
    plt.tight_layout()
    plt.show()
    plt.savefig(top_dir + '/Andrews_inversion/error_analysis/site' + '/' + outfilename + '.png')


main()


y = site_terms
x1 = vs30
x2 = lnvs30
site = []
x1_var = []
x2_var = []
y_var = []


msite = []
mx1_var = []
mx2_var = []
my_var = []

for j in range(len(y)):
    for k in range(len(x1)):
        if y[j][0] == x1[k][0]:
            site.append(y[j][0])
            x1_var.append(x1[k][1])
            x2_var.append(x2[k][1])
            y_var.append(y[j][1])
            if x1[j][2] == 'MASW':
                msite.append(y[j][0])
                mx1_var.append(x1[k][1])
                mx2_var.append(x2[j][1])
                my_var.append(y[j][1])

print 'GMPE site term vs vs30'
rval, pval = pearsonr(x1_var, y_var)
power = statsmodels.stats.power.tt_solve_power(effect_size = rval, nobs = len(x1_var), alpha = 0.05)
label1 = 'Pearson R : ' + "{:.4f}".format(rval)
label3 = 'power : ' + "{:.4f}".format(power)
label2 = 'pvalue: ' + "{:.4f}".format(pval)
print label1 + '\t' + label2 + '\t' + label3

print 'GMPE site term vs lnvs30'
rval, pval = pearsonr(x2_var, y_var)
power = statsmodels.stats.power.tt_solve_power(effect_size = rval, nobs = len(x2_var), alpha = 0.05)
label1 = 'Pearson R : ' + "{:.4f}".format(rval)
label3 = 'power : ' + "{:.4f}".format(power)
label2 = 'pvalue: ' + "{:.4f}".format(pval)
print label1 + '\t' + label2 + '\t' + label3


print 'GMPE site term vs vs30 MASW only'
rval, pval = pearsonr(mx1_var, my_var)
power = statsmodels.stats.power.tt_solve_power(effect_size = rval, nobs = len(mx1_var), alpha = 0.05)
label1 = 'Pearson R : ' + "{:.4f}".format(rval)
label3 = 'power : ' + "{:.4f}".format(power)
label2 = 'pvalue: ' + "{:.4f}".format(pval)
print label1 + '\t' + label2 + '\t' + label3

print 'GMPE site term vs lnvs30 MASW only'
rval, pval = pearsonr(mx2_var, my_var)
power = statsmodels.stats.power.tt_solve_power(effect_size = rval, nobs = len(mx2_var), alpha = 0.05)
label1 = 'Pearson R : ' + "{:.4f}".format(rval)
label3 = 'power : ' + "{:.4f}".format(power)
label2 = 'pvalue: ' + "{:.4f}".format(pval)
print label1 + '\t' + label2 + '\t' + label3 + '\n'


l = [tstar, l_levels, m_levels, h_levels]
names = ['kappa', 'l_levels', 'm_levels', 'h_levels']
for i in range(len(l)):
    x = l[i]
    y1 = vs30
    y2 = lnvs30
    site = []
    x_var = []
    y1_var = []
    y2_var = []
    
    msite = []
    mx_var = []
    my1_var = []
    my2_var = []
    for j in range(len(y1)):
        for k in range(len(x)):
            if y1[j][0] == x[k][0]:
                site.append(y1[j][0])
                x_var.append(x[k][1])
                y1_var.append(y1[j][1])
                y2_var.append(y2[j][1])
                if y1[j][2] == 'proxy':
                    msite.append(y[j][0])
                    mx_var.append(x[k][1])
                    my1_var.append(y1[j][1])
                    my2_var.append(y2[j][1])

    print 'vs30 vs '+ names[i] 
    rval, pval = pearsonr(x_var, y1_var)
    power = statsmodels.stats.power.tt_solve_power(effect_size = rval, nobs = len(x_var), alpha = 0.05)
    label1 = 'Pearson R : ' + "{:.4f}".format(rval)
    label3 = 'power : ' + "{:.4f}".format(power)
    label2 = 'pvalue: ' + "{:.4f}".format(pval)
    print label1 + '\t' + label2 + '\t' + label3
    
    print 'lnvs30 vs ' + names[i]
    rval, pval = pearsonr(x_var, y2_var)
    power = statsmodels.stats.power.tt_solve_power(effect_size = rval, nobs = len(x_var), alpha = 0.05)
    label1 = 'Pearson R : ' + "{:.4f}".format(rval)
    label3 = 'power : ' + "{:.4f}".format(power)
    label2 = 'pvalue: ' + "{:.4f}".format(pval)
    print label1 + '\t' + label2 + '\t' + label3
    
    print 'vs30 vs '+ names[i] + ' proxy only'
    rval, pval = pearsonr(mx_var, my1_var)
    power = statsmodels.stats.power.tt_solve_power(effect_size = rval, nobs = len(mx_var), alpha = 0.05)
    label1 = 'Pearson R : ' + "{:.4f}".format(rval)
    label3 = 'power : ' + "{:.4f}".format(power)
    label2 = 'pvalue: ' + "{:.4f}".format(pval)
    print label1 + '\t' + label2 + '\t' + label3

    print 'lnvs30 vs' + names[i] + ' proxy only'
    rval, pval = pearsonr(mx_var, my2_var)
    power = statsmodels.stats.power.tt_solve_power(effect_size = rval, nobs = len(mx_var), alpha = 0.05)
    label1 = 'Pearson R : ' + "{:.4f}".format(rval)
    label3 = 'power : ' + "{:.4f}".format(power)
    label2 = 'pvalue: ' + "{:.4f}".format(pval)
    print label1 + '\t' + label2 + '\t' + label3 + '\n'
    






