#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 09:53:15 2017

@author: Alexis Klimasewski

compute “B” for a M3 earthquake, and then for every ~M3 earthquake
in your dataset (call it “A”), find A/B for every frequency range you inverted for.  
Then sum up A/B over all frequency ranges and find which earthquake has the minimum value for that
because then that would suggest that earthquake is closest to a brune spectrum
"""

import glob
import numpy as np
import obspy
import matplotlib.pyplot as plt

#magnitude range for inversion constraint
mag_ub = 3.5#2.77
mag_lb = 3.0#2.75
#properties for unit
beta = 3500. #3500m/s
stressdrop = 5e6 #pascals
U = 0.63#0.63
rho = 2750. #kg/m^3

working_dir =  '/Users/aklimase/Desktop/USGS/project/'
event_spectra_dir = working_dir + '/test_codes/Andrews_inversion/'
event_spectra = glob.glob(event_spectra_dir + '[2]*.out')

writefile = 'yes'

#find events in catalog that are in mag range
catalog = working_dir + '/catalogs/all_paths_M2.5_USGS_Catalog.txt'
cat = np.genfromtxt(catalog, comments = '#', delimiter = '|', dtype = None, usecols = [1,10], encoding = None)
event = []
magl = []
for i in range(len(cat)):
    m = cat[i][1]
    if m >= mag_lb and m <= mag_ub:
        magl.append(cat[i][1])
        time = obspy.core.utcdatetime.UTCDateTime(cat[i][0].split('.')[0])
        ev = str(time.year).zfill(4) + '_' + str(time.month).zfill(2) + '_' + str(time.day).zfill(2) + '_' + str(time.hour).zfill(2) + '_' + str(time.minute).zfill(2) + '_' + str(time.second).zfill(2)
        event.append(ev)
       
#compute Brune spectra for all of the events in directory
cf_list = []
cf2_list = []
Brune_list = []
spec_list = []
ev_list = []

spec_demean_list = []
Brune_demean_list = []

spec_demean_list_log = []
Brune_demean_list_log = []
 
for i in range(len(event)):
    f = event_spectra_dir + event[i] + '.out'
    print f
    if f in event_spectra:
        ev_list.append(event[i])
        data = np.genfromtxt(f, dtype = float, comments = '#', delimiter = None, usecols = (0,1))#only read in first two cols
        freq = data[:,0]
        #these record spectra are in m
        spec = (data[:,1])
        ml = magl[i]
        #if less than 3, convert local magnitude to moment magnitude
        if ml < 3.0:
            M = 0.884 + 0.754*ml#0.884 + 0.667*ml, 754
        else:
            M = magl[i]
        #compute Brune in SI units
        #Moment from the moment magnitude
        M0 = 10.**((3./2.)*M + 9.1)
        #corner frequency
        fc = beta*(stressdrop/(8.47*M0))**(1./3.)
        omega0 = (M0*U)/(4.*rho*np.pi*(beta**(3.0)))
        #brune spectra over all frequencies
        Brune = (2.*np.pi*(freq)*omega0)/(1.+((1./fc)*freq)**2.)
        #stay in meters
        cf_list.append(np.log10(spec)-np.log10(Brune))

        Brune_list.append(Brune)
        spec_list.append(spec)

#for each event, find A/B for all other events
#sum up all A/B over freqencies we are fitting
cfarray = np.array(cf_list)
#ind = cfarray.index(0.5)
sum_list = map(sum,cfarray[:,np.arange(27,70)]**2.0) ###found the best fit from 1-32.7Hz

##find the minimum in log space
ind = sum_list.index(min(sum_list))
print(ev_list[ind])
print(min(sum_list))        

print magl[ind]
fig = plt.figure(figsize = (12,10))
plt.ylabel('Velocity amplitude (m)', fontsize = 16)
plt.xlim(0.5,70)
plt.loglog(freq , spec_list[ind], color = 'green', label = 'event spectra')
plt.grid()
plt.loglog(freq, Brune_list[ind], color = 'blue', label = 'Brune spectra')
plt.legend(loc = 'lower right', fontsize = 16)
plt.xlabel('Frequency (Hz)', fontsize = 16)
plt.title(ev_list[ind], fontsize = 16)
plt.tick_params(axis='both', which='major', labelsize=15)
plt.tick_params(axis='both', which='both', length = 5, width = 1)
plt.text(0.7, .1, 'Median log(diff) 1-32.7 Hz (demeaned): ' + str(round(sum_list[ind],3)), fontsize = 16)
plt.show()



#write the constraint file in linear space to agree with the event and station spectra
if writefile == 'yes':
    outfile = open(working_dir + 'test_codes/constraint_' + ev_list[ind] + '.out', 'w')
    out = (np.array([freq, (10.**(cf_list[ind]))]).T)
    outfile.write('#freq_bins \t cf_m \n')
    np.savetxt(outfile, out, fmt=['%E', '%E'], delimiter='\t')
    outfile.close()

