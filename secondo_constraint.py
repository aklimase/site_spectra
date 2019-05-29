#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 15:14:47 2017

@author: Alexis Klimasewski

make a constraint file and put in the box directory
then constrain each event and station

if constraint is an event, divide event and multipy station
if constraint is station, multiply event and divide station

rename the outfile for the event or station

all in m/s
"""

import numpy as np
import glob
import os.path as path

working_dir = '/Users/aklimase/Desktop/USGS/project/test_codes'

#fill in constraint event
secondo_dir = 'Andrews_inversion'
#fill in constraint event name here
constraint ='2010_01_20_11_51_46'
constraint_file =  working_dir + '/constraint_' + constraint + '.out'
outfile_path = working_dir + '/Andrews_inversion_constrained'

con = np.genfromtxt(constraint_file)
cf_spec = con.T[1] #second col
print(cf_spec)
secondo_ev =  glob.glob(working_dir + '/' + secondo_dir + '/2*.out')
secondo_stn = glob.glob(working_dir + '/'+ secondo_dir + '/[!2]*.out')
freq_list = con.T[0] 

##not in log space anymore
for i in range(len(secondo_ev)):#for each event
    #make each row into an array
    event = np.genfromtxt(secondo_ev[i])
    eventid = path.basename(secondo_ev[i]).split('.')[0]
    amp = event.T[1]/cf_spec
    std = event.T[2]/cf_spec
    outfile = open(outfile_path + '/' + eventid + '.out', 'w')
    out = (np.array([freq_list, amp, std])).T
    outfile.write('#freq_bins \t vel_spec_NE_m \t stdev_m \n')

    outfile.close()
    
    
for i in range(len(secondo_stn)):#for each station
    #make each row into an array
    stn = np.genfromtxt(secondo_stn[i])
    stnid = path.basename(secondo_stn[i]).split('.')[0]
    amp = stn.T[1]*cf_spec
    std = stn.T[2]*cf_spec
    outfile = open(outfile_path + '/' + stnid + '.out', 'w')
    out = (np.array([freq_list, amp, std])).T
    outfile.write('#freq_bins \t vel_spec_NE_m \t stdev_m \n')
    np.savetxt(outfile, out, fmt=['%E', '%E', '%E'], delimiter='\t')

    outfile.close()