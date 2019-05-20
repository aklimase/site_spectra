#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 16:20:24 2017

@author: escuser

inputs: path to uncorrected data sac files

cuts for S waves based off arrival time calculated from 3.5 km/s velocity.
add 2 seconds to beginning of sample so that the taper doesn't clip data and cuts 60 seconds after the s arrival time

outputs: writes to cutdata_s directory
"""

from spec_func import cut_swave
import glob
import os.path as path
import os

location = '*'
tsunit = 'VEL'
channel = 'HH*'

top_dir = '/Volumes/USGS_Data/project'
boxpath = top_dir + '/Andrews_inversion/all_paths'

event_dirs = glob.glob(boxpath + '/uncorrected/Event_*')

eventpaths = glob.glob(boxpath + '/uncorrected/Event_*/*.SAC')#full path
print 'Number of files: ', len(eventpaths)

cut_dir = boxpath + '/cutdata_s/'

events = [os.path.basename(x) for x in event_dirs]

#make a directory for each event
for i in range(len(events)):
    if not path.exists(cut_dir + '/' + events[i]):
        os.makedirs(cut_dir + '/' + events[i])

hard_start = 0
#hard_start = event_dirs.index(boxpath + '/uncorrected/Event_2010_01_23_13_17_17')
print(hard_start)

#loop through event directories
for i in range(hard_start, len(event_dirs)):
    #loop through all sac files in directory
    files = glob.glob(event_dirs[i] + '/*.SAC')
    print 'cutting files in event: ' + event_dirs[i]
    for j in range(len(files)):
        base = path.basename(files[j])
        network, stn, channel, space, yyyy, month, day, hh, mm, sssac =  base.split('_')
        ss = sssac.split('.')[0]
        cutfile = cut_dir + events[i] + '/' + network + '_' + stn + '_' + channel + '__' + yyyy + '_'  + month + '_' +  day + '_' +  hh + '_' + mm + '_' + ss + '.SAC'
        #calling cut function
        #cuts 2 sec before s wave arrival and 60 seconds after
        cuttime = cut_swave(files[j], cutfile, 3, 60)