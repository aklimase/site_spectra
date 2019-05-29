compute_spectra_meters
Uses mtspec to calculate spectra binned into 75 bins equally spaced in log frequency space
Input: cut and corrected Sac files 
Output: writes output files and event directories with the frequency bins, spectral amplitude bins, and 1 sigma standard deviation

Notes: Cut and corrected sac files should be stored in directories named for event time (Event_yyyy_mm_dd_hh_mm_ss). Files in each directory should be named with network, station, channel, and event time (NN_STN_CHN__yyyy_mm_dd_hh_mm_ss.SAC. If files are named differently modify naming conventions in lines 50-66. Spectral binning starts 0.1 Hz and bins into 75 bins.
Spectra files are saved in record_spectra/Event*/ directories.
Uses functions from spec_func.py. add to python path or keep in same directory as run directory

secondo_meters
Runs Andrews inversion in each frequency bin
Input: reads in each binned record spectra and uncertainties 
Output: writes event and station spectra and uncertainties in Andrews_inversion directory

findBrune
Finds the most "Brune-like" spectra of all events between chosen magnitudes
Input: reads in magnitudes from USGS catalog file and event spectra from 'Andrews_inversion'
Output: change write file parameter to 'yes' to write a constraint file in the main directory for choose event

Note: we choose to constrain for best Brune event between 1 and 35 Hz (index 27, 70) determined by noise and filtering. Change indices in line 92 to change frequencies. 

secondo_constraint
constraints site and event spectra from Andrews_inversion
Input: reads in files from 'Andrews_inversion' and the constraint file from findBrune
Output: writes site and event files to 'Andrews_inversion_constraint'

tstar_site
Reads in constrained site spectra and invest for kappa and average spectral levels
Inputs: reads in site spectra from 'Andrews_inversion_constrained'
Outputs: writes files kappa_site.out and spec_levels.out

spec_func
Contains functions bin_spec and bin_max_err used in compute_spectra_meters. Add to python path of put in the run directory