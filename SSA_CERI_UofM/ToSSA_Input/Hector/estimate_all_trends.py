# -*- coding: utf-8 -*-
#
# This program find all files in ./obs_files and estimate all trends.
#
#
# 21/2/2021 Machiel Bos, Santa Clara
##===============================================================================##
#                               Libraries                                        ##
##===============================================================================##
import matplotlib.pyplot as plt
import numpy as np
import os
import time
from datetime import datetime, timedelta
import json
import sys
import re
import pandas as pd
import argparse
from glob import glob
from pathlib import Path
import subprocess
import multiprocessing
import threading
import concurrent.futures
from concurrent.futures import ProcessPoolExecutor, as_completed

def create_step_function(time_vector, jump_times, jump_amplitudes):
    """
    Creates a step function based on the given time vector, multiple jump times, and their corresponding amplitudes.

    Parameters:
    time_vector (list): A list of time points.
    jump_times (list): A list of times at which the step jumps occur.
    jump_amplitudes (list): A list of amplitudes for each step jump.

    Returns:
    list: A list representing the step function.
    """
        # Ensure jump_dates and jump_amps are lists
    if not isinstance(jump_times, list):
        jump_times = [jump_times]
    if not isinstance(jump_amplitudes, list):
        jump_amplitudes = [jump_amplitudes]
    # if len(jump_times) != len(jump_amplitudes):
    #     raise ValueError("The lengths of jump_times and jump_amplitudes must be the same.")
   
    step_function = [0] * len(time_vector)
   
    for jump_time, jump_amplitude in zip(jump_times, jump_amplitudes):
        for i, t in enumerate(time_vector):
            if t >= jump_time:
                step_function[i] += jump_amplitude
                
    return step_function   

       
##===============================================================================##
#                   create Control file for removeoutliers                       ##
##===============================================================================##
def create_removeoutliers_ctl_file(station):
    """ Create ctl file for removeoutlier
    Args:
        station : station name (including _0, _1 or _2) of the mom-file
    """
    # Construct the directory path for the station-specific folder
    directory = Path('./pre_files')
    fname = str(directory / '{0:s}.mom'.format(station))

    #--- Create control.txt file for removeoutliers
    with open("removeoutliers.ctl", "w") as fp:
        fp.write("DataFile              {0:s}.mom\n".format(station))
        fp.write("DataDirectory         ./obs_files\n")
        fp.write("OutputFile            ./{0:s}\n".format(fname))
        fp.write("periodicsignals       365.25 182.625\n")
        fp.write("estimateoffsets       yes\n")
        fp.write("estimatepostseismic   yes\n")
        fp.write("estimateslowslipevent yes\n")
        fp.write("TimeUnit              days\n")
        fp.write("ScaleFactor           1.0\n")
        fp.write("PhysicalUnit          mm\n")
        fp.write("IQ_factor             3\n")
        fp.write("Verbose               no\n")
        fp.write("interpolate           no\n")
##===============================================================================##
#                   create Control file for estimate trend                       ##
##===============================================================================##
def create_estimatetrend_ctl_file(station,noisemodels,useRMLE,noseasonal,phi):
    """ Create estimatetrend.ctl

    Args:
        station (string) : name of station
        noisemodels (string) : PLWN, GGMWN, ...
        useRMLE (boolean): use or not use RMLE option
        noseasonal (boolean): do not include seasonal signal in estimation
        phi (float): some models have phi parameter
    """


    directory = Path('./fin_files')
    fname = str(directory / '{0:s}.mom'.format(station))

    #--- Create control.txt file for EstimateTrend
    fp = open("estimatetrend.ctl", "w")
    fp.write("DataFile            {0:s}.mom\n".format(station))
    fp.write("DataDirectory       ./pre_files\n")
    fp.write("OutputFile          ./{0:s}\n".format(fname))
    #fp.write("estimatemultivariate  yes\n")
    #fp.write("MultiVariateFile    OKGU_SSA.mom\n")
    fp.write("interpolate         no\n")
    fp.write("TimeUnit              days\n")
    fp.write("estimatepostseismic yes\n")
    fp.write("PhysicalUnit        mm\n")
    fp.write("JSON                yes\n")
    fp.write("ScaleFactor         1.0\n")
    fp.write("RandomiseFirstGuess  yes\n")
    if noseasonal==False:
        fp.write("periodicsignals     365.25 182.625\n")
    fp.write("estimateoffsets     yes\n")

#--- Create string with all requested noise models
    combination = ''
    add_small_1mphi = False
    m = re.search('PL',noisemodels)
    if m:
        combination += ' GGM'
        add_small_1mphi = True
    m = re.search('FN',noisemodels)
    if m:
        combination += ' FlickerGGM'
        add_small_1mphi = True
    m = re.search('RW',noisemodels)
    if m:
        combination += ' RandomWalkGGM'
        add_small_1mphi = True
    m = re.search('GGM',noisemodels)
    if m:
        combination += ' GGM'
    m = re.search('WN',noisemodels)
    if m:
        combination += ' White'
    m = re.search('VA',noisemodels)
    if m:
        combination += ' VaryingAnnual'
    m = re.search('AR1',noisemodels)
    if m:
        combination += ' AR1'
    m = re.search('MT',noisemodels)
    if m:
        combination += ' Matern'

    fp.write("NoiseModels         {0:s}\n".format(combination))
    if add_small_1mphi==True:
        fp.write("GGM_1mphi           6.9e-06\n")
    elif phi>0.0:
        fp.write("GGM_1mphi           {0:f}\n".format(phi))
        
    if useRMLE==True:
        fp.write("useRMLE             yes\n")
    else:
        fp.write("useRMLE             no\n")
    fp.write("Verbose               no\n")
    fp.write("OffsetEpochsFile      OffsetEpochs.txt\n")
    fp.close()


##===============================================================================##
#                create Control file for estimate spectrum                       ##
##===============================================================================##
def create_estimatespectrum_ctl_file(station):
    """ Create ctl file for estimatespectrum

    Args:
        station : station name (including _0, _1 or _2) of the mom-file
    """

    #--- Create control.txt file for removeoutliers
    with open("estimatespectrum.ctl", "w") as fp:
        fp.write("DataFile              {0:s}.mom\n".format(station))
        fp.write("DataDirectory         ./fin_files\n")
        fp.write("interpolate           no\n")
        fp.write("ScaleFactor           1.0\n")
        fp.write("TimeUnit              days\n")
        fp.write("PhysicalUnit          mm\n")
        fp.write("WindowFunction       Hann\n")
        fp.write("Fraction               0.1\n")
        fp.write("Verbose               yes\n")
        
##===============================================================================##
#                   create Control file for model spectrum                       ##
##===============================================================================##
def create_model_spectrum_ctl_file(station,noisemodels):
    """ Create ctl file for modelspectrum

    Args:
        station : station name (including _0, _1 or _2) of the mom-file
    """
    noise_models_str = ' '.join(noisemodels)
    #--- Create control.txt file for removeoutliers
    with open("modelspectrum.ctl", "w") as fp:
        fp.write("DataFile              {0:s}.mom\n".format(station))
        fp.write("DataDirectory         ./fin_files\n")
        fp.write("NoiseModels           {0:s}\n".format(noise_models_str))
        fp.write("ScaleFactor           1.0\n")
        fp.write("TimeUnit              days\n")
        fp.write("PhysicalUnit          mm\n")
        fp.write("interpolate           no\n")


##===============================================================================##
#                   Modified Joulian Date (MJD) to decimal year                  ##
##===============================================================================##
def mjd_to_fractional_year(mjd):
    # MJD reference date
    mjd_ref = 2400000.5  # JD for MJD 0
    jd = mjd + mjd_ref
    date = datetime(1858, 11, 17) + timedelta(days=mjd)
    year_start = datetime(date.year, 1, 1)
    year_end = datetime(date.year + 1, 1, 1)
    year_length = (year_end - year_start).days
    year_fraction = (date - year_start).days / year_length
    fractional_year = date.year + year_fraction
    return fractional_year  
            
##===============================================================================##
#          Main funtion for Remove outliers, Estimate trend, and Spectrum        ##
##===============================================================================##
def run_estimate_all_trends(directory):

    print("\n*******************************************")
    print("    estimate_all_trends, version 0.1.9")
    print("*******************************************\n")

    #--- Parse command line arguments in a bit more professional way
    parser = argparse.ArgumentParser(description= 'Estimate all trends')

    #--- List arguments that can be given 
    parser.add_argument('-n', dest='noisemodels', action='store',default='PLWN',
       required=False, help="noisemodel combination (PLWN, FL, etc.)")
    parser.add_argument('-phi', dest='phi', action='store',default='0.0',
       required=False, help="phi parameter in GGM")
    parser.add_argument('-s', dest='station', action='store',default='',
       required=False, help="single station name (without .mom extension)")
    parser.add_argument('-useRMLE', action='store_true',
                                    required=False, help="use RMLE option")
    parser.add_argument('-nograph', action='store_true',
                                    required=False, help="do not create png graph")
    parser.add_argument('-noseasonal', action='store_true',
                                    required=False, help="No seasonal signal")

    args = parser.parse_args()

    #--- parse command-line arguments
    noisemodels = args.noisemodels
    station = args.station
    useRMLE = args.useRMLE
    phi     = float(args.phi)
    noseasonal = args.noseasonal
    nograph = args.nograph

    
    #--- Start the clock!
    start_time = time.time()
    # Remove the last component
    one_level_up = os.path.dirname(directory)
    
    # Remove the next component
    main_dirct = os.path.dirname(one_level_up)
    #--- Read station names in directory ./obs_files
    if len(station)==0:
        # Curnt_directory = os.getcwd()
        # New_directory=os.path.join(Curnt_directory,'Stns_Dir',stn_name)
        os.chdir(directory)
        directory = Path('obs_files')
        fnames = glob(os.path.join(directory, '*.mom'))
       
        #--- Did we find files?
        if len(fnames)==0:
            print('Could not find any mom-file in obs_files')
            sys.exit()
    
        #--- Extract station names
        stations = []
        for fname in sorted(fnames):
            station = Path(fname).stem
            stations.append(station)
    
    else:
        stations = [station]
    
    #--- Does the pre-directory exists?
    if not os.path.exists('pre_files'):
        os.makedirs('pre_files')
    
    #--- Does the mom-directory exists?
    if not os.path.exists('fin_files'):
        os.makedirs('fin_files') 
        
    #------------- Analyse station------#
    output = {}
    # for station in stations:
    
    print(station)

#----------------------------------------------------------------------------------#
#                              Remove Outliers                                     #
#----------------------------------------------------------------------------------#
    #--- Remove outliers    
    create_removeoutliers_ctl_file(station)
   #Get the full path to removeoutliers.py
    
    #removeoutliers_path = os.path.join(main_dirct, 'removeoutliers')
    # if os.path.exists(removeoutliers_path):
    #     # Execute removeoutliers.py using subprocess
    #     subprocess.run([removeoutliers_path])
    # else:
    #     print("removeoutliers.py not found in the same directory as estimate_all_trends.py")
    #path_to_ctlFile = main_dirct + '/' + 'Stns_Dir' + '/'+ station 
    subprocess.run(['removeoutliers'])

#----------------------------------------------------------------------------------#
#                               Estimate Trend                                     #
#----------------------------------------------------------------------------------#
    create_estimatetrend_ctl_file(station,'RWWNFN',False,False,0)
    # Get the full path to removeoutliers.pyestimatespectrum.ctl
    #estimatetrend_path =os.path.join(main_dirct, 'estimatetrend')
    # if os.path.exists(estimatetrend_path):
    #     # Execute removeoutliers.py using subprocess
    #     subprocess.run([estimatetrend_path])
    # else:
    #     print("estimatetrend.py not found in the same directory as estimate_all_trends.py")
    subprocess.run(['estimatetrend'])
    
    #--- parse output
    if os.path.exists('estimatetrend.json')==False:
        print('There is no estimatetrend.json')
        sys.exit()
    try:
        fp_dummy = open('estimatetrend.json','r')
        results = json.load(fp_dummy)
        fp_dummy.close()
    except:
        print('Could not read estimatetrend.json')
        sys.exit()
    output[station] = results
    noisemodels = list(results['NoiseModel'].keys())
    
#----------------------------------------------------------------------------------#
#                            Estimate Spectrum                                     #
#----------------------------------------------------------------------------------#
    #--- Estimate Spectrum
    create_estimatespectrum_ctl_file(station)
    estimatespectrum_path = os.path.join(main_dirct, 'estimatespectrum.py')
    if os.path.exists(estimatespectrum_path):
        # Execute removeoutliers.py using subprocess
        subprocess.run(['python', estimatespectrum_path,'-model','-png'])
    else:
        print("estimatespectrum.py not found in the same directory as estimate_all_trends.py")
    #subprocess.run(['estimatespectrum'])

#----------------------------------------------------------------------------------#
#                              Model spectrum                                     #
#----------------------------------------------------------------------------------#
    #--- Model spectrum  
    #create_model_spectrum_ctl_file(station,'W|F|RW')
    #subprocess.run(['modelspectrum'])
    
#----------------------------------------------------------------------------------#
#                               Writing into Jason                                 #
#----------------------------------------------------------------------------------#
    #--- Save dictionary 'output' as json file
    with open('hector_estimatetrend.json','w') as fp:
        json.dump(output, fp, indent=4)

#----------------------------------------------------------------------------------#
#         SSA Input (taking out the jump and trend from the observation)           #
#----------------------------------------------------------------------------------#
    # Construct the file path
    file_path = f'fin_files/{station}.mom'
    
    # Initialize empty lists for time, observations, and fitted model data
    Mjd_Time = []
    observation_data = []
    fitted_model_data = []
    jump_epch =[]
    header_lines = []  # Store lines starting with '#'
    
    # Read data from .mom file using with open
    with open(file_path, 'r') as file:
        next(file)  # Skip header line if exists
        for line in file:
            # Skip lines that start with '#'
            if line.startswith('#'):
                header_lines.append(line.strip())  # Store full line
                jump_epch.append(float(line.split()[2]))
                continue
            # Split each line into elements based on whitespace (adjust delimiter as needed)
            elements = line.split()
            # Convert elements to floats and append to respective lists
            Mjd_Time.append(float(elements[0]))            # Assuming first column is time mjd
            observation_data.append(float(elements[1]))     # Assuming second column is observation
            fitted_model_data.append(float(elements[2]))    # Assuming third column is fitted model

    # Convert lists to numpy arrays for easier manipulation
    # Apply conversion to all time data
    current_directory = os.getcwd()
    Bob_file = os.path.join(current_directory,"..","..","Bob_my_decyr.txt")
    
    df = pd.read_csv((Bob_file), sep=' ')
    #df = pd.read_csv('Bob_my_decyr.txt', sep=' ')
    Time_fracYr_Bob = df['Var13']
    MjD_Bob = df['Var8']
    start_date = min(Mjd_Time)
    end_date = max(Mjd_Time)
    complete_dates = np.arange(start_date, end_date + 1)
    #indices = np.where(complete_dates == Mjd_Time)[0]
    indices = np.where(np.isin(MjD_Bob, Mjd_Time))[0]
    #indices = np.searchsorted(complete_dates, Mjd_Time)
    filtred_FracYr = Time_fracYr_Bob[indices]
    filtred_FracYr = np.array(filtred_FracYr)


    # filtered_MJD = MjD_Bob[(MjD_Bob >= start_date) & (MjD_Bob <= end_date)]
    # filtered_MJD=np.array(filtered_MJD)
    # indices = np.searchsorted(MjD_Bob, filtered_MJD)
    # filtred_FracYr = Time_fracYr_Bob[indices]
    # filtred_FracYr = np.array(filtred_FracYr)

    time_data_fractional_year = [mjd_to_fractional_year(mjd) for mjd in Mjd_Time]
    if len(time_data_fractional_year) != len(filtred_FracYr):
        print('DECYEAR file: The length of two time vectors are not the same.')
    # Convert lists to numpy arrays for easier manipulation
    timee = filtred_FracYr
    #timee = timee-timee[0]
    y_observed = np.array(observation_data)
    y_fitted = np.array(fitted_model_data)


    #StepFunc_Jump = np.zeros(len(complete_dates)) # Create as a NumPy array

        # Step 1: Extract values from JSON file
    with open('hector_estimatetrend.json', 'r') as file:
        data = json.load(file)
        trend = data[f'{station}']["trend"]*np.array(timee-timee[0])
        jump_sizs = data[f'{station}']['jump_sizes']
        step_funcction = create_step_function(Mjd_Time, jump_epch, jump_sizs)
        #StepFunc_Jump[complete_dates > jump_epch] += jump_sizs[:, np.newaxis]

        Obs_trend_jump = y_observed - trend- step_funcction
    
    # Open the file and write the data
    with open(f'{station}_SSA.dat', 'w') as f:
        for header in header_lines:
            f.write(header + '\n')
        
        for t, td, obs, v in zip(Mjd_Time, timee, y_observed, Obs_trend_jump):
            f.write(f"{t}\t{td:.6f}\t{obs:.6f}\t{v:.6f}\n")
    
    plt.figure(figsize=(12, 10))
    plt.plot(timee, y_observed, 'bo', label='Observations')  # Blue circles for observed data
    plt.plot(timee, y_fitted, 'r-', label='Model',linewidth=3)  # Blue circles for observed data
    for mjd in jump_epch:
        plt.axvline(mjd_to_fractional_year(mjd), color='c', linestyle='--', label='Jump')
    station_base = station[:-2]
    plt.title(f'Observations for {station_base}',fontsize=20)
    plt.xlabel('Time[yr]',fontsize=20)
    # Determine y-axis label based on the last character of the station name
    if station.endswith('_0'):
        ylabel = 'East [mm]'
    elif station.endswith('_1'):
        ylabel = 'North [mm]'
    elif station.endswith('_2'):
        ylabel = 'Up [mm]'
    else:
        ylabel = 'Value'  # Default label if none of the conditions are met
    plt.ylabel(ylabel,fontsize=20)
    plt.legend()
    plt.grid(True)
    plt.savefig(f'Observation and model {station}.png')
    plt.close()

    plt.figure(figsize=(12, 10))
    plt.plot(timee, Obs_trend_jump, 'bo', label='Obs-trend-jump')   # Red line for fitted model
    station_base = station[:-2]
    plt.title(f'Observation without trend and jumps for {station_base}',fontsize=20)
    plt.xlabel('Time[yr]',fontsize=20)
    # Determine y-axis label based on the last character of the station name
    if station.endswith('_0'):
        ylabel = 'East [mm]'
    elif station.endswith('_1'):
        ylabel = 'North [mm]'
    elif station.endswith('_2'):
        ylabel = 'Up [mm]'
    else:
        ylabel = 'Value'  # Default label if none of the conditions are met
    plt.ylabel(ylabel,fontsize=20)
    plt.legend()
    plt.grid(True)
    plt.savefig(f'Obs-trend-jump-{station}.png')
    plt.close()

    plt.figure(figsize=(12, 10))
    plt.plot(timee, y_observed, 'bo', label='Observations')  # Blue circles for observed data
    plt.plot(timee, Obs_trend_jump, 'r-', label='Obs-trend-jump',linewidth=1.5)   # Red line for fitted model
    station_base = station[:-2]
    plt.title(f'Observation Vs Residual without trend and jumps for {station_base}',fontsize=20)
    plt.xlabel('Time[yr]',fontsize=20)
    # Determine y-axis label based on the last character of the station name
    if station.endswith('_0'):
        ylabel = 'East [mm]'
    elif station.endswith('_1'):
        ylabel = 'North [mm]'
    elif station.endswith('_2'):
        ylabel = 'Up [mm]'
    else:
        ylabel = 'Value'  # Default label if none of the conditions are met
    plt.ylabel(ylabel,fontsize=20)
    plt.legend()
    plt.grid(True)
    plt.savefig(f'Obs Vs Obs-trend-jump-{station}.png')
    plt.close()

##===============================================================================##
#                           Running the program                                  ##
##===============================================================================##
t1=time.perf_counter()
os.environ['OMP_NUM_THREADS'] = '1'
#os.environ['MKL_NUM_THREADS'] = '1'
current_directory = os.getcwd()
stn_dir = os.path.join(current_directory, 'Stns_Dir')
station_directories = [os.path.join(stn_dir, d) for d in os.listdir(stn_dir) if os.path.isdir(os.path.join(stn_dir, d))]
station_directories = [d for d in station_directories if not os.path.exists(os.path.join(d, 'hector_estimatetrend.json'))]

#Iterate over each station directory
#for station_direct in station_directories:

#    run_estimate_all_trends(station_direct)


#current_directory = os.getcwd()
#failed_stations_file = os.path.join(current_directory, "failed_stations.txt")

#def run_estimate_all_trends_safe(station_directory):
#    try:
#        run_estimate_all_trends(station_directory)
#    except Exception as e:
#        with open(failed_stations_file, "a") as f:
#            f.write(f"{station_directory} - Error: {str(e)}\n")
#        print(f"Error processing {station_directory}. Logged to file.")

def estimatetrend_parallel():
    # Create a process pool with one process per station
    pool = multiprocessing.Pool(processes=1)
    
    #Map each station directory to a process using the safe wrapper
    
    pool.map(run_estimate_all_trends, station_directories)
    
    # Close the pool to prevent further tasks from being submitted
    pool.close()
    
    # Wait for all processes to finish
    pool.join()

# Call the function to run processes in parallel
if __name__ == '__main__':
    estimatetrend_parallel()
    
t2=time.perf_counter()
print('elapsed time:',t2-t1)