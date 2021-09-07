import pandas as pd
import allantools as AT
import numpy as np
import matplotlib.pyplot as plt

###########################################################
#                                                         #
#  ~~~~~~~~~~~~~~~ Instructions for use ~~~~~~~~~~~~~~~   #                                  
#                                                         #
# The only function that matters is the allan_plot func.  #
# Run this program in an ipython terminal:                # 
#                (run hoh_variance.py)                    #
# Then call allan_plot with inputs of:                    #         
# (csv file, channel number, tau resolution, sample rate) #
# Read the info above the function for more details       #
#                                                         #           
# The allan_plot function will output:                    #
#       1. An allan variance plot compared to tau^-1      #
#       2. A value for the allan time defined by when     #
#       the allan variance deviates from the white        # 
#       noise line by 10%.                                #
#                                                         #
# There is also spectro_allan_plot which does the same as #
# allan_plot except takes the allan variance and time of  #
# a differential between two bins.                        #
#                                                         # 
# Inputs to spectro_allan_plot are:                       #
# (csv file, central channel number, off center channel,  #
#  tau resolution, sample rate)                           #  
#                                                         #
# PATCH 1.0:                                              #
# Added smoothing functions to take spikes out            #
# simply call smooth(_spectro)_allan_plot to get          #
# an allan plot from de-spiked data                       #  
###########################################################


#########################################
#                                       #
# Data formatting for readout:          #
# timestamp, runtime, count(secs), bins #
# 4096 bins in spectrometer             #
# But we throw away 2 edge bins         # 
# Total of 4094 bins in data            #
#                                       #
#########################################

#### Initialize global variables ####
skips = 100 # Rows are skipped to throw out junk data from adrians spectrometer
dead_cols = 3 # number of useless columns in dataset

#########################################
#                                       #
# Function to load in csv file:         #
# Input: csv file name in quotes ('')   #
#                                       #
#########################################

def load_data(csv_file):

    #### Read in csv file and chop off useless rows and colums for counting ####
    df = pd.read_csv(csv_file, skiprows=[i for i in range(skips+2)], usecols=range(dead_cols,4094)
    )
    num_bin = len(df.columns)
    num_sec = len(df)

    # Total of N rows of useful info 
    # Columns labeled and referenced as 'Bin XXXX' where X is between 1 and N
    #### Give all of the columns a title ####
    titles = []
    for i in range(num_bin):
        titles.append('Bin %d'%(i))
    df.columns = titles
    return(df, num_bin, num_sec)


#### Plotting a time series ####
def plot_timestream(csv_file, chan):
    df = load_data(csv_file)[0]
    col = df['Bin %d'%(chan)]
    col.plot()
    #### Uncomment these lines to set plot limits for higher resolution #####
    #plt.ylim(0,2*10**12)
    #plt.xlim(0, 2000)
    plt.title('Timestream of Data(Bin %d)'%(chan))
    plt.show()    

#### Data Smoothing Function ####
def smooth_data(dataframe, chan):
    t_stream = dataframe['Bin %d'%(chan)].values
    sig = np.std(t_stream)
    med = np.median(t_stream)
    count = 1
    print('sigma equals %d'%(sig))
    while count < (len(t_stream) - 1):
        dev = np.abs(t_stream[count] - med)
        #print('time %d has a dev of %d'%(count, dev))
        if (dev > (0.5*sig)):
            t_stream[count] = med
        count += 1
    return t_stream    
    
#### Time Stream but with Smoothing ####            
def plot_smoothstream(csv_file, chan):
    df = load_data(csv_file)[0]
    t_stream = smooth_data(df, chan)
    plt.plot(t_stream)
    plt.show()
        
    
    

# Make sure to use overlapping alan variance to make all of the data pretty

##################################################
#                                                #
# Allan Variance plot function:                  #
# Inputs: (csv file name in quotes (''),         #    
#          channel number [something central],   #   
#          tau resolution [usually around 30],   #
#          sample rate [1 sec in tests])         #
# Outputs: (taus used, deviations, error, pairs) #
#                                                #
################################################## 

def allan_plot(csv_file, chan, res, rate):

    (df, num_bin, num_sec) = load_data(csv_file)
    tau = np.logspace(0, (np.log10(num_sec/5)), res)  # Create a list of taus from 1 to N seconds
    #cen_bin = smooth_data(df, chan)
    cen_bin = df['Bin %d'%(chan)].values
                                     
    ##################################################
    #                                                #
    # Overlapping Allan deviation function:          #
    # Inputs: (data, sample rate, data type, taus)   #
    # Outputs: (taus used, deviations, error, pairs) #
    #                                                #
    ################################################## 
    (tau2, adevs, adev_err, n) = AT.oadev(cen_bin, rate, data_type="freq", taus=tau)
    avars = np.square(adevs)  #square allan dev to get allan var

    # Now make white noise set for comparison line with power density equivalent to the data values
    # actually we finna cheat and just plot t^-1 line :P
    white_line = (avars[0]*(tau2**-1))
    
    #### Plotting the Allan Variance ####
    plot = plt.loglog(tau2, avars) #Plot the allan variance data
    plt.loglog(tau2,white_line)    #Plot the white noise line for comparison
    #### Uncomment the next line if you want to see (kinda sketchy) error bars ####
    #plt.errorbar(tau2, avars, yerr = (np.square(adev_err)), ecolor='g')
    plt.title('Allan Variance for SVII (RX on)(Bin %d)(%d samples tossed)'%(chan, skips))
    plt.xlabel('Integration Times (s)')
    plt.ylabel('(Probably) Power^2')
    plt.show()
    
    ##### Finding the Allan Time ####
    # Find the deviation from white line percentage and see when it reaches 10%
    
    log_diff_pct = ((np.log10(avars) - np.log10(white_line))/np.log10(white_line))
    diff_pct = ((avars - white_line)/white_line)
    allan_loc = np.argmin(np.abs(log_diff_pct-0.10))
    allan_loc2 = np.argmin(np.abs(diff_pct-10))
    allan_time = tau2[allan_loc]
    allan_time2 = tau2[allan_loc2]
    avg_allan_time = (allan_time + allan_time2)/2

    
    print('The system has an Allan Time of %d seconds'%(avg_allan_time))
    
def smooth_allan_plot(csv_file, chan, res, rate):

    (df, num_bin, num_sec) = load_data(csv_file)
    tau = np.logspace(0, (np.log10(num_sec/5)), res)  # Create a list of taus from 1 to N seconds
    cen_bin = smooth_data(df, chan)
    #cen_bin = df['Bin %d'%(chan)].values
                                     
    ##################################################
    #                                                #
    # Overlapping Allan deviation function:          #
    # Inputs: (data, sample rate, data type, taus)   #
    # Outputs: (taus used, deviations, error, pairs) #
    #                                                #
    ################################################## 
    (tau2, adevs, adev_err, n) = AT.oadev(cen_bin, rate, data_type="freq", taus=tau)
    avars = np.square(adevs)  #square allan dev to get allan var

    # Now make white noise set for comparison line with power density equivalent to the data values
    # actually we finna cheat and just plot t^-1 line :P
    white_line = (avars[0]*(tau2**-1))
    
    #### Plotting the Allan Variance ####
    plot = plt.loglog(tau2, avars) #Plot the allan variance data
    plt.loglog(tau2,white_line)    #Plot the white noise line for comparison
    #### Uncomment the next line if you want to see (kinda sketchy) error bars ####
    #plt.errorbar(tau2, avars, yerr = (np.square(adev_err)), ecolor='g')
    plt.title('Allan Variance for SVII (Smoothed Data)(RX on)(Bin %d)(%d samples tossed)'%(chan, skips))
    plt.xlabel('Integration Times (s)')
    plt.ylabel('(Probably) Power^2')
    plt.show()
    
    ##### Finding the Allan Time ####
    # Find the deviation from white line percentage and see when it reaches 10%
    
    log_diff_pct = ((np.log10(avars) - np.log10(white_line))/np.log10(white_line))
    diff_pct = ((avars - white_line)/white_line)
    allan_loc = np.argmin(np.abs(log_diff_pct-0.10))
    allan_loc2 = np.argmin(np.abs(diff_pct-10))
    allan_time = tau2[allan_loc]
    allan_time2 = tau2[allan_loc2]
    avg_allan_time = (allan_time + allan_time2)/2

    
    print('The system has an Allan Time of %d seconds'%(avg_allan_time))    
    
def spectro_allan_plot(csv_file, cen_chan, off_chan, res, rate):
    (df, num_bin, num_sec) = load_data(csv_file)
    tau = np.logspace(0, (np.log10(num_sec/5)), res)  # Create a list of taus from 1 to N seconds
    cen_bin = df['Bin %d'%(cen_chan)].values
    off_bin = df['Bin %d'%(off_chan)].values   
    #cen_bin = smooth_data(df, cen_chan)
    #off_bin = smooth_data(df, off_chan)
                                     
    ##################################################
    #                                                #
    # Overlapping Allan deviation function:          #
    # Inputs: (data, sample rate, data type, taus)   #
    # Outputs: (taus used, deviations, error, pairs) #
    #                                                #
    ################################################## 
    (tau2, adevs, adev_err, n) = AT.oadev((off_bin/cen_bin), rate, data_type="freq", taus=tau)
    avars = np.square(adevs)  #square allan dev to get allan var

    # Now make white noise set for comparison line with power density equivalent to the data values
    # actually we finna cheat and just plot t^-1 line :P
    white_line = (avars[0]*(tau2**-1))
    
    #### Plotting the Allan Variance ####
    plot = plt.loglog(tau2, avars) #Plot the allan variance data
    plt.loglog(tau2,white_line)    #Plot the white noise line for comparison
    #### Uncomment the next line if you want to see (kinda sketchy) error bars ####
    #plt.errorbar(tau2, avars, yerr = (np.square(adev_err)), ecolor='g')
    plt.title('Allan Variance for SVII (Bin %d) divided by (Bin %d)(%d samples tossed)'%(off_chan, cen_chan, skips))
    plt.xlabel('Integration Times (s)')
    plt.ylabel('(Probably) Power^2')
    plt.show()
    
    ##### Finding the Allan Time ####
    # Find the deviation from white line percentage and see when it reaches 10%
    
    log_diff_pct = ((np.log10(avars) - np.log10(white_line))/np.log10(white_line))
    diff_pct = ((avars - white_line)/white_line)
    allan_loc = np.argmin(np.abs(log_diff_pct-0.10))
    allan_loc2 = np.argmin(np.abs(diff_pct-10))
    allan_time = tau2[allan_loc]
    allan_time2 = tau2[allan_loc2]
    avg_allan_time = (allan_time + allan_time2)/2

    
    print('The system has an Allan Time of %d seconds'%(avg_allan_time))  


def smooth_spectro_allan_plot(csv_file, cen_chan, off_chan, res, rate):
    (df, num_bin, num_sec) = load_data(csv_file)
    tau = np.logspace(0, (np.log10(num_sec/5)), res)  # Create a list of taus from 1 to N seconds
    #cen_bin = df['Bin %d'%(cen_chan)].values
    #off_bin = df['Bin %d'%(off_chan)].values   
    cen_bin = smooth_data(df, cen_chan)
    off_bin = smooth_data(df, off_chan)
                                     
    ##################################################
    #                                                #
    # Overlapping Allan deviation function:          #
    # Inputs: (data, sample rate, data type, taus)   #
    # Outputs: (taus used, deviations, error, pairs) #
    #                                                #
    ################################################## 
    (tau2, adevs, adev_err, n) = AT.oadev((off_bin/cen_bin), rate, data_type="freq", taus=tau)
    avars = np.square(adevs)  #square allan dev to get allan var

    # Now make white noise set for comparison line with power density equivalent to the data values
    # actually we finna cheat and just plot t^-1 line :P
    white_line = (avars[0]*(tau2**-1))
    
    #### Plotting the Allan Variance ####
    plot = plt.loglog(tau2, avars) #Plot the allan variance data
    plt.loglog(tau2,white_line)    #Plot the white noise line for comparison
    #### Uncomment the next line if you want to see (kinda sketchy) error bars ####
    #plt.errorbar(tau2, avars, yerr = (np.square(adev_err)), ecolor='g')
    plt.title('Allan Variance for SVII (Bin %d) divided by (Bin %d)(Smoothed Data)(Rx On)(%d samples tossed)'%(off_chan, cen_chan, skips))
    plt.xlabel('Integration Times (s)')
    plt.ylabel('(Probably) Power^2')
    plt.show()
    
    ##### Finding the Allan Time ####
    # Find the deviation from white line percentage and see when it reaches 10%
    
    log_diff_pct = ((np.log10(avars) - np.log10(white_line))/np.log10(white_line))
    diff_pct = ((avars - white_line)/white_line)
    allan_loc = np.argmin(np.abs(log_diff_pct-0.10))
    allan_loc2 = np.argmin(np.abs(diff_pct-10))
    allan_time = tau2[allan_loc]
    allan_time2 = tau2[allan_loc2]
    avg_allan_time = (allan_time + allan_time2)/2

    
    print('The system has an Allan Time of %d seconds'%(avg_allan_time))      
    