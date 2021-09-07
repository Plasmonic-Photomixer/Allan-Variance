import csv
import allantools as allan
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import allantools as AT

#### Read in spectrometer data and put all rows into a single 1D vector ####
#### This test had 9999 rows of 160 data points at a rate of 16Hz ####
#### Gives us a total of 99990 seconds of data ####
def load_data(csv_file):
    df=pd.read_csv(csv_file, #skiprows = [1]
    )
    new_df=pd.DataFrame([df.iloc[0].to_list()],columns=df.columns)
    for x in range(int((df.shape[0]-1))):
        sample_df=pd.DataFrame([df.iloc[(x+1)].to_list()],columns=df.iloc[x+1].to_list())
        new_df=pd.concat([new_df,sample_df ],axis=1)
    num_sec = len(new_df.columns)/16        
    timestream = new_df.values[0]
    return (timestream, num_sec)

def plot_timestream(csv_file):
    (timestream, num_sec) = load_data(csv_file)
    plt.plot(timestream)
    plt.show()
##################################################
#                                                #
# Allan Variance Plotting Function:              #
# Inputs: (CSV data file name (in quotes ''),    #
#          A.V. resolution (~30 usually))        #
# Outputs: (Allan vriance plot, Allan time)      #  
#                                                #
##################################################
def allan_plot(csv_file, res):
    (timestream, num_sec) = load_data(csv_file)
    #tau = np.logspace(0, (np.log10(num_sec/5)), res)
    tau = np.logspace(0, 3, res)
    rate = 16 # 1/16 seconds of integration per sample
    
    ##################################################
    #                                                #
    # Overlapping Allan deviation function:          #
    # Inputs: (data, sample rate, data type, taus)   #
    # Outputs: (taus used, deviations, error, pairs) #
    #                                                #
    ##################################################
    
    (tau2, adevs, adev_err, n) = AT.oadev(timestream, rate, data_type="freq", taus=tau)
    avars = np.square(adevs) # square allan dev to get allan var
    # Make white noise t^-1 line
    white_line = (avars[0]*(tau2**-1))      
    correction = np.zeros(len(tau2))
    correction[:len(tau2)-(round(len(tau2)/5-2))] = 1 
    correction[len(tau2)-(round(len(tau2)/5)-2):] = correction[len(tau2)-(round(len(tau2)/5)-2):] + np.arange(1,2.5,1.5/len(correction[len(tau2)-(round(len(tau2)/5)-2):]))
    #Plotting the Allan Variance ####
    plot = plt.loglog(tau2, avars) #Plot the allan variance data
    plt.loglog(tau2,white_line)    #Plot the white noise line for comparison
    #plt.errorbar(tau2, avars, yerr = 2*correction[::]*(avars[::]/np.sqrt(num_sec/tau2[::])), ecolor='g')
    plt.errorbar(tau2, avars, yerr = 2*(avars[::]/np.sqrt((num_sec/tau2[::]))), ecolor='g')
    plt.title('Allan Variance for Lock-in Spectrometer (Bin %d)'%(713))
    plt.xlabel('Integration Times (s)')
    plt.ylabel('Power^2 (arbitrary(?) units)')
    plt.show()
    
    ##### Finding the Allan Time ####
    # Find the deviation from white line percentage and see when it reaches 10%
    
    log_diff_pct = ((np.log10(avars) - np.log10(white_line))/np.log10(white_line))
    diff_pct = ((avars - white_line)/white_line)
    if min(np.abs(log_diff_pct[:] - 0.10)) < 0.1 :
        allan_loc = np.argmin(np.abs(log_diff_pct-0.10))
        allan_time = tau2[allan_loc]
    else: allan_time = np.inf
    if min(np.abs(diff_pct[:]-10)) < 1 :
        allan_loc2 = np.argmin(np.abs(diff_pct-10))
        allan_time2 = tau2[allan_loc2]
        avg_allan_time = (allan_time + allan_time2)/2
        print('The system has an Allan Time of %d seconds'%(avg_allan_time))
    else: 
        allan_time2 = np.inf 
        print('Allan time is E N D L E S S')
            
        
    
    
    
    
    
    