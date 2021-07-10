import orderly_correlations_v6 as correl
import csv
import os
import warnings
import numpy as np
import time

times = []
times.append(time.time())

variables_3d = ['QI','QL','QV','T','CLOUD','OMEGA']
#variables_3d = ['QI']
variables_2d = ['TQI','TQL','TQV','TROPT']
#variables_2d = ['TQI']
lowest_merra_freqs = [0,250]
#lowest_merra_freqs = [0]
highest_merra_freqs = [750,float('inf')]
#highest_merra_freqs = [float('inf')]

spt_variables = ['T','Q','U','QUratio','Q-U']
#spt_variables = ['QUratio']
spt_frequencies = [90,150,220]
spt_frequencies = [220]
lowest_ells = [110,70,30,10]
#lowest_ells = [70]
highest_ells = [150,250,350]
#highest_ells = [250]

stat_tests = ['linregress','Spearman']
#stat_tests = ['Spearman']
FOLDER_NAME = 'MainPagerPlots2'
statLabels = {}
statLabels['linregress'] = 'r^2'
statLabels['Spearman'] = 'rho'
statLabels['ks'] = 'ks'

variable_list = {}
variable_list[2] = variables_2d
variable_list[3] = variables_3d

#This creates the folder where everything will go
dir_path = os.path.dirname(os.path.realpath(__file__))
folderPath = os.path.join(dir_path,FOLDER_NAME)
if not os.path.isdir(folderPath):
    os.mkdir(folderPath)

for STAT_TEST in stat_tests:
    with open(os.path.join(folderPath,STAT_TEST+'.csv'),mode='w') as f:
        w = csv.writer(f,delimiter=',')
        w.writerow(['Y Descriptor','Y Label', 'X descriptor', 'X Label',statLabels[STAT_TEST],'p-value'])

times.append(time.time())
#Generate the data
print("Getting SPT datasets")
big_spt = {}
for SPT_VAR in spt_variables:
    for SPT_FREQ in spt_frequencies:
        for LOWEST_ELL in lowest_ells:
            for HIGHEST_ELL in highest_ells:
                for T_INTERVAL in [1,3]:
                    print(SPT_VAR,T_INTERVAL,FOLDER_NAME,SPT_FREQ,LOWEST_ELL,HIGHEST_ELL)
                    SPT_for_correlation,SPT_data,SPT_time_numbers,rangeString,yLabel = correl.get_spt_data(SPT_VAR,T_INTERVAL,FOLDER_NAME,
                        SPT_FREQ = SPT_FREQ, LOWEST_ELL = LOWEST_ELL, HIGHEST_ELL = HIGHEST_ELL, PRINT_MESS = False)
                    big_spt[SPT_VAR,T_INTERVAL,FOLDER_NAME,SPT_FREQ,LOWEST_ELL,HIGHEST_ELL] = SPT_for_correlation,SPT_data,SPT_time_numbers,rangeString,yLabel

times.append(time.time())
print("Getting MERRA-2 datasets")
big_merra = {}
DIM = 2
for MERRA_VAR in variable_list[DIM]:
    DELTA = False
    for time_type in ['tavg','inst']:
        print(DIM,MERRA_VAR,DELTA,time_type)
        merra_data,rangeString,xLabel = correl.get_merra_data(MERRA_VAR,DIM,FOLDER_NAME, #FINAL_T_INTERVAL only has an effect if going from 1-hour intervals to 3-hour intervals
            DELTA = DELTA,time_type=time_type,
            PRINT_MESS = False)
        big_merra[DIM,MERRA_VAR,DELTA,time_type] = merra_data,rangeString,xLabel
    DELTA = True
    print(DIM,MERRA_VAR,DELTA)
    merra_data,rangeString,xLabel = correl.get_merra_data(MERRA_VAR,DIM,FOLDER_NAME, #FINAL_T_INTERVAL only has an effect if going from 1-hour intervals to 3-hour intervals
        DELTA = DELTA,PRINT_MESS = False)
    big_merra[DIM,MERRA_VAR,DELTA] = merra_data,rangeString,xLabel
times.append(time.time())
DIM = 3
for MERRA_VAR in variable_list[DIM]:
    for lowest_merra_freq in lowest_merra_freqs:
        for highest_merra_freq in highest_merra_freqs:
            DELTA = False
            for time_type in ['tavg','inst']:
                print(DIM,MERRA_VAR,lowest_merra_freq,highest_merra_freq,DELTA,time_type)
                merra_data,rangeString,xLabel = correl.get_merra_data(MERRA_VAR,DIM,FOLDER_NAME, #FINAL_T_INTERVAL only has an effect if going from 1-hour intervals to 3-hour intervals
                    MIN_FREQ = lowest_merra_freq,MAX_FREQ=highest_merra_freq,DELTA = DELTA,time_type=time_type, PRINT_MESS = False)
                big_merra[DIM,MERRA_VAR,lowest_merra_freq,highest_merra_freq,DELTA,time_type] = merra_data,rangeString,xLabel
            DELTA = True
            print(DIM,MERRA_VAR,lowest_merra_freq,highest_merra_freq,DELTA)
            merra_data,rangeString,xLabel = correl.get_merra_data(MERRA_VAR,DIM,FOLDER_NAME, #FINAL_T_INTERVAL only has an effect if going from 1-hour intervals to 3-hour intervals
                MIN_FREQ = lowest_merra_freq,MAX_FREQ=highest_merra_freq,DELTA = DELTA, PRINT_MESS = False)
            big_merra[DIM,MERRA_VAR,lowest_merra_freq,highest_merra_freq,DELTA] = merra_data,rangeString,xLabel

times.append(time.time())
print("Doing statistics")
for STAT_TEST in stat_tests:
    for ABSX in [True,False]:
        USE_RANK = False
        for LOGX in [True,False]:
            for LOGY in [True,False]:
                for SPT_VAR in spt_variables:
                    for SPT_FREQ in spt_frequencies: 
                        for LOWEST_ELL in lowest_ells: 
                            for HIGHEST_ELL in highest_ells:
                                DIM = 2
                                T_INTERVAL = 1
                                SPT_for_correlation,SPT_data,SPT_time_numbers,yRangeString,yLabel = big_spt[SPT_VAR,T_INTERVAL,FOLDER_NAME,SPT_FREQ,LOWEST_ELL,HIGHEST_ELL]
                                for MERRA_VAR in variable_list[DIM]:
                                    DELTA = False
                                    for time_type in ['tavg','inst']:
                                        merra_data,xRangeString,xLabel = big_merra[DIM,MERRA_VAR,DELTA,time_type]
                                        print(STAT_TEST,ABSX,USE_RANK,LOGX,LOGY,SPT_VAR,SPT_FREQ,LOWEST_ELL,HIGHEST_ELL,MERRA_VAR,DELTA,time_type)
                                        newyRangeString,newyLabel,newxRangeString,newxLabel,statNumber,pvalue = correl.process_data(merra_data,merra_data,T_INTERVAL*np.arange(merra_data.size),SPT_for_correlation,SPT_data,SPT_time_numbers,
                                            xLabel,yLabel,STAT_TEST,FOLDER_NAME,
                                            USE_RANK=USE_RANK, ABSX = ABSX, LOGX = LOGX, LOGY = LOGY,
                                            xRangeString=xRangeString,yRangeString=yRangeString,
                                            PRINT_MESS = False)
                                        with open(os.path.join(folderPath,STAT_TEST+'.csv'),mode='a') as f:
                                                w = csv.writer(f,delimiter=',')
                                                w.writerow([newyRangeString,newyLabel,newxRangeString,newxLabel,statNumber,pvalue])
                                    DELTA = True
                                    merra_data,xRangeString,xLabel = big_merra[DIM,MERRA_VAR,DELTA]
                                    print(STAT_TEST,ABSX,USE_RANK,LOGX,LOGY,SPT_VAR,SPT_FREQ,LOWEST_ELL,HIGHEST_ELL,MERRA_VAR,DELTA)
                                    newyRangeString,newyLabel,newxRangeString,newxLabel,statNumber,pvalue = correl.process_data(merra_data,merra_data,T_INTERVAL*np.arange(merra_data.size),SPT_for_correlation,SPT_data,SPT_time_numbers,
                                            xLabel,yLabel,STAT_TEST,FOLDER_NAME,
                                            USE_RANK=USE_RANK, ABSX = ABSX, LOGX = LOGX, LOGY = LOGY,
                                            xRangeString=xRangeString,yRangeString=yRangeString,
                                            PRINT_MESS = False)
                                    with open(os.path.join(folderPath,STAT_TEST+'.csv'),mode='a') as f:
                                            w = csv.writer(f,delimiter=',')
                                            w.writerow([newyRangeString,newyLabel,newxRangeString,newxLabel,statNumber,pvalue])
                                DIM = 3
                                T_INTERVAL = 3
                                SPT_for_correlation,SPT_data,SPT_time_numbers,yRangeString,yLabel = big_spt[SPT_VAR,T_INTERVAL,FOLDER_NAME,SPT_FREQ,LOWEST_ELL,HIGHEST_ELL]
                                for MERRA_VAR in variable_list[DIM]:
                                    for lowest_merra_freq in lowest_merra_freqs:
                                        for highest_merra_freq in highest_merra_freqs:
                                            DELTA = False
                                            for time_type in ['tavg','inst']:
                                                merra_data,xRangeString,xLabel = big_merra[DIM,MERRA_VAR,lowest_merra_freq,highest_merra_freq,DELTA,time_type]
                                                print(STAT_TEST,ABSX,USE_RANK,LOGX,LOGY,SPT_VAR,SPT_FREQ,LOWEST_ELL,HIGHEST_ELL,MERRA_VAR,lowest_merra_freq,highest_merra_freq,DELTA,time_type)
                                                newyRangeString,newyLabel,newxRangeString,newxLabel,statNumber,pvalue = correl.process_data(merra_data,merra_data,T_INTERVAL*np.arange(merra_data.size),SPT_for_correlation,SPT_data,SPT_time_numbers,
                                                    xLabel,yLabel,STAT_TEST,FOLDER_NAME,
                                                    USE_RANK=USE_RANK, ABSX = ABSX, LOGX = LOGX, LOGY = LOGY,
                                                    xRangeString=xRangeString,yRangeString=yRangeString,
                                                    PRINT_MESS = False)
                                            with open(os.path.join(folderPath,STAT_TEST+'.csv'),mode='a') as f:
                                                w = csv.writer(f,delimiter=',')
                                                w.writerow([newyRangeString,newyLabel,newxRangeString,newxLabel,statNumber,pvalue])
                                            DELTA = True
                                            merra_data,xRangeString,xLabel = big_merra[DIM,MERRA_VAR,lowest_merra_freq,highest_merra_freq,DELTA]
                                            print(STAT_TEST,ABSX,USE_RANK,LOGX,LOGY,SPT_VAR,SPT_FREQ,LOWEST_ELL,HIGHEST_ELL,MERRA_VAR,lowest_merra_freq,highest_merra_freq,DELTA)
                                            newyRangeString,newyLabel,newxRangeString,newxLabel,statNumber,pvalue = correl.process_data(merra_data,merra_data,T_INTERVAL*np.arange(merra_data.size),SPT_for_correlation,SPT_data,SPT_time_numbers,
                                                xLabel,yLabel,STAT_TEST,FOLDER_NAME,
                                                USE_RANK=USE_RANK, ABSX = ABSX, LOGX = LOGX, LOGY = LOGY,
                                                xRangeString=xRangeString,yRangeString=yRangeString,
                                                PRINT_MESS = False)
                                            with open(os.path.join(folderPath,STAT_TEST+'.csv'),mode='a') as f:
                                                w = csv.writer(f,delimiter=',')
                                                w.writerow([newyRangeString,newyLabel,newxRangeString,newxLabel,statNumber,pvalue])
                                                
        USE_RANK = True
        for SPT_VAR in spt_variables:
            for SPT_FREQ in spt_frequencies: 
                for LOWEST_ELL in lowest_ells: 
                    for HIGHEST_ELL in highest_ells:
                        DIM = 2
                        T_INTERVAL = 1
                        SPT_for_correlation,SPT_data,SPT_time_numbers,yRangeString,yLabel = big_spt[SPT_VAR,T_INTERVAL,FOLDER_NAME,SPT_FREQ,LOWEST_ELL,HIGHEST_ELL]
                        for MERRA_VAR in variable_list[DIM]:
                            DELTA = False
                            for time_type in ['tavg','inst']:
                                merra_data,xRangeString,xLabel = big_merra[DIM,MERRA_VAR,DELTA,time_type]
                                print(STAT_TEST,ABSX,USE_RANK,LOGX,LOGY,SPT_VAR,SPT_FREQ,LOWEST_ELL,HIGHEST_ELL,MERRA_VAR,DELTA,time_type)
                                newyRangeString,newyLabel,newxRangeString,newxLabel,statNumber,pvalue = correl.process_data(merra_data,merra_data,T_INTERVAL*np.arange(merra_data.size),SPT_for_correlation,SPT_data,SPT_time_numbers,
                                    xLabel,yLabel,STAT_TEST,FOLDER_NAME,
                                    USE_RANK=USE_RANK, ABSX = ABSX, LOGX = LOGX, LOGY = LOGY,
                                    xRangeString=xRangeString,yRangeString=yRangeString,
                                    PRINT_MESS = False)
                                with open(os.path.join(folderPath,STAT_TEST+'.csv'),mode='a') as f:
                                        w = csv.writer(f,delimiter=',')
                                        w.writerow([newyRangeString,newyLabel,newxRangeString,newxLabel,statNumber,pvalue])
                            DELTA = True
                            merra_data,xRangeString,xLabel = big_merra[DIM,MERRA_VAR,DELTA]
                            print(STAT_TEST,ABSX,USE_RANK,LOGX,LOGY,SPT_VAR,SPT_FREQ,LOWEST_ELL,HIGHEST_ELL,MERRA_VAR,DELTA)
                            newyRangeString,newyLabel,newxRangeString,newxLabel,statNumber,pvalue = correl.process_data(merra_data,merra_data,T_INTERVAL*np.arange(merra_data.size),SPT_for_correlation,SPT_data,SPT_time_numbers,
                                    xLabel,yLabel,STAT_TEST,FOLDER_NAME,
                                    USE_RANK=USE_RANK, ABSX = ABSX, LOGX = LOGX, LOGY = LOGY,
                                    xRangeString=xRangeString,yRangeString=yRangeString,
                                    PRINT_MESS = False)
                            with open(os.path.join(folderPath,STAT_TEST+'.csv'),mode='a') as f:
                                    w = csv.writer(f,delimiter=',')
                                    w.writerow([newyRangeString,newyLabel,newxRangeString,newxLabel,statNumber,pvalue])
                        DIM = 3
                        T_INTERVAL = 3
                        SPT_for_correlation,SPT_data,SPT_time_numbers,yRangeString,yLabel = big_spt[SPT_VAR,T_INTERVAL,FOLDER_NAME,SPT_FREQ,LOWEST_ELL,HIGHEST_ELL]
                        for MERRA_VAR in variable_list[DIM]:
                            for lowest_merra_freq in lowest_merra_freqs:
                                for highest_merra_freq in highest_merra_freqs:
                                    DELTA = False
                                    for time_type in ['tavg','inst']:
                                        merra_data,xRangeString,xLabel = big_merra[DIM,MERRA_VAR,lowest_merra_freq,highest_merra_freq,DELTA,time_type]
                                        print(STAT_TEST,ABSX,USE_RANK,LOGX,LOGY,SPT_VAR,SPT_FREQ,LOWEST_ELL,HIGHEST_ELL,MERRA_VAR,lowest_merra_freq,highest_merra_freq,DELTA,time_type)
                                        newyRangeString,newyLabel,newxRangeString,newxLabel,statNumber,pvalue = correl.process_data(merra_data,merra_data,T_INTERVAL*np.arange(merra_data.size),SPT_for_correlation,SPT_data,SPT_time_numbers,
                                            xLabel,yLabel,STAT_TEST,FOLDER_NAME,
                                            USE_RANK=USE_RANK, ABSX = ABSX, LOGX = LOGX, LOGY = LOGY,
                                            xRangeString=xRangeString,yRangeString=yRangeString,
                                            PRINT_MESS = False)
                                    with open(os.path.join(folderPath,STAT_TEST+'.csv'),mode='a') as f:
                                        w = csv.writer(f,delimiter=',')
                                        w.writerow([newyRangeString,newyLabel,newxRangeString,newxLabel,statNumber,pvalue])
                                    DELTA = True
                                    merra_data,xRangeString,xLabel = big_merra[DIM,MERRA_VAR,lowest_merra_freq,highest_merra_freq,DELTA]
                                    print(STAT_TEST,ABSX,USE_RANK,LOGX,LOGY,SPT_VAR,SPT_FREQ,LOWEST_ELL,HIGHEST_ELL,MERRA_VAR,lowest_merra_freq,highest_merra_freq,DELTA)
                                    newyRangeString,newyLabel,newxRangeString,newxLabel,statNumber,pvalue = correl.process_data(merra_data,merra_data,T_INTERVAL*np.arange(merra_data.size),SPT_for_correlation,SPT_data,SPT_time_numbers,
                                        xLabel,yLabel,STAT_TEST,FOLDER_NAME,
                                        USE_RANK=USE_RANK, ABSX = ABSX, LOGX = LOGX, LOGY = LOGY,
                                        xRangeString=xRangeString,yRangeString=yRangeString,
                                        PRINT_MESS = False)
                                    with open(os.path.join(folderPath,STAT_TEST+'.csv'),mode='a') as f:
                                        w = csv.writer(f,delimiter=',')
                                        w.writerow([newyRangeString,newyLabel,newxRangeString,newxLabel,statNumber,pvalue])

times.append(time.time())
#Calculate time differences
time_diffs = [times[i+1] - times[i] for i in range(len(times)-1)]
print(time_diffs)