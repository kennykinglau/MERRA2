import orderly_correlations_v8 as correl
import csv
import os
import warnings
import numpy as np
import time

times = []
times.append(time.time())

variables_3d = ['QI','QL','QV','T','CLOUD','OMEGA','RH','TOTAL_WIND']
#variables_3d = ['QI','TOTAL_WIND']
variables_2d = ['TQI','TQL','TQV','TROPT','TS','WIND2M','WIND10M','WIND50M','PS','TS']
#variables_2d = ['TQI','WIND2M']
lowest_merra_freqs = [0,250]
lowest_merra_freqs = [0]
highest_merra_freqs = [750,1000,float('inf')]
highest_merra_freqs = [float('inf')]

stat_tests = ['linregress','Spearman']
#stat_tests = ['Spearman']
FOLDER_NAME = 'InternalPagerTestPlots5'
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
        w.writerow(['Y Descriptor','Y Label', 'X descriptor', 'X Label','Condition Label',statLabels[STAT_TEST],'p-value'])
times.append(time.time())
print("Getting condition dataset")
QUratios = {}
conditions = {}
for T in [1,3]:
    QUratios[T] = correl.get_spt_data('QUratio',T,FOLDER_NAME,
    SPT_FREQ = 220, LOWEST_ELL = 70, HIGHEST_ELL = 250,
    RANDOM = False, PRINT_MESS = True)[0]
    #QUratios[T] = np.ma.masked_invalid(QUratios[T])
    cutoff = np.nanpercentile(QUratios[T],97.5)
    print(cutoff)
    conditions[T] = np.where(QUratios[T] >= cutoff,True,False)

times.append(time.time())
print("Getting MERRA-2 datasets")
big_merras = {}
big_merras[1] = {}
big_merras[3] = {}
DIM = 2
for FINAL_T_INTERVAL in [1,3]:
    for MERRA_VAR in variable_list[DIM]:
        DELTA = False
        for time_type in ['tavg']:
            print(DIM,MERRA_VAR,DELTA,time_type,FINAL_T_INTERVAL)
            merra_data,rangeString,xLabel = correl.get_merra_data(MERRA_VAR,DIM,FOLDER_NAME, FINAL_T_INTERVAL = FINAL_T_INTERVAL, #FINAL_T_INTERVAL only has an effect if going from 1-hour intervals to 3-hour intervals
                DELTA = DELTA,time_type=time_type,
                PRINT_MESS = False)
            big_merras[FINAL_T_INTERVAL][DIM,MERRA_VAR,DELTA,time_type] = merra_data,rangeString,xLabel
            length = merra_data.size
    big_merras[FINAL_T_INTERVAL][DIM,'TIME',DELTA,time_type] = FINAL_T_INTERVAL*np.arange(length),'','TIME'
variables_2d.append('TIME')
variable_list[2] = variables_2d
times.append(time.time())
DIM = 3
for MERRA_VAR in variable_list[DIM]:
    for lowest_merra_freq in lowest_merra_freqs:
        for highest_merra_freq in highest_merra_freqs:
            DELTA = False
            for time_type in ['tavg']:
                print(DIM,MERRA_VAR,lowest_merra_freq,highest_merra_freq,DELTA,time_type)
                merra_data,rangeString,xLabel = correl.get_merra_data(MERRA_VAR,DIM,FOLDER_NAME, #FINAL_T_INTERVAL only has an effect if going from 1-hour intervals to 3-hour intervals
                    MIN_FREQ = lowest_merra_freq,MAX_FREQ=highest_merra_freq,DELTA = DELTA,time_type=time_type, PRINT_MESS = False)
                big_merras[3][DIM,MERRA_VAR,lowest_merra_freq,highest_merra_freq,DELTA,time_type] = merra_data,rangeString,xLabel
                length = merra_data.size
big_merras[3][DIM,'TIME',DELTA,time_type] = 3*np.arange(length),'','TIME'
variables_3d.append('TIME')

print("Doing statistics")
for T_INTERVAL in [1,3]:
    times.append(time.time())
    for STAT_TEST in stat_tests:
        for conditionCode in [0,1,2]:
            if conditionCode==0:
                condition = conditions[T_INTERVAL]
                conditionLabel = 'HighQUratio'
            if conditionCode==1:
                condition = ~conditions[T_INTERVAL]
                conditionLabel = 'LowQUratio'
            if conditionCode==2:
                condition = 'None'
                conditionLabel = ''
            for xKey in big_merras[T_INTERVAL]:
                x_data,xRangeString,xLabel = big_merras[T_INTERVAL][xKey]
                for yKey in big_merras[T_INTERVAL]:
                    if not (T_INTERVAL==3 and xKey[0]==2 and yKey[0]==2):
                        y_data,yRangeString,yLabel = big_merras[T_INTERVAL][yKey]
                        if x_data.size != y_data.size:
                            print("HEY! Sizes don't match: this is weird.")
                        for ABSX in [False]:
                            for ABSY in [False]:
                                USE_RANK = False
                                for LOGX in [True]:
                                    if xLabel=='TIME':
                                        LOGX = False
                                    for LOGY in [True]:
                                        if yLabel=='TIME':
                                            LOGY = False
                                        print(T_INTERVAL,conditionLabel,xRangeString,xLabel,yRangeString,yLabel,USE_RANK,LOGX,LOGY,ABSX,ABSY)
                                        with warnings.catch_warnings():
                                            warnings.simplefilter('ignore')
                                            newyRangeString,newyLabel,newxRangeString,newxLabel,statNumber,pvalue = correl.process_data(x_data,x_data,T_INTERVAL*np.arange(x_data.size),y_data,y_data,T_INTERVAL*np.arange(y_data.size),
                                                xLabel,yLabel,STAT_TEST,FOLDER_NAME,
                                                condition=condition,conditionLabel = conditionLabel,
                                                USE_RANK=USE_RANK, ABSX = ABSX, ABSY = ABSY, LOGX = LOGX, LOGY = LOGY,
                                                xRangeString=xRangeString,yRangeString=yRangeString,PRINT_MESS = False)
                                            with open(os.path.join(folderPath,STAT_TEST+'.csv'),mode='a') as f:
                                                w = csv.writer(f,delimiter=',')
                                                w.writerow([newyRangeString,newyLabel,newxRangeString,newxLabel,conditionLabel,statNumber,pvalue])
                                USE_RANK = True
                                print(T_INTERVAL,conditionLabel,xRangeString,xLabel,yRangeString,yLabel,USE_RANK,LOGX,LOGY,ABSX,ABSY)
                                with warnings.catch_warnings():
                                    warnings.simplefilter('ignore')
                                    newyRangeString,newyLabel,newxRangeString,newxLabel,statNumber,pvalue = correl.process_data(x_data,x_data,T_INTERVAL*np.arange(x_data.size),y_data,y_data,T_INTERVAL*np.arange(y_data.size),
                                        xLabel,yLabel,STAT_TEST,FOLDER_NAME,
                                        condition=condition,conditionLabel = conditionLabel,
                                        USE_RANK=USE_RANK,
                                        xRangeString=xRangeString,yRangeString=yRangeString,PRINT_MESS = False)
                                    with open(os.path.join(folderPath,STAT_TEST+'.csv'),mode='a') as f:
                                        w = csv.writer(f,delimiter=',')
                                        w.writerow([newyRangeString,newyLabel,newxRangeString,newxLabel,conditionLabel,statNumber,pvalue])

times.append(time.time())
#Calculate time differences
time_diffs = [times[i+1] - times[i] for i in range(len(times)-1)]
print(time_diffs)