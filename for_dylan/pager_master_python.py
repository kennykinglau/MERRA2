import orderly_correlations_v3 as correl
import csv
import os
import warnings

variables_3d = ['QI','QL','QV']
#variables_3d = ['QI']
variables_2d = ['TQI','TQL','TQV']
#variables_2d = ['TQI']
spt_variables = ['T','Q','U','QUratio','Q-U']
#spt_variables = ['QUratio']
spt_frequencies = [90,150,220]
spt_frequencies = [220]
lowest_ells = [110,70,30,10]
#lowest_ells = [70]
highest_ells = [150,250,350]
#highest_ells = [250]
stat_tests = ['linregress','Spearman','ks']
#stat_tests = ['Spearman']
FOLDER_NAME = 'PagerPlots2'
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
        w.writerow(['Lowest ell','Highest ell', 'Frequency [GHz]','Y Label','X Label',statLabels[STAT_TEST],'p-value'])
        
#Generate the data
print("Getting SPT datasets")
big_spt = {}
for SPT_VAR in spt_variables:
    for SPT_FREQ in spt_frequencies:
        for LOWEST_ELL in lowest_ells:
            for HIGHEST_ELL in highest_ells:
                for DIM in [2,3]:
                    print(SPT_VAR,SPT_FREQ,LOWEST_ELL,HIGHEST_ELL,DIM)
                    SPT_for_correlation,SPT_data,SPT_time_numbers = correl.get_spt_data(SPT_VAR,DIM,FOLDER_NAME,
                    SPT_FREQ, LOWEST_ELL, HIGHEST_ELL,
                    RANDOM=False, PRINT_MESS=False, #Change
                    )
                    big_spt[SPT_VAR,SPT_FREQ,LOWEST_ELL,HIGHEST_ELL,DIM] = SPT_for_correlation,SPT_data,SPT_time_numbers
print("Getting MERRA-2 datasets")
big_merra = {}
for DIM in [2,3]:
    for DELTA in [True,False]:
        for MERRA_VAR in variable_list[DIM]:
            print(DIM,DELTA,MERRA_VAR)
            big_merra[DELTA,DIM,MERRA_VAR] = correl.get_merra_data(MERRA_VAR,DIM,FOLDER_NAME,
                DELTA=DELTA,
                PRINT_MESS=False
                ) #Change

print("Doing statistics")
for SPT_VAR in spt_variables: #U
    for SPT_FREQ in spt_frequencies: 
        for LOWEST_ELL in lowest_ells: #50
            for HIGHEST_ELL in highest_ells: #150
                for DIM in [2,3]: #2
                    for DELTA in [True,False]: #True
                        for MERRA_VAR in variable_list[DIM]: #TQV
                            for STAT_TEST in stat_tests:
                                for LOGX in [True,False]:
                                    for LOGY in [True,False]:
                                        print(SPT_VAR,SPT_FREQ,LOWEST_ELL,HIGHEST_ELL,DIM,DELTA,MERRA_VAR,STAT_TEST,LOGX,LOGY)
                                        merra_data = big_merra[DELTA,DIM,MERRA_VAR]
                                        with warnings.catch_warnings():
                                            warnings.simplefilter("ignore")
                                            (SPT_for_correlation,SPT_data,SPT_time_numbers) = big_spt[SPT_VAR,SPT_FREQ,LOWEST_ELL,HIGHEST_ELL,DIM]
                                            yLabel,xLabel,statNumber,pvalue = correl.process_data(merra_data,SPT_for_correlation,SPT_data,SPT_time_numbers,
                                                SPT_VAR,MERRA_VAR,DIM,STAT_TEST,FOLDER_NAME,
                                                color_data = False, INVY = False, LOGX = LOGX, LOGY = LOGY, DELTA = DELTA,
                                                FREQ_FILTER='None',SPT_FREQ=SPT_FREQ,LOWEST_ELL=LOWEST_ELL,HIGHEST_ELL=HIGHEST_ELL,
                                                RANDOM = False, MAKE_IMAGES = False, NUM_OF_SAMPS=300, PRINT_MESS = False, #Change
                                                )
                                        with open(os.path.join(folderPath,STAT_TEST+'.csv'),mode='a') as f:
                                            w = csv.writer(f,delimiter=',')
                                            w.writerow([LOWEST_ELL,HIGHEST_ELL,SPT_FREQ,yLabel,xLabel,statNumber,pvalue])