#!/usr/bin/env python
import time as tU
import datetime
from optparse import OptionParser
import os,sys
import shlex

def StringtoOption(string_in):
    # After trying some homemade tricks, here is an efficient way, taken from : 
    # https://stackoverflow.com/questions/12010126/how-to-parse-an-arbitrary-option-string-into-a-python-dictionary
    args = shlex.split(string_in)
    options = {k: True if v.startswith('-') else v
           for k,v in zip(args, args[1:]+["--"]) if k.startswith('-')}
    return options

def farmit(options_pt, n=10):
    """
    options_pt: dictionary containing the predict_Tsky.py options
    n: number of days in a batch
    """

    # get arguments that will be used for the farmfiles, not sure why the original script chose gndData
    site=options_pt['-l']
    gndData=options_pt['-g']
    dates=options_pt.get('-d').split(',')
    
    # reproduce the argument line for predict_Tsky.py
    arg_str=' '.join([options_pt.keys()[i]+' '+options_pt.values()[i] for i in range(len(options_pt.keys()))])

    if len(dates)==1:
        start=dates[0]
        end=start 
    elif len(dates)==2:
        start=dates[0]
        end=dates[1]
    else:
        print "Too many arguments in the time string after -d"
        sys.exit()
    
    tstart = datetime.datetime.strptime(start, '%Y%m%d')
    tend = datetime.datetime.strptime(end, '%Y%m%d')
    tdelta = tend - tstart
    dateList = [tstart + datetime.timedelta(days = i) for i in range(tdelta.days + 1)]
    dateBatches = [dateList[i:i+n] for i in range(0,len(dateList),n)]

    dataDir = 'merra2_products'
    
    # creates farmfiles directory if it doesn't exist
    farmdir='../farmfiles'
    if not os.path.exists(farmdir):
        os.mkdir(farmdir)
    maindir='../'
    if os.path.isfile(maindir+'merra2_python/predict_Tsky.py')==False:
        print "Seems like %s doesn't exist"%(maindir+'merra2_python/predict_Tsky.py')
        sys.exit()

    count = 0 #ensure that not too many MERRA2 requests are being made at once--run 100 at a time, at most.
    for batch in dateBatches:
        count += 1
        s = batch[0].strftime("%Y%m%d") #start date of a batch
        e = batch[-1].strftime("%Y%m%d") #end date of a batch
        dstr = "%s,%s"%(s,e)
        f = open("%s/%s_%s_%s_gndData%s.sh"%(farmdir,site, s, e,gndData),'w')
        f.write("#!/bin/bash\n")
        f.write('#SBATCH -n 1   # Number of cores \n')
        f.write('#SBATCH -N 1   # all cores are on one machine \n')
        f.write('#SBATCH -t %d # Runtime in minutes \n'%(n*20))
        f.write('#SBATCH -p serial_requeue,itc_cluster,kovac    # Partition to submit \n')
        f.write('#SBATCH --mem=1000  # Memory pool for all cores \n')
        f.write('#SBATCH -o %s/%s_%s_%s_gndData%s_%s.out  #STDOUT to file \n'%(farmdir,site,s,e,gndData,'%j'))
        f.write('#SBATCH -e %s/%s_%s_%s_gndData%s_%s.err  #STDERR to file \n'%(farmdir,site,s,e,gndData,'%j'))
        f.write("#SBATCH --job-name=MERRA_%s_%s \n"%(site,s))
        f.write('\n')
        f.write('\n')
        f.write('cd %s \n'%maindir)
        
        f.write('module load python/2.7.14-fasrc02\n')        
        f.write('source activate ~dbarkats/.conda/envs/mypyenv2 \n')
        f.write('./merra2_python/predict_Tsky.py %s'%arg_str) 
        f.close()
        cmd = "sbatch "+f.name
        print cmd
        os.system(cmd)
        # The original had a 5 second wait here, not sure why. I know I had some issues with submitting many jobs too fast
        # before, so I think it is safe to leave some time for the file to be saved to disk
        tU.sleep(5)
        if count%20==0:
            print "20 jobs currently running... paused for 5 minutes"
            tU.sleep(600)





if __name__ == "__main__":
    usage = '''
    ./farmit.py -n NumDayPerBatch -c 'predict_Tsky.py usual_predict_Tsky_options'
    For instance :
    ./farmit.py -n 2 -c 'predict_Tsky.py -l SouthPole -d 20180121,20180129 -g 2'
    '''

    parser = OptionParser(usage=usage)
    
    parser.add_option("-n",
                      dest="n",
                      default = 10,
                      type=int,
                      help="-n number of days to process in each batch. default = 10")
    
    parser.add_option("-c",
                      dest="cmd",
                      type=str,
                      help="-c option: predict_Tsky command. For instance 'predict_Tsky.py -l SouthPole -d 20180129,20180129 -g 2', see help from predict_Tsky.py")


    (options, args) = parser.parse_args()
    n=options.n

    if options.cmd == None:
        print " ### Error: You must specify a start datetime with the -c option."
        exit()
    else:   
        options_pt=StringtoOption(options.cmd)

    if '-d' not in options_pt:
        print " ### Error: You must specify a datetime with the -d option in the command."
        exit()
    
    print n
    farmit(options_pt, n)
