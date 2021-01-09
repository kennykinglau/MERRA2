#! /usr/bin/env python


from optparse import OptionParser
import dateutil.parser as dparser
import datetime
from pylab import *
import urllib2
import merra2Player as m2p
import pickle
import glob

#This is the poor man's fourier transform, to see the signal power at different timescales. 
if __name__ == '__main__':
    usage = '''
  
    
    '''

    import merra2Player as m2p
    from pylab import *
    parser = OptionParser(usage=usage)
    
    parser.add_option("-l",
                      dest="location",
                      default = 'None',
                      type=str,
                      help="-l option comma-separated string to define site. Can be either -s SouthPole or -s NewSite, lat, lon, alt")
    
    parser.add_option("-s",
                      dest="start",
                      type=str,
                      help="-s option: Year to produce a plot for.")
    ##TODO: Detect when the daterange goes out of the range of MERRA2, and switch to using GFS. 
    
    parser.add_option("-e",
                      dest="end",
                      type=str,
                      help="-e option: year end.")

    parser.add_option("-g",
                      dest="groundData",
                      default = 2,
                      type=int,
                      help="-g Type of surface level data. 1: historical measured data. 2 (default): MERRA2 surface data. 3: Extrapolation from MERRA2 3d data.")

    parser.add_option("-c",
                      dest="cloud",
                      default = "vapor",
                      type=str,
                      help="-c Type of cloud data in the am profile. vapor: only h2o vapor in the am layers. icewater: ice and water content included.")

    (options, args) = parser.parse_args()
    currentSite=options.location
    gndData = options.groundData
    cloud = options.cloud
    siteInitials = {"SouthPole":"SP","ChajnantorPlateau":"CP","ChajnantorCerro":"CC","MaunaKea":"MK"}
    try:
        siteInitial = siteInitials[currentSite]
    except:
        print "Invalid site. Choose among SouthPole, ChajnantorPlateau, ChajnantorCerro, and MaunaKea."
        exit()

    tstart = options.start
    tend = options.end
    if tstart == None:
        print " ### Error: You must specify a start datetime with the -s option."
        exit()
    if tend == None:
        tend = tstart

    transmission = 0.8

    for year in range(int(tstart),int(tend)+1):
        filename = glob.glob("/n/home10/anwang16/merra2_work/merra2_products/%s/final/%s_%s*_gndData%s_%s*"%(currentSite,currentSite,year,gndData,cloud))[0]
        data = pickle.load(open(filename,'rb'))
        for key in data.keys():
            data[key] = np.array(data[key])

        #All analysis using resampled data.
        datetimes = data['t_tipper']
        residualTsky = data['Tsky_tipper']/transmission - data['tipper850_interp'] #Tsky residual. window transmission must be accounted for
        averages = [] #collection of the intervalAverages below.

        averages += [[mean(residualTsky)]]
        residualTsky = residualTsky - mean(residualTsky)

        tstart = datetime.datetime(year,1,1)
        tend = datetime.datetime(year+1,1,1)
        intervals = [90.,30.,10.,1.,1.0/4]
        #intervals = [90.,30.,10.,1.,1.0/4,1.0/24,1.0/144] #interval sizes in number of days

        for interval in intervals:
            numIntervals = int(360./interval)
            dateList = [tstart + datetime.timedelta(days = i*interval) for i in range(numIntervals+1)] # creates a list of separating dates.
            intervalAverages = [] #for each specific interval: over the year, over 3 months, over 1 day, etc.
            if type(data['t_tipper'][0]) != datetime.datetime:
                continue
            for i in range(numIntervals):
                intervalFilter = (data['t_tipper'] > dateList[i]) & (data['t_tipper'] < dateList[i+1])
                residuals = residualTsky[intervalFilter] #only get the residuals in the interval.
                intervalAverages += [mean(residuals)]
                residualTsky[intervalFilter] = residualTsky[intervalFilter]-mean(residuals)
            print len(dateList[:-1]), len(intervalAverages)
            averages += [intervalAverages]

        with open("/n/home10/anwang16/merra2_work/merra2_products/analysis/TimescaleAnalysis_%s_%s_gndData%s_%s.csv"%(currentSite,year,gndData,cloud),"w") as f:
            for line in averages:
                line = ['%2.3f'%i for i in line]
                f.write(','.join(line)+"\n")
            f.close()

