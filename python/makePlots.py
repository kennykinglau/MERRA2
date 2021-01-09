#! /usr/bin/env python


from optparse import OptionParser
import dateutil.parser as dparser
from collections import Counter
import datetime
import pickle
import os
import glob

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

    m = m2p.merra2Player(site = currentSite)
    #for i in range(1996,2018):
    #    print i, m.findFailedSpectra(str(i))

    #Deleting empty MERRA2 files and re-downloading them
    #delta = datetime.timedelta(days = 1)
    #start = datetime.datetime(1996,1,1)
    #end = datetime.datetime(2017,6,1)
    #date = start
    #while date < end:
    #    datestr = date.strftime("%Y%m%d")
    #    try:
    #        merrafile = glob.glob("/n/home10/anwang16/merra2_work/merra2_products/%s/merra2_raw_data/*inst3_3d*%s*"%(currentSite,datestr))[0]
    #        if os.stat(merrafile).st_size < 10000:
    #            os.remove(merrafile)
    #            m.retrieve_merra2_data_for_date(datestr)
    #            m.retrieve_merra2_data_for_date(datestr, dataset = 'singleLevel')
    #            time.sleep(10)
    #        date += delta
    #    except:
    #        m.retrieve_merra2_data_for_date(datestr)
    #        m.retrieve_merra2_data_for_date(datestr, dataset = 'singleLevel')
    #        time.sleep(10)
    #        date += delta


    #monthly plots.
    #year = 1998
    #for month in range(5,6):
    #    if month < 10:
    #        month = "0"+str(month)
    #    else:
    #        month = str(month)

    #annual plots.
    month = ""
    for year in range(int(tstart),int(tend)+1):
        data = {}
        keys = ['t','pwv', 'BK100','BK150','BK220','BK270','tipper850','t_tipper','averaged_Tsky_tipper','Tsky_sigma','Tsky_tipper','tau','Tatm','tau850','tipper850_interp','tipper_unc']
        for i in keys:
            data[i] = np.array([])

        #files = ['SouthPole_20170101_20170110_gndData2.p', 'SouthPole_20170210_20170219_gndData2.p', 'SouthPole_20170322_20170331_gndData2.p', 'SouthPole_20170501_20170510_gndData2.p', 'SouthPole_20170111_20170120_gndData2.p', 'SouthPole_20170220_20170301_gndData2.p', 'SouthPole_20170401_20170410_gndData2.p', 'SouthPole_20170511_20170520_gndData2.p', 'SouthPole_20170121_20170130_gndData2.p', 'SouthPole_20170302_20170311_gndData2.p', 'SouthPole_20170411_20170420_gndData2.p', 'SouthPole_20170521_20170530_gndData2.p', 'SouthPole_20170131_20170209_gndData2.p', 'SouthPole_20170312_20170321_gndData2.p', 'SouthPole_20170421_20170430_gndData2.p', 'SouthPole_20170531_20170601_gndData2.p']
        #files.sort()
        #Chajnantor
        for i in glob.glob("/n/home10/anwang16/merra2_work/merra2_products/%s/final/%s_%s%s*_gndData%s_%s*"%(currentSite,currentSite,year,month,gndData,cloud)):
        #for i in files:
            new = pickle.load(open(i,'rb'))
            if len(new.keys()) == 6:
                print i
                os.remove(i)
                continue
            for key in data.keys():
                try:
                    data[key] = np.append(data[key],new[key])
                except:
                    data[key] += new[key]

        #accounting for the window. propagated error is already added in merra2Player and unc. 
        transmission = 0.8
        data['Tsky_tipper'] = data['Tsky_tipper']/transmission
        data['averaged_Tsky_tipper'] = data['averaged_Tsky_tipper']/transmission
        data['Tatm'] = data['Tatm']/transmission


        #for saving to a text file (Txtdata array)
        years = [x.year for x in data['t']]
        months = [x.month for x in data['t']]
        days = [x.day for x in data['t']]
        hours = [x.hour for x in data['t']]
        
        dtype = [('year',str),('month',str),('day',str),('hour',str),('BK100 Tsky (Krj)',float),('BK150 Tsky (Krj)',float),('BK220 Tsky (Krj)',float),('BK270 Tsky (Krj)',float),('PWV (um)',float)]
        Txtdata = np.array([years,months,days,hours,data['BK100'],data['BK150'],data['BK220'],data['BK270'],data['pwv']])
        Txtdata = np.transpose(Txtdata)
        Txtdata = np.array(Txtdata)
        #np.sort(data,order = ["year","month","day","hour"])
        np.savetxt("%s_BK_MERRA2.txt"%year,Txtdata, fmt = "%.4i %.2i %.2i %.2i %3.3f %3.3f %3.3f %3.3f %3.3f", header = "Year, Month, Day, Hour (UTC), BK100 Tsky (K_rj), BK150 Tsky (K_rj), BK220 Tsky (K_rj), BK270 Tsky (K_rj), PWV (um)")

        figure(figsize = (35,15))
        scatter(data['t'],data['BK270'],color='black',label = '270 GHz')
        scatter(data['t'],data['BK220'],color='red', label = '220 GHz')
        scatter(data['t'],data['BK150'],color='blue',label='150 GHz')
        scatter(data['t'],data['BK100'],color='green',label='100 GHz')
        #scatter(data['t'],data['BK270'],color='blue',label='270')
        ylabel("Tsky (K_rj)")
        xlabel("Date")
        title("Brightness Temperature for BK Bands")
        legend()
        grid()
        #savefig("plots/tsky/BKbands%s"%year)
        #show()
        clf()

        #print len(data['Tsky_tipper']), data['Tsky_tipper']
        #print len(data['averaged_Tsky_tipper']), data['averaged_Tsky_tipper']

        original = ~np.isnan(data['Tsky_tipper'])
        new = ~np.isnan(data['averaged_Tsky_tipper'])
        print len(data['averaged_Tsky_tipper'])
        print len(data['t'])
        print len(data['tipper850'])
        print len(data['t'][new])

        matplotlib.rcParams.update({'font.size': 18})

        if sum(new) ==0: #can't plot empty datasets
            continue
        if sum(original) == 0:
            continue 

        fig = figure(figsize = (35,15))
        plot(data['t_tipper'],data['Tsky_tipper'],'b.-',label='Tipper data')
        plot(data['t_tipper'],data['tipper850_interp'],'r.-',label='Resampled MERRA2 Predictions')
        plot(data['t'],data['tipper850'],'g.',label='Original 3-hour MERRA2 Predictions')
        legend()
        grid()
        ylabel("Tsky (K_rj)")
        xlabel("Date")
        title("Zenith Brightness Temperature (Tsky): %s %s%s"%(currentSite, year,month))
        savefig("merra2_products/plots/%s_timeseries_resampled_%s_groundData%s_%s"%(siteInitial,year,gndData,cloud))
        #show()
        close()

        #figure(figsize = (15,15))
        fig = figure(figsize = (20,20))
        tipperTsky = data['Tsky_tipper']
        merraTsky = data['tipper850_interp']
        general = (tipperTsky > 0.) & (merraTsky > 0.)
        tipperTsky = tipperTsky[general]
        merraTsky = merraTsky[general]
        fit = np.polyfit(merraTsky, tipperTsky, 1)
        unsaturated = (tipperTsky < 240.) & (merraTsky < 240.) #filter out saturated data (completely opaque atmosphere)
        unsaturatedFit = np.polyfit(merraTsky[unsaturated],tipperTsky[unsaturated],1)        fit_fn = np.poly1d(fit)
        bestfit = np.poly1d(unsaturatedFit)
        fit_perfect = np.poly1d([1,0])

        limits = [40,270]        scatter(merraTsky, tipperTsky, color = 'blue', label = 'data')        plot(limits, fit_fn(limits), color='red', label = "best fit for all data (R = %.2f): y = "%(np.corrcoef(merraTsky,tipperTsky)[1][0]) + str(fit_fn)[2:])
        plot(limits, bestfit(limits), color='green', label = "best fit for Tsky < 240K (R = %.2f): y = "%(np.corrcoef(merraTsky[unsaturated],tipperTsky[unsaturated])[1][0]) + str(bestfit)[2:])
        plot(limits, fit_perfect(limits), color = 'black', label = 'y = x')
        #errorbar(merraTsky,tipperTsky, yerr = data['Tsky_sigma'][new], fmt = 'o')
        ylabel("Tipper Tsky (K_rj)")
        xlabel("Resampled MERRA Tsky (K_rj)")
        ylim(limits)
        xlim(limits)
        legend(loc = 'upper left')
        grid()
        title("Comparison of Zenith Brightness Temperature (Tsky) at %s for %s%s"%(currentSite,year,month))
        print year, np.corrcoef(merraTsky,tipperTsky)[1][0], len(tipperTsky)
        savefig("merra2_products/plots/%s_correlation_resampled_%s_groundData%s_%s"%(siteInitial,year,gndData,cloud))
        #savefig("plots/tsky/Comparison%s%s%s"%(currentSite,year,month))
        #show()
        close()

        fig = figure(figsize = (35,15))
        scatter(data['t'],data['tipper850'], color = 'blue', label = 'MERRA2 Prediction')
        scatter(data['t'][new], data['averaged_Tsky_tipper'][new], color = 'red', label = '850GHz Tipper Data')
        #errorbar(data['t'][new], data['averaged_Tsky_tipper'][new], yerr = data['Tsky_sigma'][new], color = 'red', label = '850GHz Tipper Data')
        legend()
        grid()
        ylabel("Tsky (K_rj)")
        xlabel("Date")
        title("Zenith Brightness Temperature (Tsky): %s %s%s"%(currentSite, year,month))
        savefig("merra2_products/plots/%s_timeseries_averaged_%s_groundData%s_%s"%(siteInitial,year,gndData,cloud))
        #savefig("plots/tsky/JuxtaPose %s %s-%s"%(currentSite,year,month))
        #show()
        clf()
        close(fig)

        #figure(figsize = (15,15))
        fig = figure(figsize = (20,20))
        tipperTsky = data['averaged_Tsky_tipper'][new]
        merraTsky = data['tipper850'][new]
        if len(tipperTsky) ==0:
            continue
        fit = np.polyfit(merraTsky, tipperTsky, 1)
        unsaturated = (data['averaged_Tsky_tipper'] < 240.) & (data['tipper850'] < 240.) #filter out saturated data (completely opaque atmosphere)
        unsaturatedFit = np.polyfit(data['tipper850'][unsaturated],data['averaged_Tsky_tipper'][unsaturated],1)        fit_fn = np.poly1d(fit)
        bestfit = np.poly1d(unsaturatedFit)
        fit_perfect = np.poly1d([1,0])

        limits = [40,270]        scatter(merraTsky, tipperTsky, color = 'blue', label = 'data')        plot(limits, fit_fn(limits), color='red', label = "best fit for all data (R = %.2f): y = "%(np.corrcoef(merraTsky,tipperTsky)[1][0]) + str(fit_fn)[2:])
        plot(limits, bestfit(limits), color='green', label = "best fit for Tsky < 240K (R = %.2f): y = "%(np.corrcoef(data['tipper850'][unsaturated],data['averaged_Tsky_tipper'][unsaturated])[1][0]) + str(bestfit)[2:])
        plot(limits, fit_perfect(limits), color = 'black', label = 'y = x')
        #errorbar(merraTsky,tipperTsky, yerr = data['Tsky_sigma'][new], fmt = 'o')
        ylabel("Tipper Tsky (K_rj)")
        xlabel("MERRA Tsky (K_rj)")
        ylim(limits)
        xlim(limits)
        legend(loc = 'upper left')
        grid()
        title("Comparison of Zenith Brightness Temperature (Tsky) at %s for %s%s"%(currentSite,year,month))
        print year, np.corrcoef(merraTsky,tipperTsky)[1][0], len(tipperTsky)
        savefig("merra2_products/plots/%s_correlation_averaged_%s_groundData%s_%s"%(siteInitial,year,gndData,cloud))
        #savefig("plots/tsky/Comparison%s%s%s"%(currentSite,year,month))
        #show()
        clf()
        close()
        close()
        close(fig)

        fig = figure(figsize = (25,10))
        residual = merraTsky - tipperTsky
        residual = residual - mean(residual)
        chisquared = np.sum((residual/data['Tsky_sigma'][new])**2)
        errorbar(data['t'][new],residual,yerr=data['Tsky_sigma'][new],label='Averaged, mean-subtracted residual (%s points): '%(len(residual)) + r'$\chi^2=%s$'%chisquared)
        ylabel("Residual Tsky [MERRA - Tipper] (K)")
        xlabel("Date")
        grid()
        legend()
        title("Residual of Tsky (MERRA > Tipper) at %s for %s%s"%(currentSite,year,month))
        savefig("merra2_products/plots/%s_residual_averaged_%s_groundData%s_%s"%(siteInitial,year,gndData,cloud))
        #savefig("plots/tsky/Residual%s%s%s"%(currentSite,year,month))
        #show()
        clf()
        close(fig)

        fig = figure(figsize = (25,10))
        residual = data['tipper850_interp']-data['Tsky_tipper']
        residual = residual - mean(residual)
        chisquared = np.sum((residual/data['tipper_unc'])**2)
        errorbar(data['t_tipper'],residual,yerr=data['tipper_unc'],label = 'Resampled, mean-subtracted residual (%s points): '%(len(residual))+r'$\chi^2=%s$'%chisquared)
        grid()
        legend()
        ylabel("Residual Tsky [MERRA - Tipper] (K)")
        xlabel("Date")
        title("Mean-Subtracted Residual of Tsky (MERRA > Tipper) at %s for %s%s"%(currentSite,year,month))
        savefig("merra2_products/plots/%s_residual_resampled_%s_groundData%s_%s"%(siteInitial,year,gndData,cloud))
        #show()
        close()

        fig = figure(figsize = (25,10))
        fractionalresidual = 2*(merraTsky-tipperTsky)/(merraTsky+tipperTsky)
        scatter(data['t'][new],fractionalresidual)
        ylabel("Fractional Residual")
        xlabel("Date")
        grid()
        title("Fractional Residual of Tsky at %s for %s%s"%(currentSite,year,month))
        savefig("merra2_products/plots/%s_fractional_averaged_%s_groundData%s_%s"%(siteInitial,year,gndData,cloud))
        #savefig("plots/tsky/FracResidual%s%s%s"%(currentSite,year,month))
        #show()
        clf()
        close(fig)

        fig = figure(figsize = (25,10))
        fractionalresidual = 2*(merraTsky-tipperTsky)/(merraTsky+tipperTsky)
        scatter(data['t_tipper'],data['tipper850_interp']-data['Tsky_tipper'])
        grid()
        ylabel("Fractional Residual")
        xlabel("Date")
        title("Fractional Residual of Tsky (MERRA > Tipper) at %s for %s%s"%(currentSite,year,month))
        savefig("merra2_products/plots/%s_fractional_resampled_%s_groundData%s_%s"%(siteInitial,year,gndData,cloud))
        #show()
        close()

        fig = figure(figsize = (30,15))
        tauFilter = data['tau'] > 0.
        scatter(data['t_tipper'][tauFilter], data['tau'][tauFilter], color = 'red', label = 'Tipper opacity')
        scatter(data['t'],data['tau850'],color='blue',label = 'MERRA2-predicted opacity')
        legend()
        grid()
        title("Optical Depth (Tau) at %s for %s%s"%(currentSite,year,month))
        savefig("merra2_products/plots/%s_tau_averaged_%s_groundData%s_%s"%(siteInitial,year,gndData,cloud))
        #show()
        clf()
        close(fig)

        fig = figure(figsize = (15,15))
        tau_m, tau_s = m.averageTipper(data['t'], data['t_tipper'], data['tau'])
        tau_m = np.array(tau_m)
        generalFilter = tau_m < 8. #get rid of clear outliers
        smallTauFilter = (tau_m < 4.) & (data['tau850'] < 4.)  #only do small tau. The fit should be better in this case, so it refers to bestfit variable. 
        allTau =  data['tau850'][generalFilter]
        smallTau = data['tau850'][smallTauFilter]
        scatter(allTau,tau_m[generalFilter])
        generalfit = np.polyfit(allTau, tau_m[generalFilter], 1)
        bestfit = np.polyfit(smallTau, tau_m[smallTauFilter], 1)        taufit = np.poly1d(generalfit)
        BestTauFit = np.poly1d(bestfit)
        print year, np.corrcoef(merraTsky,tipperTsky)[1][0], len(tipperTsky)
        plot(allTau, taufit(allTau), color='red', label = "best fit for all data (R = %.2f): y = "%(np.corrcoef(allTau,tau_m[generalFilter])[1][0]) + str(taufit)[2:])
        plot(allTau, BestTauFit(allTau), color='green', label = "best fit for tau < 4 (R = %.2f): y = "%(np.corrcoef(smallTau,tau_m[smallTauFilter])[1][0]) + str(BestTauFit)[2:])
        plot(allTau, fit_perfect(allTau), color='black', label = 'y = x')
        legend(loc = 'upper left')
        xlabel("MERRA2 Tau")
        ylabel("Tipper Tau")
        grid()
        title("Comparison of Tau at %s for %s%s"%(currentSite,year,month))
        savefig("merra2_products/plots/%s_taucorrelation_averaged_%s_groundData%s_%s"%(siteInitial,year,gndData,cloud))
        #show()
        clf()
        close(fig)

        fig = figure(figsize = (25,10))
        tauresidual = allTau - tau_m[generalFilter]
        scatter(data['t'][generalFilter],tauresidual)
        ylabel("Residual tau [MERRA - Tipper]")
        xlabel("Date")
        grid()
        title("Residual of tau (MERRA > Tipper) at %s for %s%s"%(currentSite,year,month))
        savefig("merra2_products/plots/%s_tauresidual_averaged_%s_groundData%s_%s"%(siteInitial,year,gndData,cloud))
        #savefig("plots/tsky/Residual%s%s%s"%(currentSite,year,month))
        #show()
        clf()
        close(fig)

        fig = figure(figsize = (25,10))
        fractionalresidual = 2*(allTau-tau_m[generalFilter])/(allTau+tau_m[generalFilter])
        scatter(data['t'][generalFilter],fractionalresidual)
        ylabel("Fractional Residual")
        xlabel("Date")
        grid()
        title("Fractional Residual of tau at %s for %s%s"%(currentSite,year,month))
        savefig("merra2_products/plots/%s_taufractional_averaged_%s_groundData%s_%s"%(siteInitial,year,gndData,cloud))
        #savefig("plots/tsky/FracResidual%s%s%s"%(currentSite,year,month))
        #show()
        clf()
        close(fig)

        #Plots that are currently not used
        #figure(figsize = (15,15))
        freqResidual = np.fft.fft(residual)
        n = len(residual)
        k = np.arange(n)
        Fs = 1./(60*180)
        T = n/Fs
        frq = k/T
        frq = frq[range(n/2)]
        freqResidual = freqResidual[range(n/2)]
        #plot(frq,abs(freqResidual))
        #show()
        #clf()

        #figure(figsize = (25,15))
        #scatter(data['t'],data['pwv'])
        #title("PWV")
        #show()
        

#plot(f850,band850rj, color = 'black', label = 'bandpass filter')
#plot(fs, trj/mean(trj), color = 'red', label = 'Rayleigh-Jeans Tsky')
#plot(fs, tb/mean(tb), color = 'blue', label = 'Brightness Temp')
#plot(fs, tau/mean(tau), color = 'green', label = 'Opacity')
#legend()
#show()
