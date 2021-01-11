import os
import urllib
import datetime
import dateutil.rrule as rr
import dateutil.parser as dparser
#import netCDF4 # TODO figure out how to make netCDF work
import subprocess
from pylab import *
import glob
import http.cookiejar as cookielib
import traceback
import time
import base64
from scipy import interpolate, stats, signal

# WARNING
# Needs to have activated a special environment
# source activate mypyenv (which has scipy 0.18.1)
# https://www.rc.fas.harvard.edu/resources/documentation/software-on-odyssey/python/
# add some more warnings here.

class merra2Player():
    """
    class to download merra2 data and generate am profiles from it
    """
    def __init__(self):
        """

        TODO: reduce the time to run am even further 
        by checking on the changes of the integrated band

        """
        self.version = 1.0
        self.debug= True
        self.verbose = True
        self.defineSite()
        self.product = 'inst' #averaged (tavg) or instantaneous (inst)
        self.checkLinks(product=self.product)
       
    def checkLinks(self, product = 'inst'):

        if not(os.path.islink('merra2_products')):
            print("WARNING: You must have a symlink which points to \"merra2_products\"")
            print("WARNING: ln -s /n/holylfs02/LABS/kovac_lab/keck/wvr_products/merra2_analysis merra2_products")
            print("Stopping")
            exit()

        if product == 'inst':
            folder = 'merra2_products'
            self.webDir = '/n/holylfs02/LABS/kovac_lab/www/merra2_web/web_output_files/'
        elif product == 'tavg':
            folder = 'merra2_products_averaged'
            self.webDir = 'merra2_products_averaged/web_output_files/'
        self.dataDir = '%s/%s/'%(folder,self.site['name'])
        if not os.path.exists(self.dataDir):
            os.system('mkdir %s/%s'%(folder,self.site['name']))
        self.merraDir = self.dataDir+'merra2_raw_data/'
        if not os.path.exists(self.merraDir):
            os.system('mkdir %s/%s/merra2_raw_data/'%(folder,self.site['name']))
        self.tipperDir = self.dataDir+'tipper_raw_data/'
        if not os.path.exists(self.tipperDir):
            os.system('mkdir %s/%s/tipper_raw_data/'%(folder,self.site['name']))
        self.amcDir = self.dataDir+'amcFiles/'
        if not os.path.exists(self.amcDir):
            os.system('mkdir %s/%s/amcFiles/'%(folder,self.site['name']))
        self.finalDir = self.dataDir+'final/'
        if not os.path.exists(self.finalDir):
            os.system('mkdir %s/%s/final/'%(folder,self.site['name']))
        
        if not(os.path.islink('am')):
            print("WARNING: You must have a symlink which points to ''am'' execuatble to run this")
            print("Stopping")
            exit()

        self.webDir = '/n/holylfs02/LABS/kovac_lab/www/merra2_web/web_output_files/'
        self.wxDir = '%s/wx_data/'%folder #independent of holylfs directory
        self.auxDataDir = '/n/holylfs02/LABS/kovac_lab/keck/keck_aux_data/bandpass/'

    def defineSite(self, site = None):
        """
        function to define a site dictionnary
        
        """
        if site == None or site == {}:
            site ={}
            site['type']='presel'
            site['name']='SouthPole'
            site['gndData']= 2
            site['cldtype'] = 'vapor'
            
        if site['type']=='presel':
            if site['name']== 'SouthPole':
                site['sname']='sp'
                site['lat'] = -90.0
                site['lon'] = 0.0
                site['alt'] = 2835
            elif site['name'] =='ChajnantorPlateau':
                site['sname']='cp'
                site['lat'] =  -23.0285
                site['lon'] = -67.76175
                site['alt'] =  5060
            elif site['name'] == 'ChajnantorCerro':
                site['sname']='cc'
                site['lat'] =  -22.9856
                site['lon'] = -67.74027
                site['alt'] =  5612
            elif site['name'] == 'MaunaKea':
                site['sname']='mk'
                site['lat'] =  19.8
                site['lon'] = -155.45 
                site['alt'] =  4100
            elif site['name'] == 'Summit':
                site['sname']='su'
                site['lat'] = 72.5785
                site['lon'] = -38.4525 
                site['alt'] =  3200
            elif site['name'] == 'Qubic':
                site['sname']='qu'
                site['lat'] = -24.18628
                site['lon'] = -66.47817 
                site['alt'] =  4870    
            else:
                print("Only the following default sites are defined for the moment:",)
                print("SouthPole, ChajnantorPlateau, ChajnantorCerro, MaunaKea, Summit")
                print("Use default site above or define site dict with keys:'type: custom, name, sname, lat, long, alt'")
                raise ValueError('siteName error')
                
        elif site['type'] =='custom':
            ## check that name, sname, lat,long,alt are defined and within bounds
            keys = ['name','sname','lat','lon','alt']
            for k in keys:
                if k in site.keys():
                    pass
                else:
                    print("%s key missing in input custom site dictionary"%k)
                    raise ValueError('siteName error') 
               
        if 'gndData' not in site.keys():
            site['gndData']= 2
        if 'cldtype' not in site.keys(): # TODO we probably want to look at cldtype to get the ice information out
            site['cldtype'] = 'vapor'

        self.site = site
        return site


    def runMERRA(self, dateopt={'start':'20150123','end':'20150123'}, bandopt={'name':'BK95'}):
        """
        second level function
        start/end are ascii 'YYYYMMDD'
        siteopt options are already defined in self.site
        bandopt: by default, compute all BK and tipper bands. If bandopt['name']='custom'
        compute also that custom band
        """
        dateList = self.retrieve_merra2_data_for_dateRange(dateopt['start'],dateopt['end']) 
        if self.site['gndData'] == 1:
            self.retrieve_SPWx_data_for_dateRange(dateopt['start'],dateopt['end'])

        keys = ['t','pwv','BK100','BK150','BK210','BK220','BK270','tipper850','tau850','EHT_hi','EHT_lo']
        tsky={}
        for k in keys: tsky[k]=[]
        f270,band270f,band270rj = self.readBandpass(bandopt={'name':'BK270'})
        f220,band220f,band220rj = self.readBandpass(bandopt={'name':'BK220'})
        f210,band210f,band210rj = self.readBandpass(bandopt={'name':'BK210'})
        f150,band150f,band150rj = self.readBandpass(bandopt={'name':'BK150'})
        f100,band100f,band100rj = self.readBandpass(bandopt={'name':'BK95'})
        f850,band850f,band850rj = self.readBandpass(bandopt={'name':'tipper850'})
        feht,bandeht_hi,bandeht_hi = self.readBandpass(bandopt={'name':'EHT_hi'})
        feht,bandeht_lo,bandeht_lo = self.readBandpass(bandopt={'name':'EHT_lo'})

        if bandopt['name'] == 'custom':
             fc,bandcf,bandcrj = self.readBandpass('custom')
             tsky['custom']=[]

        for date in dateList:
            amcFileList = self.checkAmcFileForDate(date)
            if size(amcFileList) != 0:
                print("Found %d amc input profiles for %s. Not regenerating ... " %(size(amcFileList),date))
            else:
                profile = self.createProfile(date.strftime('%Y%m%d'), plotFig=False)
                amcFileList = self.profile2am(profile)

            for amcFile in amcFileList:
                dat = self.getDatetimeFromAmcFile(amcFile)
                tsky['t'].append(dat)
                try:
                    fs,tb,trj,pwv,tau = self.run_am(amcFile, f0 = 0.0, f1 = 1200.0, df = 1000.0)
                    tsky['pwv'].append(pwv)
                    tsky['BK100'].append(self.integBand(f100,band100rj,fs,trj))
                    tsky['BK150'].append(self.integBand(f150,band150rj,fs,trj))
                    tsky['BK210'].append(self.integBand(f210,band210rj,fs,trj))
                    tsky['BK220'].append(self.integBand(f220,band220rj,fs,trj))
                    tsky['BK270'].append(self.integBand(f270,band270rj,fs,trj))
                    tsky['tipper850'].append(self.integBand(f850,band850rj,fs,tb)) #Note: Integrate tb because its a volume-based Planck law rather than a 1d. 
                    tsky['EHT_hi'].append(self.integBand(feht,bandeht_hi,fs,trj))
                    tsky['EHT_lo'].append(self.integBand(feht,bandeht_lo,fs,trj))
                    if band['name'] == 'custom':
                        tsky['custom'].append(self.integBand(fc,bandcrj,fs,trj))
                    integratedTx = self.integBand(f850,band850rj,fs,exp(-tau))
                    tsky['tau850'].append(-log(integratedTx))
                except:
                    pass
            
        datestr_list = [d.strftime('%Y-%m-%d T%H:%M:%S') for d in tsky['t']]
        tsky['tstr']=datestr_list

        return  tsky

    def runTipper(self, tsky_merra, dateopt={'start':'20150123','end':'20150123'}):
        """
        function to process tipper data
        tsky_merra is existing dictionnary returned by runMERRA
        start/end are ascii 'YYYYMMDD'
        """
        try:
            t_merra = tsky_merra['t']   # merra2 time stamps
            (t_tipper, tipperData) = self.readTipper(dateopt['start'],dateopt['end'])  # raw tipper data (tau, Tatm) and time stamps
            Tsky_tipper, Tsky_tipper_sigma = self.calc_Tsky_sigma(tipperData[7],tipperData[8],tipperData[9],tipperData[10])  # measured Zenith tipper Tsky 
            Tsky_tipper_av, Tsky_tipper_av_sigma = self.averageTipper(t_merra, t_tipper, Tsky_tipper, Tsky_tipper_sigma,product=self.product)   # measuremd zenith tipper Tsky averaged over merra2 timestamp
            # commented out 22 Aug 2017 by dB. this needs source activate mypyephem. Will fix later.
            # tipper850_interp =  self.resample(t_merra,t_tipper,tsky_merra['tipper850']) 
            # tsky_merra['tipper850_interp'] = tipper850_interp #merra2 Tsky resampled to the tipper time.
            tsky_merra['t_tipper'] = t_tipper  # raw tipper time stamps
            tsky_merra['Tsky_tipper'] = Tsky_tipper  #  raw tipper Tsky at raw time stamps
            tsky_merra['tipper_unc'] = Tsky_tipper_sigma # tipper Tsky_sigma at raw time axis
            tsky_merra['averaged_Tsky_tipper'] = Tsky_tipper_av # tipper Tsky at merra2 time axis
            tsky_merra['Tsky_sigma'] = Tsky_tipper_av_sigma   # tipper Tsky_sigma at merra2 time axis

            tsky_merra['tau'] = tipperData[7]
            tsky_merra['Tatm'] = tipperData[9]
        except:
            print("runTipper exception")
            # if there is no tipper data.
            keys = ['t_tipper','averaged_Tsky_tipper','Tsky_sigma','Tsky_tipper','tau','Tatm','tipper850_interp','tipper_unc']
            for k in keys:
                #a hack so that the arrays have reasonable lengths.
                tsky_merra[k] = np.empty(len(tsky_merra['t'])) * np.nan 

        return tsky_merra

    def readTipper(self, start='', end=''):
        """
        Reads in the raw tipper data
        input: a start and end date of format 'YYYYMMDD'
        output: raw tipper time and tipper data
        """
        print("Reading tipper data ...")
        #conv = {4: lambda s: epoch+datetime.timedelt(days=s)}

        #get a list of all the months for which we want tipper data
        ds = datetime.datetime.strptime(start,'%Y%m%d')
        if end == '':
            de = ds
        else:
            de = datetime.datetime.strptime(end,'%Y%m%d')
        months = [dt for dt in rrule(MONTHLY, dtstart=ds, until=de)]
        tipperTime = []
        tipperData = []

        initial = self.site['sname']
        for month in months:
            yearMonth = month.strftime('%Y-%m')
            tipperFile = self.tipperDir+'%s-%s-0.dat'%(initial,yearMonth)
            try:
                tipperData.append(genfromtxt(tipperFile))
            except:
                pass
        tipperData = concatenate(tipperData, axis = 0)
        tipperData = np.transpose(tipperData)
        
        #create datetime
        epoch = datetime.datetime(1995, 1, 1, 0, 0, 0)
        for f in tipperData[4,:]:
            tipperTime.append(epoch+datetime.timedelta(days=f))
        tipperTime = array(tipperTime)

        # chop tipper tod to request start and end dates
        #istart,val0 = self.index_nearest(tipperTime,ds , 1)
        #iend, val1 = self.index_nearest(tipperTime,de , 1)

        #tipperTime = tipperTime[istart[0]:iend[0]]
        #tipperData = tipperData[:,istart[0]:iend[0]]

        return (tipperTime,tipperData)

    def calc_Tsky_sigma(self,tau,sigma_tau,tatm,sigma_tatm,theta=0.8,sigma_theta=0.03):
        """
        
        calculates the propagated uncertainty onto Tsky 
        given tau, tatm, Window Tx  theta, and stdev on each.
        input can be single values or arrays. outputs will be in the same format.
        """
        print("Calculating tipper Tsky and Tsky_error")
        tsky = (tatm/theta) * (1-exp(-tau))

        dtatm = (1-exp(-tau))/theta*sigma_tatm
        dtheta = tatm*(1-exp(-tau))/theta**2 * sigma_theta
        dtau = tatm*exp(-tau)/theta * sigma_tau
        dtsky = sqrt(dtatm**2+dtheta**2+dtau**2)

        return tsky, dtsky

    def resample(self, merraTime, tipperTime, merraTsky):
        """
        Convert the 3-hourly MERRA2 prediction timestream into a  timestream samples
        a t_tipper
        input: merraTime is the 3-hour timestamps of MERRA. 
        tipperTime is the irregular timestamps of the tipper data.
        merraTsky is the 3-hour Tsky predictions of MERRA2.
        output: An array with Tsky['tipper850'] resampled to the times of tipperTime.
        """
        #print "Resampling Merra2 data onto tipper time"

        # helper function to convert datetime array into seconds array
        datetime2seconds = vectorize(lambda x: x.total_seconds())

        minutes = 1
        # merra2 is exactly every 3 hours = 365*24/3. 2920 pts per year
        # 365 * 24 * 60 / 1  = 525600 points per year at 1mn interval
        Nhi = len(merraTsky)*3*60 / minutes

        # create regular grid of time at 1mn interval
        t = list(rr.rrule(rr.MINUTELY,interval = minutes, dtstart=merraTime[0], count = Nhi))
        dt = array(t) - merraTime[0]
        dt = datetime2seconds(dt)

        #convert merra2 time into array of seconds since start of merra2 year.
        dt_merra2 = array(merraTime) - merraTime[0]
        dt_merra2 = datetime2seconds(dt_merra2)
        # resample merra2 data to new 10mn time interval
        tipper850_resamp = signal.resample(merraTsky,Nhi)

        # convert tipper time into array of seconds since start of merra2 year.
        dt_tipper = array(tipperTime) - merraTime[0]
        dt_tipper = datetime2seconds(dt_tipper)

        f  = interpolate.interp1d(dt,tipper850_resamp, fill_value = 'extrapolate')
        tipper850_interp = f(dt_tipper)
        
        return tipper850_interp

    def averageTipper(self, newdt, dt, Tsky, unc, product = 'inst', plotFig = False):
        """
        input: dt, Tsky, unc from readTipper. product defines the way the tipper is interpolated to 3-hour intervals. inst takes the nearest tipper point and its uncertainty, 
        while aver takes the weighted averaged of all tipper points within the 60 minute range, and a weighted error on the mean. inst is default
        newdt is a datetime array of merra2 times (3hour interval)
        tipper data is collected every 12.75 minutes, merra2 every3 hours so we average ~14 tipper points to get 1.
        option to plot the original data alongside the averaged data is available.
        
        output: 3-hour tipper averages and tipper uncertainty measurements.
        """
        #print "Averaging tipper onto MERRA2 time"
        tipper_m = []
        tipper_s = []

        if product == 'inst':
            for i in newdt:
                #for every data point of the 3 hour MERRA2 data, find the nearest datetime in the Chajnantor data.
                #average the data from the 3 hour period centered around it
                #try-except is for edge cases when there isn't data on one side.
                #this is more or less wrong, because we don't always have data every 12.75 minutes.
                closest,value = self.index_nearest(dt, i, 1)
                try:
                    closest = closest[0]
                except:
                    print(value)
                    tipper_m.append(float("nan"))
                    tipper_s.append(float("nan"))
                    continue

                #deal with gaps in tipper data by setting them as NaN if the distance between nearest point is too far.
                interval = abs(dt[closest]-i)
                if interval > datetime.timedelta(minutes = 20):
                    tipper_m.append(float("nan"))
                    tipper_s.append(float("nan"))
                else:
                    tipper_m.append(Tsky[closest]) #take the closest Tsky and uncertainty, for a time that must be within 20 minutes.
                    tipper_s.append(unc[closest])
        
        elif product == 'tavg':
            for i in newdt:
                #for every data point of the 3 hour MERRA2 data, find the nearest datetime in the Chajnantor data.
                #average the data from the 2 hour period centered around it
                closePoints = (dt > i-datetime.timedelta(hours=1)) & (dt < i+datetime.timedelta(hours=1))

                #deal with gaps in tipper data by setting them as NaN if the distance between nearest point is too far.
                if sum(closePoints)==0:
                    tipper_m.append(float("nan"))
                    tipper_s.append(float("nan"))
                else:
                    weights = (1/unc[closePoints])**2
                    average = np.average(Tsky[closePoints],weights=weights)
                    #variance = 1/np.average()
                    tipper_m.append(average) #weighted average of Tsky. ##TODO: make it so that you actually average the 60 minutes within. 
                    tipper_s.append(0) #weighted error on the mean. 
        else:
            raise ValueError("invalid product name. choose inst or aver.")

        if (plotFig):
            allData = scatter(dt, Tsky, label = "all data, every 12.75 minutes")
            averagedData = errorbar(newdt, tipper_m, yerr = tipper_s, elinewidth = 2, color = "red", fmt = "o", label = "averaged data, every 3 hours")
            legend(handles = [allData, averagedData])

            grid()
            xlim([newdt[0], newdt[-1]])
            xlabel('DateTime')
            ylabel('Tsky (K)')
            suptitle("Tipper data: raw and averaged")
            savefig("test.png")

        return tipper_m, tipper_s

    def getDatetimeFromAmcFile(self,fn):
           # get the datetime for this file
            date = fn.split('_')[1]
            time = fn.split('_')[2]
            return datetime.datetime.strptime('%s %s'%(date,time),'%Y%m%d %H%M%S')

    def checkAmcFileForDate(self,date):

        gnd = self.site['gndData']
        cldtype = self.site['cldtype']
        if type(date) is str:
            datestr = date
        elif type(date) is datetime.datetime:
            datestr = date.strftime('%Y%m%d')
        cwd = os.getcwd()
        os.chdir(self.amcDir)
        amcFiles = glob.glob('%s_%s_*_gndData%s_%s*.amc'%(self.site['name'],datestr,gnd,cldtype))
        os.chdir(cwd)

        return sort(amcFiles)
    
    def retrieve_SPWx_data_for_dateRange(self, dstart='', dend=''):
        
        ds = datetime.datetime.strptime(dstart,'%Y%m%d')
        de  = datetime.datetime.strptime(dend,'%Y%m%d')
        for month in rrule(MONTHLY, dtstart=ds, until=de):
            date = month.strftime('%Y%m%d')
            ret = self.retrieve_SPWx_data_for_date(date)

        return ret

    def retrieve_SPWx_data_for_date(self, date = ''):
        """
        New Function by Andrew
        Download from:
        # ftp://aftp.cmdl.noaa.gov/data/meteorology/in-situ/spo/[year]/met_spo_insitu_1_obop_minute_[year]_[month].txt;
        for a given date/time string time should be in 1 hour increment
        """
        year = date[0:4]
        month = date[4:6]

        filename = 'met_spo_insitu_1_obop_minute_%s_%s.txt'%(year,month)
        if os.path.exists(self.wxDir+filename):
            print(self.wxDir+filename+" already exists, skipping... ")
            return 1
        url = "ftp://aftp.cmdl.noaa.gov/data/meteorology/in-situ/spo/%s/met_spo_insitu_1_obop_minute_%s_%s.txt"%(year,year,month)
        print("Downloading data from %s"%url)
        try:
            f = urllib2.urlopen(url)
            data = f.read()
            with open(self.wxDir+filename, "wb") as code:
                code.write(data)
            return 1
        except:
            ValueError('Download failed ...')
            return 0

    def readSPWxData(self, date='', time='',shrink = 0, verbose=True,):
        """
        Instead of reading B2/Keck Wx data, reads minutely spo met data.
        Downloaded from ftp://aftp.cmdl.noaa.gov
        for a given date/time string time should be in 1 hour increment
        shrink [hrs] is number of hours before and after we want to select. 
        if zero provides the whole month file

        """
        dat=datetime.datetime.strptime('%s %s'%(date,time),'%Y%m%d %H%M%S')
        print(dat)
        year = date[0:4]
        month = date[4:6]
        wxFile = 'met_spo_insitu_1_obop_minute_%s_%s.txt'%(year,month)

        if verbose: print("Reading %s"%wxFile)
        d = genfromtxt(self.wxDir+wxFile,delimiter='',dtype='S,i,i,i,i,i,i,f,i,f,f,f,f,f,f')
        dt = array([datetime.datetime(t[1],t[2],t[3],t[4],t[5]) for t in d])
        wx ={'time': dt, 'presmB':d['f9'],'tempC': d['f11'],'rh':d['f13'],'wsms':d['f7'],'wddeg':d['f6']}

        # exclude junk data
        q = flatnonzero((wx['presmB'] != -999.9) & (wx['tempC'] != -999.9) & (wx['rh'] != -99.))
        if q != []:
            for k in wx.keys():
                wx[k]=wx[k][q]
        # down-select to given time range
        if shrink != 0:
            f = shrink*60
            nrows = size(wx['time'])
            q = flatnonzero(wx['time'] == dat)
            while(size(q) == 0):
                dat=dat+datetime.timedelta(seconds=60)
                q = flatnonzero(wx['time'] == dat)
            for k in wx.keys():
                wx[k]=wx[k][max(q-60,0):min(q+60,nrows)]

        return wx

    def readWxData(self, date='', time='', filename='', verbose=True):
        """
        reads the Keck/B2 Wx data given a
        date, time.
        Gives generic pressure if Wxdata directory not available.
        date format: 'YYYYMMDD'
        time format: 'HHMMSS'
        time should be increment of 1 hour, ie 000000, 010000
        """

        if not(os.path.exists(self.wxDir)) or  not(os.path.exists(self.wxDir+filename)):
            print("WARNING: Keck/B2 Wx data not available, assigning default pressure...")
            if self.site['name'] == 'SouthPole':
                return {'presmB':[680]}

        if filename != '':
            wxFile = filename
        elif date != '':
            cwd = os.getcwd()
            os.chdir(self.wxDir)
            wxFile = glob.glob('*%s_%s*'%(date,time))[0]
            os.chdir(cwd)
        else:
            print("specify either filename to be read or date/time strings")
            return

        if verbose: print("Reading %s"%wxFile)
        wx = genfromtxt(self.wxDir+wxFile, delimiter='',names=True, invalid_raise = False)

        return wx

    def makeWxMeans(self):
        """
        to be moved to weather Lib
        Not really used for not for anything in here.
        """
        start=datetime.datetime.strptime('2010','%Y')
        dateList = list(rr.rrule(rr.MONTHLY,start, count=12))

        for month in dateList:
            print(month)
            cwd = os.getcwd()
            os.chdir(self.wxDir)
            dat = month.strftime('%Y%m')+'*'
            files = glob.glob('%s_wx_B2.txt'%dat)
            print(size(files))
            os.chdir(cwd)
            mwx = []
            for f in files:
                wx = self.readWxData(filename = f, verbose=False)
                mwx.append(mean(wx['presmB']))
            clf()
            plot(mwx)
            print(dat, mean(mwx), std(mwx))
            raw_input()

    def closestPoints(self, lat, lon):
        """
        latitude ranges from -90 to 90 (inclusive) 
        longitude ranges from -180 to 180 (inclusive)
        This retuns the lat/long grid MERRA2 grid points
        """
        lat0, lat1, lon0, lon1 = 0,0,0,0
        for i in arange(-90, 90, 0.5):
            if (lat >= i and lat <= i + 0.5):
                lat0 = i
                lat1 = i+0.5
                break
        for j in arange(-180, 180, 0.625):
            if (lon >= j and lon <= j + 0.625):
                lon0 = j
                lon1 = j + 0.625
                break
        return lat0, lat1, lon0, lon1

    def get_bbox(self, lat, lon):
        """
        define a bounding box to search the grid in.
        Default Bounding Box  is 1deg in lat by 1.25deg in long

        """
        lat0, lat1, lon0, lon1 = self.closestPoints(lat, lon)
        lat1 += 0.5
        lon1 += 0.625
        return [min(lat0,lat1),max(lat0,lat1),lon0,lon1]

    def get_host_name(self,dataset):
        """
        Returns  the  hostname from which to retrieve the MERRA2 data
        """
        if dataset == 'multiLevel':
            hostname ='http://goldsmr5.gesdisc.eosdis.nasa.gov/daac-bin/OTF/HTTP_services.cgi?'
        elif dataset == 'singleLevel':
            hostname = 'http://goldsmr4.gesdisc.eosdis.nasa.gov/daac-bin/OTF/HTTP_services.cgi?'
        else:
            raise ValueError("dataset must be 'singleLevel' or 'sultiLevel'. You entered %s"%dataset)
        return  hostname

    def get_url_for_date(self, date, dataset = "multiLevel"):
        """
        Returns the URL of the MERRA2 website for a given date
        Date: (string) YYYYMMDD
        dataset: "multiLevel" for 3D Merra product, "singleLevel" for 2D
        """

        [lat0,lat1,long0,long1] = self.get_bbox(self.site['lat'], self.site['lon'])

        year = date[:4]
        month = date[4:6]
        dt = datetime.datetime.strptime(date,'%Y%m%d')
        if dt >= datetime.datetime.strptime('20110101','%Y%m%d'):
            streamN = '400'
        elif dt >= datetime.datetime.strptime('20010101','%Y%m%d'):
            streamN = '300'
        elif dt >= datetime.datetime.strptime('19920101','%Y%m%d'):
            streamN = '200'
        elif dt >= datetime.datetime.strptime('19800101','%Y%m%d'):
            streamN = '100'

        if dataset == "multiLevel":
            if self.product == 'inst':
                shortname = 'M2I3NPASM'
                filenameBase = 'FILENAME=/data/s4pa/MERRA2/M2I3NPASM.5.12.4'
                filename = 'MERRA2_%s.inst3_3d_asm_Np.%s.nc4'%(streamN,date)
            elif self.product == 'tavg':
                shortname = 'M2T3NVASM'
                filenameBase = 'FILENAME=/data/s4pa/MERRA2/M2T3NVASM.5.12.4'
                filename = 'MERRA2_%s.tavg3_3d_asm_Nv.%s.nc4'%(streamN,date)
            format = '&FORMAT=bmM0Lw&BBOX=%2.1f,%2.3f,%2.1f,%2.3f&LABEL=svc_%s&FLAGS=&SHORTNAME=%s&SERVICE=SUBSET_MERRA2&LAYERS=&VERSION=1.02&VARIABLES=t,ps,qv,o3,h,ql,qi'%(lat0,long0,lat1,long1,filename,shortname)
            hostname = self.get_host_name(dataset)

        elif dataset == "singleLevel":
            if self.product == 'inst':
                shortname = 'M2I1NXASM'
                filenameBase = 'FILENAME=/data/MERRA2/M2I1NXASM.5.12.4'
                filename = 'MERRA2_%s.inst1_2d_asm_Nx.%s.nc4'%(streamN,date)
            elif self.product == 'tavg':
                shortname = 'M2T1NXSLV'
                filenameBase = 'FILENAME=/data/MERRA2/M2T1NXSLV.5.12.4'
                filename = 'MERRA2_%s.tavg1_2d_slv_Nx.%s.nc4'%(streamN,date)
            format = '&FORMAT=bmM0Lw&BBOX=%2.1f,%2.3f,%2.1f,%2.3f&LABEL=%s&FLAGS=&SHORTNAME=%s&SERVICE=SUBSET_MERRA2&LAYERS=&VERSION=1.02&VARIABLES=ps,qv2m,t2m,tqi,tql,tqv'%(lat0,long0,lat1,long1,filename,shortname)
            hostname = self.get_host_name(dataset)
        else:
            raise ValueError("dataset must be 'SingleLevel' or 'MultiLevel'. You entered %s"%dataset)

        url= '%s%s/%s/%s/%s&%s'%(hostname,filenameBase,year,month,filename,format)
        return url, filename

    def retrieve_merra2_data_for_date(self, date, dataset = 'multiLevel'):
        """
        date format: '20160301'
        if data is saved and present, returns 1
        if data is not present returns 0
        """
        print("Retrieving data for %s, %s "%(date, dataset))
        url, filename = self.get_url_for_date(date, dataset)
        if os.path.exists(self.merraDir+filename):
            if os.stat(self.merraDir+filename).st_size < 10000 :
                print(self.merraDir+filename+"  was corrupted at download, re-download ... ")
            else:
                print(self.merraDir+filename+"  already exists, skipping...  ")
                return 1
    
        try:
            if self.verbose: print("downloading MERRA2 URL: \n %s"%url)
            #https://docs.python.org/2/howto/urllib2.html#id6 documentation for urllib2 authorization
            #http://stackoverflow.com/questions/2407126/python-urllib2-basic-auth-problem for a fix
            # create a password manager
            username = 'anwang16' #TODO change this, probably
            password = 'AstroCMB1'
            request = urllib2.Request(url)
            base64string = base64.b64encode('%s:%s' % (username, password))

            cj = cookielib.CookieJar()
            opener = urllib2.build_opener(urllib2.HTTPCookieProcessor(cj))
            request.add_header("Authorization", "Basic %s" % base64string)
            furl = opener.open(request)

            meta = furl.info()
            file_size = int(meta.getheaders("Content-Length")[0])
            print("Downloading: %s Bytes: %s" % (filename, file_size))
            with open(self.merraDir+filename,'wb') as output:
                output.write(furl.read())
            print("Saved file: %s"%(self.merraDir+filename))
            return 1
        except Exception:
            print('generic exception: ' + traceback.format_exc())
            return 0

    def retrieve_merra2_data_for_dateRange(self, dateStart, dateEnd = None, dataset = "multiLevel"):
        """
        to retrieve data for a date range
        dateStart, dateEnd, strings: 'YYYYMMDD'
        """

        dstart = dparser.parse(dateStart)
        if dateEnd == None:
            dend = datetime.datetime.now()
        else:
            dend = dparser.parse(dateEnd)

        dt = dend - dstart
        count = dt.days +1
        dateList = list(rr.rrule(rr.DAILY,dstart, count= count))

        for dat in dateList:
            if self.verbose: print(dat)
            #quick patch for the problem of overloading MERRA2 servers with requests
            try:
                ret0 = self.retrieve_merra2_data_for_date(dat.strftime('%Y%m%d'), dataset)
                ret1 = self.retrieve_merra2_data_for_date(dat.strftime('%Y%m%d'),"singleLevel")
            except:
                print("MERRA2 servers not responding... trying again in 10 seconds")
                time.sleep(10)
                ret0 = self.retrieve_merra2_data_for_date(dat.strftime('%Y%m%d'),dataset)
                ret1 = self.retrieve_merra2_data_for_date(dat.strftime('%Y%m%d'),"singleLevel")
            if (ret0 == 0) or (ret1 == 0):
                dateList.pop(-1)

        print("Summary: Start: %s, End: %s : %d files available, last: %s"\
            %(dateStart, dateEnd, size(dateList), dateList[-1]))
                
        return dateList
        
    def createProfile(self, date,  simplify = 10,  plotFig=False):
        """
        date string: 'YYYYMMDD'
        simplify is a pressure in mB. at which we want to truncate and average the profile above that pressure level. if simplify = 0, do nothing.
        groundData = 1, 2, or 3. 
                     1: use groundWxData. For the moment, only available at SP
                     2: (Default) use MERRA2 2D data and interpolate to geopotential height. Available everywhere. 
                     3: Simply extrapolate the last 3 data points for the Temp, o3 and h20
        mixing ratio to the mean ground pressure. This sometimes causes negative h2o mixing ratio.
        """
        print("Converting MERRA2 data into a profile ...")

        #setup the plot
        if (plotFig): figure(1,figsize=(10,10)); clf()
        
        gnd = self.site['gndData']
        cldtype = self.site['cldtype']

        url, merraFilename_m = self.get_url_for_date(date,'multiLevel')
        url, merraFilename_s = self.get_url_for_date(date,'singleLevel')

        if os.path.exists(self.merraDir+merraFilename_m):
            print("Reading raw multiLevel MERRA2 data from %s"%merraFilename_m,)
            d = netCDF4.Dataset(self.merraDir+merraFilename_m)
            print('... Done')
        else:
            raise ValueError("Merra2 multiLevel data for %s is missing, retrieve it"%date)
        
        if os.path.exists(self.merraDir+merraFilename_s):
            print("Reading raw singleLevel MERRA2 data from %s"%merraFilename_s,)
            s = netCDF4.Dataset(self.merraDir+merraFilename_s)
            print('... Done')
        else:
            raise ValueError("Merra2 singleLevel data for %s is missing, retrieve it"%date)
        
        lat,lon=self.site['lat'],self.site['lon']
        lat0, lat1, lon0, lon1 = self.closestPoints(lat, lon)
        interpList = [(lat0,lon0), (lat0,lon1), (lat1, lon0), (lat1,lon1)]
        icoordList = [(0,0), (0,1), (1,0), (1,1)]
        timeList =  d.variables['time'][:]  # in minutes
        ntimesteps = timeList.size

        #find the lowest available common layer to all times and all 4 grid points
        q=[]
        for ti in range(ntimesteps):
            for (ilat, ilon) in icoordList:
                try:
                    q.append(min(flatnonzero(d.variables['T'][ti,:,ilat, ilon].mask == False)))
                except:
                    q.append(min(flatnonzero(d.variables['H'][ti,:,ilat, ilon] <1e10)))
        qb = max(q) # index of bottom layer

        dateStart = datetime.datetime.strptime(date,'%Y%m%d')
        pressureList =  d.variables['lev'][qb:]  # in mB
        nPressureLevels = pressureList.size

        # define the timevalues
        # 20190617: Increase the resolution here, change to use 2D time dim (1hr)
        timeList =  s.variables['time'][:]
        ntimesteps = timeList.size      
        dateTime = []
        for t in timeList:
            dateTime.append(dateStart+datetime.timedelta(seconds=t*60))

        layers = zeros([ntimesteps,nPressureLevels,7])  # 5 for P, T, xH20, xO3,  H. need to add QI (ice water), QL (liquid water)
        layers[:,:,0]= tile(pressureList,[ntimesteps,1])

        #For each time/pressure level, interpolate inwards.
        for t in range(ntimesteps):
            for l in range(nPressureLevels):
                dataPoints = [[] for i in range(6)]  #T, QV, O3, H
                for (ilat, ilon) in icoordList:
                  if t%3 == 0:
                    dataPoints[0].append(d.variables['T'][t/3,qb+l,ilat,ilon])
                    dataPoints[1].append(d.variables['QV'][t/3,qb+l,ilat,ilon])
                    dataPoints[2].append(d.variables['O3'][t/3,qb+l,ilat,ilon])
                    dataPoints[3].append(d.variables['H'][t/3,qb+l,ilat,ilon])
                    dataPoints[4].append(d.variables['QL'][t/3,qb+l,ilat,ilon])
                    dataPoints[5].append(d.variables['QI'][t/3,qb+l,ilat,ilon])
                  elif t%3 == 1:
                    dataPoints[0].append(d.variables['T'][(t-1)/3,qb+l,ilat,ilon]*s.variables['T2M'][t,ilat,ilon]/s.variables['T2M'][t-1,ilat,ilon])
                    dataPoints[1].append(d.variables['QV'][(t-1)/3,qb+l,ilat,ilon]*s.variables['TQV'][t,ilat,ilon]/s.variables['TQV'][t-1,ilat,ilon])
                    dataPoints[2].append(d.variables['O3'][(t-1)/3,qb+l,ilat,ilon]) # assume the edge height stays the same
                    dataPoints[3].append(d.variables['H'][(t-1)/3,qb+l,ilat,ilon]) # assume the edge height stays the same
                    dataPoints[4].append(d.variables['QL'][(t-1)/3,qb+l,ilat,ilon]*s.variables['TQL'][t,ilat,ilon]/s.variables['TQL'][t-1,ilat,ilon])
                    dataPoints[5].append(d.variables['QI'][(t-1)/3,qb+l,ilat,ilon]*s.variables['TQI'][t,ilat,ilon]/s.variables['TQI'][t-1,ilat,ilon])
                  elif t%3 == 2:
                    dataPoints[0].append(d.variables['T'][(t-2)/3,qb+l,ilat,ilon]*s.variables['T2M'][t,ilat,ilon]/s.variables['T2M'][t-2,ilat,ilon])
                    dataPoints[1].append(d.variables['QV'][(t-2)/3,qb+l,ilat,ilon]*s.variables['TQV'][t,ilat,ilon]/s.variables['TQV'][t-2,ilat,ilon])
                    dataPoints[2].append(d.variables['O3'][(t-2)/3,qb+l,ilat,ilon]) # assume the edge height stays the same
                    dataPoints[3].append(d.variables['H'][(t-2)/3,qb+l,ilat,ilon]) # assume the edge height stays the same
                    dataPoints[4].append(d.variables['QL'][(t-2)/3,qb+l,ilat,ilon]*s.variables['TQL'][t,ilat,ilon]/s.variables['TQL'][t-2,ilat,ilon])
                    dataPoints[5].append(d.variables['QI'][(t-2)/3,qb+l,ilat,ilon]*s.variables['TQI'][t,ilat,ilon]/s.variables['TQI'][t-2,ilat,ilon])
                TInterp = interpolate.interp2d(zip(*interpList)[0],zip(*interpList)[1],dataPoints[0])
                QVInterp = interpolate.interp2d(zip(*interpList)[0],zip(*interpList)[1],dataPoints[1])
                O3Interp = interpolate.interp2d(zip(*interpList)[0],zip(*interpList)[1],dataPoints[2])
                HInterp = interpolate.interp2d(zip(*interpList)[0],zip(*interpList)[1],dataPoints[3])
                QLInterp = interpolate.interp2d(zip(*interpList)[0],zip(*interpList)[1],dataPoints[4])
                QIInterp = interpolate.interp2d(zip(*interpList)[0],zip(*interpList)[1],dataPoints[5])
                layers[t,l,1] = TInterp(lat,lon)
                qv = QVInterp(lat,lon)
                layers[t,l,2] = 1e6 * (28.964 / 18.015) * qv / (1.0 - qv) # h20 mixing ratio PPM
                o3 = O3Interp(lat,lon)
                layers[t,l,3] = 1e6 * (28.964 / 47.997) * o3 # O3 mixing ratio PPM
                layers[t,l,4] = HInterp(lat,lon)
                layers[t,l,5] = QLInterp(lat,lon)
                layers[t,l,6] = QIInterp(lat,lon)

        # redefine a new array to hold layers
        # and include one more for extrapolated layer
        layers_ex = zeros([ntimesteps,nPressureLevels+1,7])
        last_layer=zeros(7)
        minindx = nPressureLevels #eventually to resize the layer_ex array to contain only up to surface

        #incorporating surface data
        # for each time step
        for ti in range(ntimesteps):
            date = datetime.datetime.strftime(dateTime[ti],'%Y%m%d')
            if (gnd == 1):
               #retrieve data from the relevant Wx file
                time = datetime.datetime.strftime(dateTime[ti],'%H%M%S')
                wxFilename = '%s_%s_wx_keck.txt'%(date,time)
                # wx = self.readWxData(filename=wxFilename)
                wx = self.readSPWxData(date=date, time=time, shrink=1)
                Pgnd = mean(wx['presmB'])
                Tgnd = mean(wx['tempC'])+273.15  # from C to K.
                Rhgnd = mean(wx['rh'])
                print("Mean Ground Pres: %3.2f mBar"%Pgnd)
                print("Mean Ground Temp: %3.2f K"%Tgnd)
                print("Mean Ground RH: %3.2f%%"%Rhgnd)

                last_layer[0]=Pgnd
                Vmrgnd = self.calcVmr(Pgnd,Tgnd,Rhgnd)
                last_layer[1] = Tgnd
                last_layer[2] = 1e6*Vmrgnd
                last_layer[3] = layers[ti,0,3]
                last_layer[4] = layers[ti,0,4]

            elif (gnd == 3):
                # fit a 1d polynomial to lowest 2 layersand extrapolate on it.
                for i in range(1,5):
                    fit = polyfit(log(layers[ti,0:2,0]),layers[ti,0:2,i],1)
                    f = poly1d(fit)
                    # only works with later version of scipy.
                    # f = interp1d(layers[0,0:3,0],layers[0,0:3,i],fill_value='extrapolate')
                    last_layer[i] = f(log(Pmean))

            elif (gnd == 2):
            #we first interpolate vertically and then inwards.
                Hgnd = self.site['alt']
                try:
                    P = interpolate.interp1d(layers[ti,:,4], log10(layers[ti,:,0]),fill_value = "extrapolate")
                    Pgnd = 10**P(Hgnd)
                    #print layers[ti,0:2,4]
                    #print layers[ti,0:2,0]
                    #print Hgnd, Pgnd
                except:
                    Pgnd = self.calcPgnd(Hgnd)
                   
                dataPoints = [[] for i in range(2)] #T and QV
                for (ilat, ilon) in icoordList:
                    gnd_pressures = insert(d.variables['lev'][qb:], 0, [s.variables['PS'][ti,ilat, ilon] / 100])
                    if ti%3 == 0:
                      gnd_temps = insert(d.variables["T"][ti/3, qb:, ilat, ilon], 0, [s.variables['T2M'][ti, ilat, ilon]])
                      gnd_wv = insert(d.variables["QV"][ti/3, qb:, ilat, ilon], 0, [s.variables['QV2M'][ti, ilat, ilon]])
                    if ti%3 == 1:
                      gnd_temps = insert(d.variables["T"][(ti-1)/3, qb:, ilat, ilon]*s.variables['T2M'][ti,ilat,ilon]/s.variables['T2M'][ti-1,ilat,ilon],
                                         0, [s.variables['T2M'][ti, ilat, ilon]])
                      gnd_wv = insert(d.variables["QV"][(ti-1)/3, qb:, ilat, ilon]*s.variables['TQV'][ti,ilat,ilon]/s.variables['TQV'][ti-1,ilat,ilon],
                                         0, [s.variables['QV2M'][ti, ilat, ilon]])
                    if ti%3 == 2:
                      gnd_temps = insert(d.variables["T"][(ti-2)/3, qb:, ilat, ilon]*s.variables['T2M'][ti,ilat,ilon]/s.variables['T2M'][ti-2,ilat,ilon],
                                         0, [s.variables['T2M'][ti, ilat, ilon]])
                      gnd_wv = insert(d.variables["QV"][(ti-2)/3, qb:, ilat, ilon]*s.variables['TQV'][ti,ilat,ilon]/s.variables['TQV'][ti-2,ilat,ilon],
                                         0, [s.variables['QV2M'][ti, ilat, ilon]])
                    temp = interpolate.interp1d(gnd_pressures, gnd_temps , fill_value = "extrapolate")  
                    dataPoints[0].append(temp(Pgnd))
                    qv = interpolate.interp1d(gnd_pressures,gnd_wv , fill_value = "extrapolate")
                    dataPoints[1].append(qv(Pgnd))
                TInterp = interpolate.interp2d(zip(*interpList)[0],zip(*interpList)[1],dataPoints[0])
                QVInterp = interpolate.interp2d(zip(*interpList)[0],zip(*interpList)[1],dataPoints[1])
                Tgnd = TInterp(lat,lon)
                QVgnd = QVInterp(lat,lon)
                last_layer[0] = Pgnd
                last_layer[1] = Tgnd
                last_layer[2] = 1e6 * (28.964 / 18.015) * QVgnd / (1.0 - QVgnd)
                last_layer[3] = layers[ti,0,3]
                last_layer[4] = Hgnd
                #current just leave the cloud layers at 0 on the ground. potentially we'll want to change this.

            else:
                print("ERROR: No height provided or no surface level data available")
                return 0
            
            #print Pgnd, layers[ti,:,0]
            if Pgnd > layers[ti,0,0]:
                #minindx = 0
                layers_ex[ti,:,:] = concatenate([[last_layer],layers[ti,:,:]])
            else:
                #insert the ground layer at the right position in the pressure list
                #also delete the layers that are "below" the ground layer (usually because it is on a peak).
                indx = where(layers[ti,:,0] < Pgnd)[0][0]
                #if (indx < minindx):
                #    maxindx = indx
                #layers_ex[ti,:nPressureLevels-indx+1,:] = concatenate([[last_layer],layers[ti,indx:,:]])
                layers_ex[ti,:nPressureLevels,:] = concatenate([np.tile(last_layer,(indx,1)),layers[ti,indx:,:]])
                ## introduce indx # of copies of last_layer (which doesn't affect am's results)
                ## this hack ensures that the simplifying of the layers doesn't improperly reshape the layer array.
                #layers_ex[ti,:,:]= insert(layers[ti,:,:],1,[last_layer],0)

            # create plot for the day
            if (plotFig):
                subplot(2,2,1)  # temperature
                plot(layers_ex[ti,:,1], layers_ex[ti,:,0],'.-')
                #plot([last_layer[1]],[last_layer[0]],'or')
                if ti == ntimesteps-1:
                    ylim([750,1e-1])
                    grid()
                    yscale('log')
                    xlabel('Temp [K]')
                    ylabel('Press [mB]')

                subplot(2,2,2)  # height
                plot(layers_ex[ti,:,4], layers_ex[ti,:,0],'.-')
                #plot([last_layer[4]],[last_layer[0]],'or')
                if ti == ntimesteps-1:
                    grid()
                    yscale('log')
                    ylim([750,1e-1])
                    xlabel('height [m]')
                    ylabel('Press [mB]')

                subplot(2,2,3)  # xh2o
                plot(layers_ex[ti,:,2], layers_ex[ti,:,0],'.-')
                #plot([last_layer[2]],[last_layer[0]],'or')
                if ti == ntimesteps-1:
                    #yscale('log')
                    ylim([750,1e-1])
                    grid()
                    xlabel('h2o mixing ratio [ppm]')
                    ylabel('Press [mB]')

                subplot(2,2,4)  # xO3
                plot(layers_ex[ti,:,3], layers_ex[ti,:,0],'.-')
                #plot([last_layer[3]],[last_layer[0]],'or')
                if ti == ntimesteps-1:
                    ylim([750,1e-1])
                    grid()
                    yscale('log')
                    xlabel('O3 mixing ratio [ppm]')
                    ylabel('Press [mB]')
                    suptitle('MERRA2-based amc profile for %s for %s, gndData:%d '%(self.site['name'],date,gnd))
                    savefig(self.amcDir+'%s_gndData%d_amcProfile.png'%(date,gnd))

        #layers_ex = layers_ex[layers_ex != 0]
        profile={'layers': layers_ex, 'time':dateTime, 'gndData':gnd, 'cldtype':cldtype}
        if simplify != 0 :
            profile = self.simplifyProfile(profile, simplify)

        return profile

    def calcPgnd(self,Hgnd):
        """
        calc Pgnd given Hgnd
        """
        Pgnd = 10* 101.29 * ((15.04 - .00649 * Hgnd + 273.1)/288.08)**5.256
        return Pgnd
    
    def calcVmr(self, P, T, Rh):
        """
        based on h2o_psat.c from am 9.01
        First Computes H2O saturation vapor pressure over a flat liquid  surface at temperature T.
        Then from Psat, Pgnd, and Rh, calculates the VMR of water.
        T in K
        Rh in %
        P in mB
        """

        # The following formula gives log(Psat), with Psat in Pa.
        ln_T = log(T);
        r1 = 54.842763 - 6763.22 / T - 4.210   * ln_T + 0.000367 * T;
        r2 = 53.878    - 1331.22 / T - 9.44523 * ln_T + 0.014025 * T;
        ln_Psat = r1 + r2 * tanh(0.0415 * (T - 218.8));
        # convert Psat to mB
        PsatmB = 0.01 * exp(ln_Psat);
        vmr =  0.01 * Rh * (PsatmB / P)

        return vmr

    def simplifyProfile(self, profile, lastLayerPressure=10):
        """
        we remove all layers with  Pbase < lastLayerPressure in mB
        Replace: pressure, temperature and height with those of the last layer
        replace xH2o and xO3 with their weighted mean
        weights is the pressure diff within a layer
        """

        print("Simplifying pressure down to Pbase = %s mB"%lastLayerPressure)
        l = profile['layers']
        nPressureLevels = l[0,:,0].size
        dp =[]
        for i in range(1,nPressureLevels):
            dp.append(l[0,i-1,0]-l[0,i,0])
        dp.append(l[0,-1,0])

        # find layer with Pbase = 10mB
        q = flatnonzero(l[0,:,0] >= lastLayerPressure)
        ilim = q[-1]
        w = dp[q[-1]:]
        # redefine the xh2o and xO3
        l[:,ilim,2]=average(l[:,ilim:,2],weights=w,axis=1)
        l[:,ilim,3]=average(l[:,ilim:,3],weights=w,axis=1)
        l[:,ilim,5]=average(l[:,ilim:,5],weights=w,axis=1)
        l[:,ilim,6]=average(l[:,ilim:,6],weights=w,axis=1)

        l=l[:,0:ilim+1,:]
        profile['layers'] = l

        return profile

    def profile2am(self, profile):
        """
        Make an am model configuration file from a given layers
        and save it.

        """
        print("Converting given profile into amc files ...")

        # define set parameters
        P_ground = 1000.0       # surface pressure [mbar]
        P_tropopause = 100.0    # tropopause pressure [mbar]
        P_stratopause = 1.0     # stratopause pressure [mbar]
        P_Voigt = 1.1   # Pressure [mbar] below which Voigt lineshape will be used
        H2O_SUPERCOOL_LIMIT = 228. # ice assumed below this temperature [K]

        amcList=[]
        #loop over all the timesteps
        for t_idx, t in enumerate(profile['time']):

            layers = profile['layers'][t_idx,:,:]
            nlayers = shape(layers)[0]

            tstr = datetime.datetime.strftime(t,'%Y%m%d_%H%M%S')
            amcList.append('%s_%s_gndData%d_%s.amc'%(self.site['name'],tstr, profile['gndData'],profile['cldtype']))
            f = open(self.amcDir+amcList[-1],'w')
            print("Saving amc file to:  %s"%(self.amcDir+amcList[-1]))

            # print header
            f.write('? \n')
            f.write('? Usage: \n')
            f.write('?   am this_file f_min f_max df zenith_angle trop_h2o_scale_factor\n')
            f.write('? Example: \n')
            f.write('?   am this_file  0 GHz  300 GHz  10 MHz  0 deg 1.0 \n')
            f.write('? \n')
            f.write('f %1 %2  %3 %4  %5 %6 \n')
            f.write('output f GHz tau Tb K Trj K \n')
            f.write('za %7 %8 \n')
            f.write('tol 1e-4 \n')
            f.write('\n')
            f.write('Nscale troposphere h2o %9 \n')
            f.write('\n')
            f.write('T0 2.7 K \n')

            # loop over all levels in the layer
            for l in range(nlayers):
                layer = layers[-1-l,:]

                P = layer[0]
                T = layer[1]
                if P > P_ground:
                    print ("pressure greater than P_ground, skipping layer")
                    continue
                if (P < P_stratopause):
                    layer_type = 'mesosphere'
                elif (P < P_tropopause):
                    layer_type = 'stratosphere'
                else:
                    layer_type = 'troposphere'

                # replace x_H2O and x_O3 with their layer midpoint values,
                # except for the top layer.
                if l == 0:
                    xh2o = 1e-6*layers[-1-l,2]
                    xo3 =  1e-6*layers[-1-l,3]
                    if layers[-1-l,4] > 1e-6 and layers[-1-l,4] > 1e-6:
                        xh2ol = (100. * (layers[-2-l,0] - layers[-1-l,0])/9.80665) * layers[-1-l,5]
                        xh2oi = (100. * (layers[-2-l,0] - layers[-1-l,0])/9.80665) * layers[-1-l,6]
                    else: 
                        xh2oi = 0.0
                        xh2ol = 0.0
                else:
                    xh2o = 1e-6*(layers[-1-l,2]+layers[-l,2])/2.
                    xo3 =  1e-6*(layers[-1-l,3]+layers[-l,3])/2.
                    if layers[-1-l,4] > 1e-6 and layers[-1-l,4] > 1e-6:
                        xh2ol = (100. * (layers[-1-l,0] - layers[-l,0])/9.80665) * (layers[-1-l,5]+layers[-l,5])/2.
                        xh2oi = (100. * (layers[-1-l,0] - layers[-l,0])/9.80665) * (layers[-1-l,6]+layers[-l,6])/2.
                    else: 
                        xh2oi = 0.0
                        xh2ol = 0.0

                f.write('layer %s \n'%layer_type)
                f.write('Pbase %3.1f mbar \n' %P)
                f.write('Tbase %3.2f K \n' %T )
                if (P <= P_Voigt):
                    f.write('lineshape Voigt-Kielkopf \n')
                f.write('column dry_air vmr \n' )
                f.write('column h2o vmr %1.3e \n' %xh2o)
                f.write('column o3 vmr %1.3e \n'%xo3)
                #if you want to have ice/water in the profiles. default is just  h2o vapor.
                if profile['cldtype'] == 'icewater': 
                    if (T > H2O_SUPERCOOL_LIMIT):
                        f.write('column lwp_abs_Rayleigh %2.3f kg*m^-2 \n'%xh2ol)
                    else:
                        f.write('column iwp_abs_Rayleigh %2.3f kg*m^-2 \n'%xh2oi)
                f.write('\n')

        f.close()
        return amcList

    def findFailedSpectra(self,year):
        failedSpectra = []
        outfiles = glob.glob(self.amcDir+'*'+year+'*'+'gndData2*.out')
        for outfile in outfiles:
            if os.stat(outfile).st_size == 0:
                failedSpectra.append(outfile)
                os.remove(outfile)
                amcfile = outfile.replace('.out','.amc')
                os.remove(amcfile)
        print (size(failedSpectra))
        return failedSpectra
    
    def findMissingSpectra(self,year):
        missingSpectraDays = []
        
        start = datetime.datetime.strptime('%s0101'%year,'%Y%m%d')
        dateList = [start + datetime.timedelta(days = i) for i in range(365)]
        for date in dateList:
            datestr = date.strftime('%Y%m%d')
            nspecPerDay = size(glob.glob(self.amcDir+'*%s*gndData2*.out'%datestr))
            if nspecPerDay  != 8:
                print ('%s: only %d of 8 amc spectra in that day'%(datestr,nspecPerDay))
                missingSpectraDays.append(date)

        print('%d days with missing spectra'%size(missingSpectraDays)  )
 
        return missingSpectraDays
     
    def findIceWater(self, year):
        """find files that have icewater"""
        icewater = []
        amcfiles = glob.glob(self.amcDir+'*_'+year+'*'+'.amc')
        for amcfile in amcfiles:
            if 'wp_abs_Rayleigh' in open(amcfile).read():
                icewater.append(amcfile)
                os.remove(amcfile)
                outfn = amcfile.replace('.amc','.out')
                os.remove(outfn)
        print(size(icewater))
        return icewater

    def run_am(self, amcFileName, za=0., f0=50., f1=300., df=100, h2o_scale=1.0):
        '''
        f0: starting frequency [GHz]
        f1: ending frequency [GHz]
        df: frequency step size [MHz]
        za = zenigh angle [deg] (za = 90 - el)
        h2o_scale = trop_h2o_scale_factor [unitless]

        '''
        fUnits  = 'GHz'
        dfUnits = 'MHz'
        zaUnits = 'deg'

        # run am
        outfn, errfn = self.drive_am(amcFileName, [f0,fUnits,
                                       f1,fUnits,
                                       df,dfUnits,
                                       za,zaUnits,
                                       h2o_scale])


        amOut, pwv, dry_air,o3, Tgnd, Pgnd = self.readAmResults(outfn)
        fs = amOut['f']
        tb = amOut['Tb']
        trj = amOut['Trj']
        tau = amOut['tau']
        
        # interpolate to 10MHz frequency.
        # df = f1 - f0
        # fstep = 0.01
        # nf = int(df / fstep)
        # fs2 = np.linspace(f0,f1,nf)
        # tb2 = np.interp(fs2,fs,tb)
        # trj2 = np.interp(fs2,fs,trj)

        return fs,tb,trj, pwv, tau

    def drive_am(self, infn, pars):
        """

        """
        outfn = infn.replace('.amc','.out')
        errfn = infn.replace('.amc','.err')

        if os.path.exists(self.amcDir+outfn):
            if os.path.getsize(self.amcDir+outfn) != 0:
                print (self.amcDir+outfn+" already exists, skip running drive_am on this... ")
                return outfn, errfn
        ##TODO: Remove the hard-coded path (likely by creating an alias in the predict_tsky_web).
        #cmd = 'am %s%s '%(self.amcDir,infn)
        cmd = './am %s%s '%(self.amcDir,infn) # 20190617: Change to the local am
        if self.debug:
            print ("am parameters: f0: %.2f %s, "
                   "f1: %.2f %s,  df: %.2f %s, "
                   "za: %2.2f %s, h2o_scale=%2.2f"
                   %(pars[0],pars[1],pars[2],pars[3],pars[4],
                     pars[5],pars[6],pars[7],pars[8]))
        cmd += ("%.2f %s %.2f %s %.2f %s %.2f %s %.2f"
                %(pars[0],pars[1],pars[2],pars[3],pars[4],
                  pars[5],pars[6],pars[7],pars[8]))

        cmd +=  ' >%s%s '%(self.amcDir,outfn)
        cmd += ' 2>%s%s '%(self.amcDir,errfn)
        if self.debug: print (cmd)

        subprocess.Popen(cmd,shell=True).wait()
        return outfn, errfn

    def readAmResults(self,outfn):
        """
        Reads the results of the .out as well
        as well as the total from the .err for diagnostics
        """
        
        try:
            Pbase = []
            Tbase = []

            inErrFile = open(self.amcDir+outfn.replace('.out','.err'),'r')
            for line in inErrFile:
                if line.startswith('output'):
                    outfmt = line.split()[1::2]
                if line.startswith('Pbase'):
                    Pbase.append(float(line.split()[1]))
                if line.startswith('Tbase'):
                    Tbase.append(float(line.split()[1]))   
                if line[0] == '#' and 'dry_air' in line:
                    dry_air = float(line.split()[2])
                if line[0] == '#' and 'h2o' in line:
                    pwv = float(next(inErrFile).split()[1][1:])
                if line[0] == '#' and 'o3' in line:
                    o3 = float(line.split()[2])
            Tgnd = Tbase[-1]
            Pgnd = Pbase[-1]

            amOut = genfromtxt(self.amcDir+outfn,delimiter='',names=outfmt)
            
        except:
            print( "ERROR: could not read file output spectra for %s"%outfn)

        return  amOut, pwv, dry_air, o3, Tgnd, Pgnd

    def readBandpass(self, bandopt):
        """
        bandName can be:
           - BK95, BK150, BK210, BK220, BK270 for pre-selected measured bandpasses
           - tipper850 for Simon's 850um tipper
           - custom, then it also needs v_cen and frac_bw
        """
        if 'name' not in bandopt:
            bandopt['name'] = 'BK150'
            
        # produce the tipper850 instrument bandpass
        if bandopt['name'] == 'tipper850':
            #single peaked gaussian
            c = 857 #band center in GHz 857
            s = 43.74 #FWHM of 103 GHz in total 43.74
            f=600+arange(600) #range 600,1200
            Aflat = exp(-(f-c)**2/(2*s**2))

            #two half gaussians with a flat top of width 36GHz, centered on 850GHz. low-freq FWHM = 60GHz, high-freq = 74.
            #here are some hacks to ensure that the two Gaussians reach 0.8 at the point where the flat peak starts.
            c1 = 849 #849 #  mean of low frequency Gaussian
            c3 = 847 #847 # mean of high frequency Gaussian
            s1 = 25.48 # sigma in  GHz (FWHM = 60GHz) (sigma = 25.48 GHz)
            s3 = 31.42 #FWHM = 74GHz (sigma = 31.42 GHz)
            f1=600+arange(232) #low-frequency half gaussian #600+arange(232)
            f2 = 832 + arange(36) #flat middle peak #832 + arange(36)
            f3 = 868 + arange(300) #high-frequency half gaussian #868 + arange(300)
            f = concatenate([f1,f2,f3])

            Aflat1 = exp(-(f1-c1)**2/(2*s1**2))
            Aflat2 = full(36,0.8) #full(36,0.8)
            Aflat3 = exp(-(f3-c3)**2/(2*s3**2))
            Aflat = concatenate([Aflat1,Aflat2,Aflat3])

            return f,Aflat, Aflat

        # read and produce the BK bands bandpass
        elif bandopt['name'].startswith('BK'):
            if bandopt['name'] == 'BK150':
                filename = 'B2_frequency_spectrum_20141216.txt'
                cutoff=[110,190]
            elif bandopt['name'] == 'BK95':
                filename = 'K95_frequency_spectrum_20150309.txt'
                cutoff=[68,120]
            elif bandopt['name'] == 'BK210':
                filename = 'K210_frequency_spectrum_20160101.txt'
                cutoff=[150,300]
            elif bandopt['name'] == 'BK220':
                filename = 'K220_frequency_spectrum_20160120.txt'
                cutoff=[150,300]
            elif bandopt['name'] == 'BK270':
                filename = 'K270_frequency_spectrum_20170710.txt'
                cutoff=[200,340]
            else:
                print("bandName must be 'BK95', 'BK150', 'BK210', 'BK220', or 'BK270'")
                return 0
            e = genfromtxt(self.auxDataDir+filename, delimiter=',', comments='#')
            f = e[:,0]
            if shape(e)[1] == 3:
                Aflat = e[:,1]
                Arj = e[:,2]
                if Arj[0] == 0:
                    Arj = Aflat
            else:
                Aflat = e[:,1]
                Arj = e[:,1]
            q1 = flatnonzero(f > cutoff[1])
            q2 = flatnonzero(f < cutoff[0])
            q = concatenate((q1,q2))
            Aflat[q] = 0
            Arj[q] = 0
            
            return f, Aflat, Arj
        
        elif bandopt['name'].startswith('EHT'):
            if bandopt['name'] == 'EHT_hi':
                c=228.1
                bw = 4.0
            elif bandopt['name'] == 'EHT_lo':
                c=214.1
                bw = 4.0
            f = arange(0,1200,0.1)
            Aflat = ones(12000)
            q1 = flatnonzero(f > c+bw/2.)
            q2 = flatnonzero(f < c-bw/2.)
            q = concatenate((q1,q2))
            Aflat[q] = 0

            return f,Aflat, Aflat

        # produce a custom tophat bandpass
        elif bandopt['name'] == 'custom':
            c = bandopt['v_cen']    # band center[GHz]
            bw = bandopt['frac_bw'] # fractional bandwidth
            f = arange(0,1200,0.1)
            Aflat = ones(12000)
            q1 = flatnonzero(f > c+c*bw/2.)
            q2 = flatnonzero(f < c-c*bw/2.)
            q = concatenate((q1,q2))
            Aflat[q] = 0

            return f,Aflat, Aflat

    def integBand(self,fband,band,freq,trj):
        """
        inputs are:
        fband, band: frequency axis and bandpass peak normalized
        freq, tb: frequency axis and atmospheric spectra
        """

        # interpolate bandpass onto Tb's frequency axis
        bandn = interp(freq,fband,band)

        #take the weighted average
        wint = sum(trj*bandn)/sum(bandn)

        return wint

    def index_nearest(self,items, pivot, n=1):
        """
        given a list of datetimes items 
        and a desired number = pivot 
        n is the number of nearest we want to find. n=1 is nearest
        return the index of the element closest to pivot
        copied from http://stackoverflow.com/questions/32237862/find-the-closest-date-to-a-given-date
        
        """
        import heapq
        ivalues = []
        try:
            values = heapq.nsmallest(n,items, key=lambda x: abs(x - pivot))
            for val in values:
                ivalues += [find(items == val)[0]]
            return ivalues, values
        except:
            print("failed to find nearest index of single date")
            pass

    def printSingle(self, tstart, Tsky_merra):
        """
        prints the results if single is requested
        """
        dstart = dparser.parse(tstart)
        inear,t_merra_near = self.index_nearest(array(Tsky_merra['t']),dstart,2)
        i0 = inear[0]
        i1 = inear[1]
        print("Requested Tsky at %s"%tstart)
        print('MERRA2-derived Tsky for %s and %s'%(t_merra_near[0].isoformat(),t_merra_near[1].isoformat()))
        keys = ['BK100','BK150','BK210','BK220','BK270','tipper850','EHT_lo','EHT_hi'] 
        print('Band:%10s  PWV = %5.1f um,     PWV = %5.1f um  '%('PWV',Tsky_merra['pwv'][i0],Tsky_merra['pwv'][i1]))
        for key in keys:
            print('Band:%10s  Tsky= %5.1f K,     Tsky= %5.1f K  '%(key,Tsky_merra[key][inear[0]],Tsky_merra[key][inear[1]]))
    

    def save2pickle(self, Tsky_merra, opts):
        """
        function to save results to pickle file
        """
        import pickle

        print("Saving results to pickle file...")
        now = opts['now']
        strnow = now.isoformat()
        bandopt =  opts['bandopt']
        siteopt = opts['siteopt']
        dateopt = opts['dateopt']
        
        keys = ['tstr','BK100','BK150','BK210','BK220','BK270','tipper850']
        if 'custom' in  bandopt['name']:
            keys.append('custom')
        keys=keys+['pwv']
        
        if opts['web']== False:
            self.webDir = self.finalDir
        if dateopt['type'] == 'datelist':
            dateInfo = '%s_%s'%(dateopt['start'],dateopt['end'])

        pickleFile = 'Tsky_merra2_output_%s_%s.pickle'%(siteopt['name'],dateInfo)

        print("Saving resuls to %s"%(self.webDir+pickleFile))
        with open(self.webDir+pickleFile,'w') as f:
            pickle.dump([Tsky_merra,opts], f)

        return pickleFile

    def save2csv(self, Tsky_merra, opts):
        """
        function to save the results to a csv file

        """
        print("Saving results to csv file...")
        now = opts['now']
        strnow = now.isoformat()
        bandopt =  opts['bandopt']
        siteopt = opts['siteopt']
        dateopt = opts['dateopt']

        if opts['web'] == False:
            self.webDir = self.finalDir
        
        keys = ['tstr','BK100','BK150','BK210','BK220','BK270','tipper850']
        if 'custom' in  bandopt['name']:
            keys.append('custom')
        keys=keys+['pwv']
        
        fmt = '{:8.3f}, '*(len(keys)-1)
        fmt = '{:10}, '+fmt
        fmt = fmt+'\n'

        if dateopt['type'] == 'datelist':
            dateInfo = '%s_%s'%(dateopt['start'],dateopt['end'])

        csvFile = 'Tsky_merra2_output_%s_%s.csv'%(siteopt['name'],dateInfo)
        
        print("Saving resuls to %s"%(self.webDir+csvFile))
        with open(self.webDir+csvFile, 'w') as f:
            f.write('# Output from predict_Tsky.py generated on %s \n'%strnow)
            f.write('# Summary of user inputs \n')
            f.write('# Band options:  name=\'%s\''%bandopt['name'])

            if 'v_cen' in bandopt.keys():
                f.write(', v_cen= %5.1f, frac_bw=%3.1f'%(bandopt['v_cen'],bandopt['frac_bw']))
            else:
                f.write('\n')

            f.write('# Site options:  type=\'%s\', name=\'%s\''%(siteopt['name'],siteopt['type']))
            f.write(', shortname=\'%s\', lat=%6.2f [deg], lon=%6.2f [deg], alt=%6.2f [m] \n'
                    %(siteopt['sname'],siteopt['lat'],siteopt['lon'],siteopt['alt']))

            if dateopt['type'] == 'datelist':
                f.write('# Date options:  type=\'%s\', start=\'%s\', end=\'%s\' \n'
                        %(dateopt['type'],dateopt['start'],dateopt['end']))
            f.write("#\n")
            f.write('# All outputs are in Kelvin_RJ, pwv in units of microns. \n')
            f.write('# ')
            f.write(',  '.join(keys))
            f.write('\n')
            f.write('#\n')

            output_rows = zip(*[Tsky_merra[key] for key in keys])
            for row in output_rows:
                f.write(fmt.format(*row))
            f.close()

            return csvFile

    def save2plot(self,Tsky_merra,opts):
        """
        function to plot Tsky and PWV in two saved plots 

        """
        print('Plotting results...')
        ioff()
        now = opts['now']
        strnow = now.isoformat()
        bandopt =  opts['bandopt']
        siteopt = opts['siteopt']
        dateopt = opts['dateopt']

        #only the bands that Tsky will be plotted for. PWV and time will be called separately.
        plottedBands = ['BK100','BK150','BK220','BK270'] 
        if 'custom' in  bandopt['name']:
            plottedBands.append('custom')

        dateInfo = '%s_%s'%(dateopt['start'],dateopt['end'])
        tskyPlot = self.webDir+'Tsky_merra2_plot_%s_%s.png'%(siteopt['name'],dateInfo)

        figure(figsize=(20,10))
        subplot(2,1,1)
        for band in plottedBands:
            if band == 'custom' and 'v_cen' in bandopt.keys():
                plot(Tsky_merra['t'],Tsky_merra[band],'.-',label='custom, v_cen= %5.1f, frac_bw=%3.1f'%(bandopt['v_cen'],bandopt['frac_bw'])) 
            else:
                plot(Tsky_merra['t'],Tsky_merra[band],'.-',label=band)
        legend()
        grid()
        ylim([0,max([max(Tsky_merra[band]) for band in plottedBands])+2]) #limit the y-axis at 0 and the maximum Tsky + 2. 
        ylabel(r"Brightness Temperature $T_{rj}$ [K]")
        xlabel("Date")
        title("Rayleigh-Jeans Zenith Brightness Temperature at %s from %s to %s"%(siteopt['name'],dateopt['start'],dateopt['end']))

        subplot(2,1,2)
        plot(Tsky_merra['t'],array(Tsky_merra['pwv'])/1000.,'.-',label='pwv') 
        grid()
        ylabel("PWV [mm]")
        xlabel("Date")
        title("Precipitable Water Vapor at %s from %s to %s"%(siteopt['name'],dateopt['start'],dateopt['end']))
        savefig(tskyPlot)
        close()
        
        return tskyPlot


#    data is the data read from a single netcdf4 file using d = netCDF4.Dataset("filename")
#    d.variables['T'].dimensions ---> (u'time', u'lev', u'lat', u'lon')
#levs is pressure levels in mb
#In [146]: for k in d.variables.keys():
#    print k, d.variables[k].units, d.variables[k].standard_name
#   .....:
#EPV K m+2 kg-1 s-1 ertels_potential_vorticity
#H m edge_heights
#O3 kg kg-1 ozone_mass_mixing_ratio
#OMEGA Pa s-1 vertical_pressure_velocity
#PHIS m+2 s-2 surface geopotential height
#PS Pa surface_pressure
#QI kg kg-1 mass_fraction_of_cloud_ice_water
#QL kg kg-1 mass_fraction_of_cloud_liquid_water
#QV kg kg-1 specific_humidity
#RH 1 relative_humidity_after_moist
#SLP Pa sea_level_pressure
#T K air_temperature
#U m s-1 eastward_wind
#V m s-1 northward_wind
#lat degrees_north latitude
#lev hPa air_pressure
#lon degrees_east longitude
#time minutes since 2016-01-01 00:00:00 time

#import scipy.io as sio
#
#mat_contents = sio.loadmat('data.mat')
#print mat_contents
#To write a .mat file:
#
#import numpy as np
#import scipy.io as sio
#
#vect = np.arange(10)
#print vect.shape
#sio.savemat('data.mat', {'vect':vect})

# 0.5 deg in lat and 5/8deg in lon
#nc4 is netcdf file, nco tools, http://nco.sourceforge.net/
# http://schubert.atmos.colostate.edu/~cslocum/netcdf_example.html
