import urllib2
import StringIO
import gzip
import datetime
import calendar
import numpy as np
from dateutil.rrule import *

# Run from the merra2_work directory, and save APEX wx files in merra2_products/wx_data.
# Save monthly files, with data for each minute.

keys = [
    't',
    'temp',
    'pres',
    'rh',
]  # Datetime, Temperature, Pressure, Relative Humidity. Need to find units.

for year in range(2006, 2017):
    for month in range(1, 13):
        wx = {}
        for i in keys:
            wx[i] = []
        tempData = [
            [],
            [],
        ]  # holds data for the entire month in case files are saved on wrong dates.
        presData = [[], []]
        rhData = [[], []]
        # for every month, we want to store a file with the data from all the days in that month.
        for day in range(1, calendar.monthrange(year, month)[1] + 1):
            date = datetime.datetime(year, month, day)
            datestr = date.strftime("%Y-%m-%d")

            f = urllib2.urlopen(
                "http://www.apex-telescope.org/weather/Historical_weather/temperature/%s_temperature.log.gz"
                % datestr
            )
            compressedFile = StringIO.StringIO()
            compressedFile.write(f.read())
            compressedFile.seek(0)
            decompressedFile = gzip.GzipFile(fileobj=compressedFile, mode='r')
            for line in decompressedFile.readlines():
                if len(line.split(" ")) == 2:
                    (t, temp) = line.split()
                    tempData[0] += [datetime.datetime.fromtimestamp(float(t))]
                    tempData[1] += [float(temp)]

            f = urllib2.urlopen(
                "http://www.apex-telescope.org/weather/Historical_weather/pressure/%s_pressure.log.gz"
                % datestr
            )
            compressedFile = StringIO.StringIO()
            compressedFile.write(f.read())
            compressedFile.seek(0)
            decompressedFile = gzip.GzipFile(fileobj=compressedFile, mode='r')
            for line in decompressedFile.readlines():
                if len(line.split(" ")) == 2:
                    (t, pres) = line.split()
                    presData[0] += [datetime.datetime.fromtimestamp(float(t))]
                    presData[1] += [float(pres)]

            f = urllib2.urlopen(
                "http://www.apex-telescope.org/weather/Historical_weather/humidity/%s_humidity.log.gz"
                % datestr
            )
            compressedFile = StringIO.StringIO()
            compressedFile.write(f.read())
            compressedFile.seek(0)
            decompressedFile = gzip.GzipFile(fileobj=compressedFile, mode='r')
            for line in decompressedFile.readlines():
                if len(line.split(" ")) == 2:
                    (t, rh) = line.split()
                    rhData[0] += [
                        datetime.datetime.fromtimestamp(float(t) + 18000.0)
                    ]  # 18000 seconds added to adjust for GMT time
                    rhData[1] += [float(rh)]
            for hour in range(24):
                for minute in range(60):
                    t = date + datetime.timedelta(hours=hour, minutes=minute)
                    wx['t'] += [t]
                    if t in tempData[0]:
                        wx['temp'] += [tempData[1][tempData[0].index(t)]]
                    else:
                        wx['temp'] += [float('NaN')]

                    if t in presData[0]:
                        wx['pres'] += [presData[1][presData[0].index(t)]]
                    else:
                        wx['pres'] += [float('NaN')]

                    if t in rhData[0]:
                        wx['rh'] += [rhData[1][rhData[0].index(t)]]
                    else:
                        wx['rh'] += [float('NaN')]
        wxArray = np.array([wx['t'], wx['temp'], wx['pres'], wx['rh']])
        wxArray = np.transpose(wxArray)
        print wxArray[0]
        with open(
            "/n/home10/anwang16/merra2_work/merra2_products/wx_data/APEX_wx_%s-%s.txt"
            % (year, month),
            'wb',
        ) as f:
            np.savetxt(f, wxArray, fmt="%s")
