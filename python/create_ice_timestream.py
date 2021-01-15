#! /usr/bin/env python

import argparse
import datetime
import os
import merra2Player as m2p
import netCDF4
import numpy as np
import pickle

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--location",
        default='SouthPole',
        help="Sites are: SouthPole, ChajnantorPlateau, CerroChajnantor, MaunaKea, Summit, Qubic. Default: SouthPole",
    )

    parser.add_argument(
        "--date",
        default=["20190101", "20201129"],
        nargs=2,
        help="Date range in the format YYYYMMDD, YYYYMMDD",
    )

    # parser.add_argument(
    #    "-g",
    #    dest="groundData",
    #    default=2,
    #    type=int,
    #    help="-g. Type of surface level data. 1: historical measured data when available. 2(default): MERRA2 surface data. 3: Extrapolation from MERRA2 3d data.",
    # )

    # parser.add_argument(
    #    "-c",
    #    dest="cloud",
    #    action="store_true",
    #    default=False,
    #    help="-c. adds ice and liquid water layers in the am profile. Default (false) is for am profile to only contain water-vapor data.",
    # )

    # parser.add_argument(
    #    "-s",
    #    dest="saveResults",
    #    action="store_false",
    #    default=True,
    #    help="-s. Default=True, Forces the output of a run to be saved in a csv and pickle file",
    # )

    args = parser.parse_args()

    import merra2Player as m2p

    m = m2p.merra2Player()

    # This should load the cloud data (?)
    siteopt = {}
    siteopt['type'] = 'presel'
    siteopt['name'] = args.location
    siteopt['cldtype'] = 'icewater'
    siteopt['groundData'] = 2
    siteopt = m.defineSite(site=siteopt)

    # Simplify the datelist interface
    dateopt = {}
    dateopt['type'] = 'datelist'
    dateopt['start'] = args.date[0]
    dateopt['end'] = args.date[1]
    dateopt['single'] = False

    # Links and directory names have to be defined first as they are used throughout
    print(dateopt, siteopt)
    m.checkLinks()

    # Now let's load the MERRA2 data for the correct date range
    datelist = m.retrieve_merra2_data_for_dateRange(dateopt['start'], dateopt['end'])

    # Here is where we will store the timestreams
    merradata = {}

    # merradata['QI'] = np.zeros((8 * len(datelist), 42, 3, 3), dtype=float)
    for key in ['TQI', 'TQL', 'TQV']:
        merradata[key] = np.zeros((24 * len(datelist), 3, 3), dtype=float)

    merradata['start'] = dateopt['start']
    merradata['end'] = dateopt['end']

    # Let's chuck all of this into a dictionary
    for i, date in enumerate(datelist):

        print(date)

        # Let's load the netDCF4 filenames
        url, merraFilename_m = m.get_url_for_date(date.strftime('%Y%m%d'), 'multiLevel')
        url, merraFilename_s = m.get_url_for_date(
            date.strftime('%Y%m%d'), 'singleLevel'
        )
        del url

        # Let's load in the CDF4 data
        data_s = netCDF4.Dataset(os.path.join(m.merraDir, merraFilename_s))
        data_m = netCDF4.Dataset(os.path.join(m.merraDir, merraFilename_m))
        print(data_s.variables.keys())
        print(data_m.variables.keys())

        # Stick it in the TODs
        # merradata['QI'][i * 8 : (i + 1) * 8, :, :] = np.array(data_m.variables['QI'])
        for key in ['TQI', 'TQL', 'TQV']:
            merradata[key][i * 24 : (i + 1) * 24, :, :] = np.array(
                data_s.variables[key]
            )

    # Now save it out
    with open('merra2_timestream.pkl', 'wb') as handle:
        pickle.dump(merradata, handle, protocol=pickle.HIGHEST_PROTOCOL)
