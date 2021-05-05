#! /usr/bin/env python

from spt3g.util import files
from spt3g.std_processing import obsid_to_g3time
import numpy as np

if __name__ == '__main__':

    d = files.load_pickle('polarized_atmosphere_ratios_2019.pkl')

    data = [
        (d[obsid]['220'], obsid_to_g3time(obsid).isoformat())
        for obsid in sorted(d.keys())
    ]

    with open('ratios_dates.txt', 'w') as f:
        for ratio, obsid in data:
            f.write('%f %s\n' % (ratio, obsid))
