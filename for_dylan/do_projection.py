#! /usr/bin/env python

import numpy as np
from scipy import interpolate
import time
from matplotlib.pyplot import figure, imshow, gca, cm, savefig

def reproject(data):
    
    # Ok, so the idea here is that we want to always do math on arrays, not on individual floats.
    # The other idea is that we want to always use canned methods, and not rewrite things that 
    # some computer scientist has already optimized. If you gave me a set of 2d grid points in 
    # lon, lat this is how I would do the flat sky projection.
    
    # This should be whatever your sample points are, here I'm guessing at the bounds.
    # Both of these should be monotonically increasing, or else the scipy interpolate method later will error
    lat = np.linspace(-90., -85., num=data.shape[0])
    lon = np.linspace(-180., 180., num=data.shape[1])
      
    # Next, make the cartesian grid
    x = np.linspace(-5., 5., num=100)
    y = np.linspace(-5, 5., num=100)
    
    # Now let's turn these into 2d arrays
    xx, yy = np.meshgrid(x, y)
    
    # Ok, now here's the really important part. We're going to do the math ON ARRAYS, not on individual floats.
    # This allows us to call the numpy array math routines under the hood, which a) aren't slowed down by the python interpreter
    # and b) can use CPU vectorization (https://stackoverflow.com/questions/1422149/what-is-vectorization) to speed things up dramatically

    # Get the radial and angular coordinates, nothing that everything in these expressions is an array
    llat = -90. + np.sqrt(xx**2 + yy**2)
    # Make sure the -180 -> +180 versus 0 -> 360 convention is preserved here if you change the data bounds
    llon = np.arctan(xx/yy) * (180. / np.pi)

    # Now, let's use someone else's interpolation library to rebin the data
    # We have a rectangular grid in r, theta so use that fact to pick a faster interpolation function.
    # We'll use scipy.interpolate.RectBivariateSpline since it's going to be way faster than anything we'd write ourselves
    interpfn = interpolate.RectBivariateSpline(lat, lon, data)
    
    # Here the return array inherits the shape of rr and tt, which inherits the shape of xx and yy
    return interpfn(llat, llon, grid=False)

if __name__ == '__main__':

    # Run the interpolation
    begin=time.time()
    out = reproject(np.random.normal(loc=0., scale=1., size=(50, 180)))
    end=time.time()

    print('Elapsed time')
    print(end - begin)

    figure()
    imshow(out, interpolation='nearest', origin='lower', cmap=cm.gray)
    gca().grid(True)
    savefig('interpolate_test.png')


