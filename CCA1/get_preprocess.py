##### Codes written by Hang-Hyun Jo for the manuscript "Waveforms of molecular oscillations reveal circadian timekeeping mechanisms" by Jo et al. (submitted to Communications Biology)

import numpy as np
from scipy.interpolate import interp1d
import math

### from pixel to real number
def read_data(fn, prefix):
    file0 = open(fn, "r")

    if prefix[0] == 'p':
        ymin = 732.
        ymax = 207.
        yreal = 450.
    elif prefix[0] == 'r':
        ymin = 449.
        ymax = 133.
        yreal = 0.7
    elif prefix[0] == 'd':
        ymin = 713.
        ymax = 227.
        yreal = 100.

    x = []
    y = []
    for line in file0:
        items = line.split()
        t = float(items[0])
        x0 = float(items[1])
        x0 = (x0 - ymin) / (ymax - ymin) * yreal
        x.append(t)
        y.append(x0)
    file0.close()

    return x, y

def print_data(x, y, fn):
    file0 = open(fn, "w")
    for index, x0 in enumerate(x):
        file0.write("%f %f\n" % (x0, y[index]))
    file0.close()

def print_spline(x, f, fn):
    xmin = x.min()
    xmax = x.max()
    dx = 0.1

    file0 = open(fn, "w")
    for x0 in np.arange(xmin, xmax, dx):
        y = f(x0)
        if x0 >= xmin + dx:
            dy = (f(x0) - f(x0 - dx)) / dx
            if dy < 0:
                file0.write("%f %f %f %f\n" % (x0, y, dy, -dy/y))
            else:
                file0.write("%f %f %f 0\n" % (x0, y, dy))
    file0.close()

########## main
if __name__ == "__main__":
    prefixes = ["protein", "rt", "degradation"]
    for prefix in prefixes:
        fn0 = "pixel_" + prefix + ".txt"
        fn1 = "real_" + prefix + ".txt"
        if prefix[0] == 'p':
            fn2 = "spline_" + prefix + ".txt"

        x, y = read_data(fn0, prefix)
        print_data(x, y, fn1)

        if prefix[0] == 'p':
            x0 = np.asarray(x)
            y0 = np.asarray(y)
            f = interp1d(x0, y0, kind='cubic')
            print_spline(x0, f, fn2)
