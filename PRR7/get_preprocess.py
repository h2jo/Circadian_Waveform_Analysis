##### Codes written by Hang-Hyun Jo for the manuscript "Waveforms of molecular oscillations reveal circadian timekeeping mechanisms" by Jo et al. (submitted to Communications Biology)

import numpy as np
import math
from scipy.interpolate import interp1d

##### constants and conditions

p_ymin = 122.5
p_ymax = 187.5
p_yreal = 1.2
ymax_inverse = 0.0055

m_ymin = 191.75 # 1
m_ymax = 277.75 # 1000
m_yreal = math.log(1000)

### protein: from pixel to real number
def read_data_protein(fn):
    file0 = open(fn, "r")

    x = []
    y = []
    for line in file0:
        items = line.split()
        t = float(items[0])
        x0 = float(items[1])
        x0 = (x0 - p_ymin) / (p_ymax - p_ymin) * p_yreal
        x.append(t)
        y.append(x0)
    file0.close()

    return x, y

### mRNA: from pixel to real number
def read_data_mRNA(fn):
    file0 = open(fn, "r")

    x = []
    y = []
    for line in file0:
        items = line.split()
        t = float(items[0])
        x0 = float(items[1])
        x0 = math.exp((x0 - m_ymin) / (m_ymax - m_ymin) * m_yreal)
        x.append(t)
        y.append(x0)
    file0.close()

    return x, y

### print data
def print_data(x, y, fn, common):
    file0 = open(fn, "w")
    if common[0] == "m":
        for index, x0 in enumerate(x):
            file0.write("%f\t%f\n" % (x0, y[index] * ymax_inverse))
    else:
        for index, x0 in enumerate(x):
            file0.write("%f\t%f\n" % (x0, y[index]))
    file0.close()

### print splined curve
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
                file0.write("%f\t%f\t%f\t%f\n" % (x0, y, dy, -dy/y))
            else:
                file0.write("%f\t%f\t%f\t0\n" % (x0, y, dy))
    file0.close()

########## main
if __name__ == "__main__":
    commons = ["protein_origin", "protein_smooth", "mRNA"]
    for index, common in enumerate(commons):
        fn0 = "pixel_" + common + ".txt"
        fn1 = "real_" + common + ".txt"
        fn2 = "spline_" + common + ".txt"

        if index <= 1:
            x, y = read_data_protein(fn0)
        elif index == 2:
            x, y = read_data_mRNA(fn0)
        print_data(x, y, fn1, common)

        x0 = np.asarray(x)
        y0 = np.asarray(y)
        f = interp1d(x0, y0, kind='cubic')
        print_spline(x0, f, fn2)
