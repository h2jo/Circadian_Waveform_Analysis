##### Codes written by Hang-Hyun Jo for the manuscript "Waveforms of molecular oscillations reveal circadian timekeeping mechanisms" by Jo et al. (submitted to Communications Biology)

import numpy as np
import math
from scipy.interpolate import interp1d

##### constants

r_min = 0.473573
T = 23.4
h0 = 21.86

def read_data_halflife(fn0, fn1):
    file0 = open(fn0, "r")
    file1 = open(fn1, "w")

    x = []
    y = []
    for line in file0:
        items = line.split()
        t = float(items[0])
        halflife = float(items[1])
        r = math.log(2) / halflife
        file1.write("%f\t%f\n" % (t, r))
        x.append(t)
        y.append(r)
    file0.close()
    file1.close()

    return x, y

def print_data(x, y, fn):
    file0 = open(fn, "w")
    for index, x0 in enumerate(x):
        file0.write("%f\t%f\n" % (x0, y[index]))
    file0.close()

def print_spline_rt(x, y, f, fn):
    xmin = x.min()
    xmax = x.max()
    dx = 0.1

    T = 234
    t_begin = 330 - T
    t_end = 330

    r0 = y.min()
    r1 = y.max()

    file0 = open(fn, "w")
    for t in range(t_begin, 190):
        if t <= 149:
            rd = r_min
        elif t <= 190:
            rd = r0
        file0.write("%f\t%f\n" % (t * 0.1, rd))

    for x0 in np.arange(xmin, xmax, dx):
        y0 = f(x0)
        file0.write("%f\t%f\n" % (x0, y0))

    for t in range(300, t_end):
        rd = r1
        file0.write("%f\t%f\n" % (t * 0.1, rd))
    file0.close()

def get_cost(fn0, fn1):
    file0 = open(fn0, "r")
    file1 = open(fn1, "r")

    rs = {}
    for line in file1:
        items = [float(x) for x in line.split()]
        t = int(items[0] * 10. + 0.5)
        r = items[1]
        rs[t] = r
    file1.close()

    t_begin = int((33.-T) * 10.)
    t_end = 330
    cost_local = 0.
    cost_global = 0.
    num = 0
    x_max = 121.563552
    for line in file0:
        items = [float(x) for x in line.split()]
        t = int(items[0] + 0.5)
        x = items[1] / x_max
        dx = items[2] / x_max
        if t >= t_begin and t < t_end:
            cost_global += x
            cost_local += rs[t] * x
            num += 1
    file0.close()

    cost_global = cost_global / float(num) * r_min
    cost_local = cost_local / float(num)
    print "c_g:", cost_global, "c:", cost_local

########## main
if __name__ == "__main__":
    fn0 = "halflifet_empirics.txt"
    fn1 = "rt_empirics.txt"
    fn2 = "rt.txt"
    fn3 = "real_protein.txt"

    x, y = read_data_halflife(fn0, fn1)
    x = np.asarray(x)
    y = np.asarray(y)
    f = interp1d(x, y, kind='linear')

    print_spline_rt(x, y, f, fn2)

    get_cost(fn3, fn2)
