##### Codes written by Hang-Hyun Jo for the manuscript "Waveforms of molecular oscillations reveal circadian timekeeping mechanisms" by Jo et al. (submitted to Communications Biology)

import numpy as np
import math
from scipy.interpolate import interp1d

##### constants

T = 23.4
T10 = 234

### protein
def read_data_protein(fn):
    file0 = open(fn, "r")

    ts = []
    xs = []

    num = 0. # for measuring h0
    avg = 0.
    for line in file0:
        items = line.split()
        t = int(float(items[0]) * 10. + 0.5) # integer
        x = float(items[1])
        x_19 = float(items[2])
        x_22 = float(items[3])
        x_25 = float(items[4])
        x_28 = float(items[5])
        x_30 = float(items[6])

        ts.append(t)
        xs.append(x)

        if t > 400 and t < 500: # for measuring h0
            avg += (x_19 + x_22 + x_25 + x_28 + x_30)
            num += 5
    file0.close()

    h0 = avg / num
    return ts, xs, h0

### moving average
def print_data_moving(ts, xs, h0, fn1):
    file0 = open(fn1, "w")

    window10 = 30 # window size = 3
    index_max = int(window10 * 1.5)
    num_t = len(ts)
    t_back = -1
    x_back = 0
    check = 0
    ts_R = []
    Rs = []
    for i, t in enumerate(ts):
        x0 = xs[i]
        avg = x0
        num = 1
        for j in range(1, index_max):
            if i + j >= num_t:
                break
            tj = ts[i + j]
            if tj < t + window10:
                xj = xs[i + j]
                avg += xj
                num += 1
        avg /= float(num)
        if check == 0:
            check = 1
        else:
            mid_t = (t + t_back + window10) / 2.
            x = avg - h0
            dx = (avg - x_back) / (t - t_back) * 10.
            R = -dx / x
            if R < 0:
                R = 0
            file0.write("%f\t%f\t%f\t%f\n" % (mid_t, x, dx, R))
            ts_R.append(mid_t)
            Rs.append(R)
        x_back = avg
        t_back = t
    file0.close()

    return ts_R, Rs

def print_data_moving_Rt(ts, xs, fn):
    file0 = open(fn, "w")

    window10 = 10
    index_max = int(window10 * 1.5)
    num_t = len(ts)
    t_back = -1
    x_back = 0
    for i, t in enumerate(ts):
        x0 = xs[i]
        avg = x0
        num = 1
        for j in range(1, index_max):
            if i + j >= num_t:
                break
            tj = ts[i + j]
            if tj < t + window10:
                xj = xs[i + j]
                avg += xj
                num += 1
        avg /= float(num)
        mid_t = (t + t_back + window10) / 2.
        file0.write("%f\t%f\n" % (mid_t, avg))

        x_back = avg
        t_back = t

    file0.close()

def read_data_mRNA(fn):
    file0 = open(fn, "r")

    ymin = 191.75 # 1
    ymax = 277.75 # 1000
    yreal = math.log(1000)

    x = []
    y = []
    for line in file0:
        items = line.split()
        t = float(items[0])
        x0 = float(items[1])
        x0 = math.exp((x0 - ymin) / (ymax - ymin) * yreal)
        x.append(t)
        y.append(x0)
    file0.close()

    return x, y

def print_data(x, y, fn):
    file0 = open(fn, "w")
    for index, x0 in enumerate(x):
        file0.write("%f\t%f\n" % (x0, y[index]))
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
            file0.write("%f\t%f\t%f\n" % (x0, y, dy))
    file0.close()

########## main
if __name__ == "__main__":
    fn0 = "raw_data.txt"
    fn1 = "real_protein.txt"
    fn2 = "moving_Rt.txt"

    ts, xs, h0 = read_data_protein(fn0)
    ts_R, Rs = print_data_moving(ts, xs, h0, fn1)
    print_data_moving_Rt(ts_R, Rs, fn2)
