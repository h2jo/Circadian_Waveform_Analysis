##### Codes written by Hang-Hyun Jo for the manuscript "Waveforms of molecular oscillations reveal circadian timekeeping mechanisms" by Jo et al. (submitted to Communications Biology)

import math

rt_base = 0.15

def get_cost():
    fn3 = "real_protein.txt"
    fn0 = "spline_protein.txt"
    fn1 = "real_rt.txt"
    fn2 = "real_degradation.txt"

    file3 = open(fn3, "r")
    file0 = open(fn0, "r")
    file1 = open(fn1, "r")
    file2 = open(fn2, "r")

    xt = {}
    for line in file3:
        items = [float(x) for x in line.split()]
        t = int(items[0] * 10)
        if t < 0 or t > 240:
            continue
        xt[t] = items[1]
    file3.close()

    Rt = {}
    for line in file0:
        items = [float(x) for x in line.split()]
        t = int(items[0] * 10)
        Rt[t] = items[3]
    file0.close()

    rt = {}
    for line in file1:
        items = [float(x) for x in line.split()]
        t = int(items[0] * 10)
        rt[t] = items[1] - rt_base
    file1.close()

    degt = {}
    for line in file2:
        items = [float(x) for x in line.split()]
        t = int(items[0] * 10)
        degt[t] = items[1]
    file1.close()

    r_min = min(max(Rt.values()), max(rt.values()))

    x_avg = sum(xt.values()) / len(xt.values())
    cost_global = sum(xt.values()) / len(xt.values()) * r_min
    cost_local = sum(degt.values()) / len(degt.values())

    print "<x>:", x_avg, "c_g:", cost_global, "c:", cost_local

########## main
if __name__ == "__main__":
    get_cost()
