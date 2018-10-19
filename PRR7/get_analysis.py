##### Codes written by Hang-Hyun Jo for the manuscript "Waveforms of molecular oscillations reveal circadian timekeeping mechanisms" by Jo et al. (submitted to Communications Biology)

import math

##### constants and conditions

rt_4h = 0.0891687
rt_12h = 0.341079
#rt_18h = 0.446286

r_min_origin = 0.876871
r_min_smooth = 0.687622

rt_max_origin = 1.019178
rt_max_smooth = 0.740464

type_xt = ["_origin", "_smooth"] # original or smoothened x(t)
type_rt = ["_corrected", "_uncorrected"] # corrected or uncorrected r(t)
type_r18h = ["_r18avg", "_r18low", "_r18high"] # experimental r(18h)
r18hs = [0.446286, 0.342691, 0.553561]
    # the value of r(18h): average, average - std, average + std
type_kt = ["_kc", "_kt"] # constant or time-varying k(t)

##### functions

def get_kt_origin(t): # fitted to sine wave
    kmin = 0.00117188958654226
    kta = -0.0171683
    ktb = -0.280267
    ktc = 0.0130855
    kt = kta * math.sin(t * 3.141592 / 12. + ktb) + ktc
    if kt < kmin:
        kt = kmin
    return kt

def get_cost(index_xt, index_rt, index_r18h, index_kt):
    postfix_xt = type_xt[index_xt]
    postfix_rt = type_rt[index_rt]
    postfix_r18h = type_r18h[index_r18h]
    rt_18h = r18hs[index_r18h]
    postfix_kt = type_kt[index_kt]

    fn0 = "spline_protein" + postfix_xt + ".txt"
    fn1 = "real_mRNA.txt"
    fn2 = "rt" + postfix_xt + postfix_rt + postfix_r18h + postfix_kt + ".txt"

    file0 = open(fn0, "r")
    file1 = open(fn1, "r")
    file2 = open(fn2, "w")

    xt = {}
    dxt = {}
    for line in file0:
        items = [float(x) for x in line.split()]
        t = int(items[0] * 10)
        xt[t] = items[1]
        dxt[t] = items[2]
    file0.close()

    yt = {}
    for line in file1:
        items = [float(x) for x in line.split()]
        t = int(items[0] * 10)
        yt[t] = items[1]
    file1.close()

    k_4h = (dxt[40] + rt_4h * xt[40]) / yt[40]
    k_12h = (dxt[120] + rt_12h * xt[120]) / yt[120]
    k_18h = (dxt[180] + rt_18h * xt[180]) / yt[180]
    k_const = (k_4h + k_12h + k_18h) / 3. # when k(t) is constant

    if index_xt == 0: # original x(t)
        r_min = r_min_origin
        rt_max = rt_max_origin
    elif index_xt == 1: # smoothened x(t)
        r_min = r_min_smooth
        rt_max = rt_max_smooth

    cost_global = 0.
    cost_local = 0.
    num = 0.
    for t in sorted(yt.keys()):
        if t == -120 or t == 360:
            continue

        if index_kt == 0: # const k
            kk = k_const
        elif index_kt == 1: # time-varying k(t), only for original x(t)
            kk = get_kt_origin(t * 0.1)

        rt = (kk * yt[t] - dxt[t]) / xt[t]

        # common procedure
        if t == 0 or t == 40 or t == 240:
            rt = rt_4h
        if t == 120:
            rt = rt_12h
        if t == 180:
            rt = rt_18h
        if rt < rt_4h:
            rt = rt_4h
        if rt > rt_max:
            rt = rt_max

        # correction for r(t): interpolation and extrapolation
        if index_rt == 0:
            if t == 20:
                rt = rt_4h
            if t >= 60 and t <= 100:
                rt = rt_4h + (rt_12h - rt_4h) / 80. * (float(t) - 40.)

        file2.write("%f\t%f\n" % (t * 0.1, rt)) # print r(t)

        # calculate costs
        if t >= 0 and t < 240:
            cost_local += rt * xt[t]
            cost_global += xt[t]
            num += 1.
    file2.close()

    cost_local = cost_local / num
    cost_global = cost_global / num * r_min

    print(fn2)
    print "c_g:", cost_global, "c:", cost_local


########## main
if __name__ == "__main__":
    # usage: get_cost(index_xt, index_rt, index_r18h, index_kt)
    # 1st argument: index of type_xt
    # 2nd argument: index of type_rt
    # 3rd argument: index of type_r18h
    # 4th argument: index of type_kt
    get_cost(0, 0, 0, 0)
    get_cost(0, 1, 0, 0)
    get_cost(1, 0, 0, 0)
    get_cost(1, 1, 0, 0)
    get_cost(0, 0, 1, 0)
    get_cost(0, 0, 2, 0)
    get_cost(0, 0, 0, 1)
