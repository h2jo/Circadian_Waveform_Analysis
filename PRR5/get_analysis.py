##### Codes written by Hang-Hyun Jo for the manuscript "Waveforms of molecular oscillations reveal circadian timekeeping mechanisms" by Jo et al. (submitted to Communications Biology)

import math

##### constants and conditions

rt_12h = 0.232185
rt_19h = 0.414324

r_min_origin = 1.689532
r_min_smooth = 0.546538

rt_max_origin = 2.467566
rt_max_smooth = 0.854944

type_xt = ["_origin", "_smooth"] # original or smoothened x(t)
type_rt = ["_corrected", "_uncorrected"] # corrected or uncorrected r(t)
type_kt = ["_kc", "_kt"] # constant or time-varying k(t)

##### functions

def get_kt_origin(t, phase):
    kmin = 0.003542206643591266
    if phase == 0:
        kta = 0.00744587
        ktb = 0.0109881
    elif phase == 3:
        kta = -0.00589984
        ktb = 0.00309332
    elif phase == 6:
        kta = -0.00744588
        ktb = 0.0035422
    elif phase == 9:
        kta = -0.00526503
        ktb = 0.00726514
    elif phase == 12:
        kta = -0.00744586
        ktb = 0.0109881
    elif phase == 15:
        kta = -0.00709143
        ktb = 0.0122795
    elif phase == 18:
        kta = 0.00744589
        ktb = 0.00354219
    elif phase == 21:
        kta = 0.00526503
        ktb = 0.00726514
    kt = kta * math.sin((t + phase) * 3.141592 / 12.) + ktb
    if kt < kmin:
        kt = kmin
    return kt

def get_cost(index_xt, index_rt, index_kt, phase):
    postfix_xt = type_xt[index_xt]
    postfix_rt = type_rt[index_rt]
    postfix_kt = type_kt[index_kt]

    fn0 = "spline_protein" + postfix_xt + ".txt"
    fn1 = "real_mRNA.txt"
    if index_kt == 0:
        fn2 = "rt" + postfix_xt + postfix_rt + postfix_kt + ".txt"
    else:
        fn2 = "rt" + postfix_xt + postfix_rt + postfix_kt + str(phase) + ".txt"

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

    k_12h = (dxt[120] + rt_12h * xt[120]) / yt[120]
    k_19h = (dxt[180] + rt_19h * xt[180]) / yt[180]
    k_const = (k_19h + k_12h) * 0.5 # when k(t) is constant

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
            kk = get_kt_origin(t * 0.1, phase)

        rt = (kk * yt[t] - dxt[t]) / xt[t]

        # common procedure
        if t == 0 or t == 120 or t == 240:
            rt = rt_12h
        if t == 180:
            rt = rt_19h
        if rt < rt_12h:
            rt = rt_12h
        if rt > r_min:
            rt = r_min

        # correction for r(t): interpolation and extrapolation
        if index_rt == 0:
            if t >= 0 and t <= 120:
                rt = rt_12h

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
    # usage: get_cost(index_xt, index_rt, index_kt, phase)
    # 1st argument: index of type_xt
    # 2nd argument: index of type_rt
    # 3rd argument: index of type_kt
    # 4th argument: phase in the case with k(t)
    get_cost(0, 0, 0, 0)
    get_cost(0, 1, 0, 0)
    get_cost(1, 0, 0, 0)
    get_cost(1, 1, 0, 0)
    for phase in [0, 3, 6, 9, 12, 15, 18, 21]:
        get_cost(0, 0, 1, phase)
