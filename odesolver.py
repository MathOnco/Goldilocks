import numpy as np

# A dormand prince solver

# returns one step of the Dormand Prince solver analysis

def _dopri45(dydt, h, y0, t0, t1, e_abs):


    c2 = 1/5.0
    c3 = 3/10.0
    c4 = 4/5.0
    c5 = 8/9.0
    c6 = 1.0
    c7 = 1.0

    a21 = 1/5.0

    a31 = 3/40.0
    a32 = 9/40.0

    a41 = 44/45.0
    a42 = -56/15.0
    a43 = 32/9.0

    a51 = 19372/6561.0
    a52 = -25360/2187.0
    a53 = 64448/6561.0
    a54 = -212/729

    a61 = 9017/3168.0
    a62 = -355/33.0
    a63 = 46732/5247.0
    a64 = 49/176.0
    a65 = -5103/18656.0

    a71 = 35/384.0
    a72 = 0
    a73 = 500/1113.0
    a74 = 125/192.0
    a75 = -2187/6784.0
    a76 = 11/84.0

    b11 = a71
    b12 = a72
    b13 = a73
    b14 = a74
    b15 = a75
    b16 = a76
    b17 = 0.0

    b21 = 5179/57600.0
    b22 = 0
    b23 = 7571/16695.0
    b24 = 393/640.0
    b25 = -92097/339200
    b26 = 187/2100.0
    b27 = 1/40.0

    #print(y0)

    k1 = dydt(t0, y0)

    k2 = dydt(t0 + c2 * h, y0 + h * (np.multiply(a21,k1)))

    k3 = dydt(t0 + c3 * h, y0 + h * (np.multiply(a31, k1) + np.multiply(a32, k2)))

    k4 = dydt(t0 + c4 * h, y0 + h * (np.multiply(a41, k1) + np.multiply(a42, k2) + np.multiply(a43, k3)))

    k5 = dydt(t0 + c5 * h, y0 + h * (np.multiply(a51, k1) + np.multiply(a52, k2) + np.multiply(a53, k3) + np.multiply(a54, k4)))

    k6 = dydt(t0 + c6 * h, y0 + h * (np.multiply(a61, k1) + np.multiply(a62, k2) + np.multiply(a63, k3) + np.multiply(a64, k4) + np.multiply(a65, k5)))

    last = y0 + h * (np.multiply(a71, k1) + np.multiply(a72, k2) + np.multiply(a73, k3) + np.multiply(a74, k4) + np.multiply(a75, k5) + np.multiply(a76, k6))

    k7 = dydt(t0 + c7 * h, last)

    f1 = last

    f2 = np.multiply(b21, k1) + np.multiply(b22, k2) + np.multiply(b23, k3) + np.multiply(b24, k4) + np.multiply(b25, k5) + np.multiply(b26, k6) + np.multiply(b27, k7)

    error = f1 - f2
    error = np.max(error)
    # Calculate s
    s = ((h * e_abs)/(2 * (t1 - t0) * np.abs(error)))**0.25
    return f1, 1.5

def dopri45(dydt, dt, y0, t0, t1, e_abs=0.001, maxsteps = 100):
    h = dt

    new_row = np.hstack((t0, y0))

    results = new_row

    nsteps = 0
    while t0 < t1:
        success = False
        y0, s = _dopri45(dydt, h, y0, t0, t0 + h, e_abs = e_abs)
        if s >= 2.0:
            h = 2 * h
            success = True
        elif s >= 1:
            success = True
        elif nsteps >= maxsteps:
            success = True
        else:
            success = False
            h = h/2.0
            nsteps += 1

        if success:

            t0 += h
            new_row = np.hstack([[t0,],y0])
            results = np.vstack([results, new_row])
            nsteps = 0

    return np.array(results)
