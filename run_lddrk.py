
import argparse

import numpy
import matplotlib.pyplot as P

import ricker

LDDRK_ALPHA = [0.0, -0.737101392796, -1.634740794341, -0.744739003780, -1.469897351522, -2.813971388035]
LDDRK_BETA = [0.032918605146, 0.823256998200, 0.381530948900, 0.200092213184, 1.718581042715, 0.27]
LDDRK_C = [0.0, 0.032918605146, 0.249351723343, 0.466911705055, 0.582030414044, 0.847252983783]
LDDRK_STAGES = len(LDDRK_ALPHA)

def d2udt2(Minv, C, K, stf, source, t, s, v):

    a = -K.dot(s) - C.dot(v)
    a[source] += stf(t)
    a = numpy.multiply(Minv, a)

    return a

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input', type = str, required = True, help = 'Input prefix')
    parser.add_argument('-o', '--output', type = str, required = True, help = 'Output file')

    parser.add_argument('-s', '--sourcex', type = float, default = 0.0, help = 'Source x')
    parser.add_argument('-S', '--sourcey', type = float, default = 0.0, help = 'Source y')
    parser.add_argument('-v', '--vertical', action = 'store_true', default = False, help = 'Vertical source')

    parser.add_argument('-r', '--receiverx', type = float, default = 1.0, help = 'Receiver x')
    parser.add_argument('-R', '--receivery', type = float, default = 1.0, help = 'Receiver y')

    #
    # Source
    #
    parser.add_argument('-t', '--t0', type = float, default = 1.0, help = 'Start time of source')
    parser.add_argument('-A', '--amplitude', type = float, default = 1.0, help = 'Amplitude of source')
    parser.add_argument('-f', '--frequency', type = float, default = 0.4, help = 'Dominant frequency of source')

    #
    # Simulation
    #
    parser.add_argument('-d', '--dt', type = float, default = 1.0e-3, help = 'Time step')
    parser.add_argument('-N', '--nsteps', type = int, default = 1000, help = 'No. steps')

    
    
    args = parser.parse_args()

    #
    # Load Matrices (M is a vector)
    #
    M = numpy.load('%sM.npy' % args.input)
    K = numpy.load('%sK.npy' % args.input)
    C = numpy.load('%sC.npy' % args.input)
    points = numpy.load('%sPoints.npy' % args.input)

    Minv = 1.0/M

    #
    # Locate source and receiver
    #
    nglob, _ = points.shape
    source_ix = numpy.argmin((points[:,0] - args.sourcex)**2 + (points[:,1] - args.sourcey)**2)
    source_iy = source_ix + nglob

    source = source_ix
    if args.vertical:
        source = source_iy

    print('Located source  : %4d %4d (requested: %10.6f %10.6f, actual %10.6f %10.6f)' % (source_ix, source_iy,
                                                                                          args.sourcex, args.sourcey,
                                                                                          points[source_ix,0], points[source_ix,1]))

    receiver_ix = numpy.argmin((points[:,0] - args.receiverx)**2 + (points[:,1] - args.receivery)**2)
    receiver_iy = receiver_ix + nglob
    print('Located receiver: %4d %4d (requested: %10.6f %10.6f, actual %10.6f %10.6f)' % (receiver_ix, receiver_iy,
                                                                                          args.receiverx, args.receivery,
                                                                                          points[receiver_ix,0], points[receiver_ix,1]))

    #
    # State variables
    #
    seismogram = numpy.zeros((args.nsteps, 4))

    s = numpy.zeros((2*nglob,))
    v = numpy.zeros((2*nglob,))
    a = numpy.zeros((2*nglob,))

    s_lddrk = numpy.zeros((2*nglob,))
    v_lddrk = numpy.zeros((2*nglob,))

    #
    # Source time function
    #
    stf = ricker.mkricker(args.amplitude, args.frequency, args.t0)
        
    #
    # Loop
    #
    
    for i in range(args.nsteps):

        t = float(i + 1) * args.dt

        for j in range(LDDRK_STAGES):

            a = d2udt2(Minv, C, K, stf, source, t + args.dt*LDDRK_C[j], s, v)

            v_lddrk = LDDRK_ALPHA[j]*v_lddrk + args.dt*a
            s_lddrk = LDDRK_ALPHA[j]*s_lddrk + args.dt*v

            v += LDDRK_BETA[j]*v_lddrk
            s += LDDRK_BETA[j]*s_lddrk

        seismogram[i, 0] = t
        seismogram[i, 1] = stf(seismogram[i, 0])
        seismogram[i, 2] = s[receiver_ix]
        seismogram[i, 3] = s[receiver_iy]


    numpy.save(args.output, seismogram)
    
    fig, ax = P.subplots(2, 1)

    ax[0].plot(seismogram[:, 0], seismogram[:, 2])
    ax[1].plot(seismogram[:, 0], seismogram[:, 3])

    P.show()
        
    
    
    
    
    
