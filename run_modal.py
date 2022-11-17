
import argparse

import numpy
import matplotlib.pyplot as P

import ricker

def fftconvolve(gf, stf):

    N = gf.size
    Ap = numpy.zeros((2*N,))
    Bp = numpy.zeros((2*N,))

    o = N//2
    Ap[o:o + N] = gf
    Bp[o:o + N] = stf

    return numpy.fft.fftshift(numpy.fft.irfft(numpy.fft.rfft(Ap) * numpy.fft.rfft(Bp)))[:N]

def replace_rb_modes(indices, M, C, K, points):

    if len(indices) != 3:
        raise Exception('Expected 3 RB modes')
    
    N, _ = K.shape
    nglob = N//2
    phi = numpy.zeros((N, 3))

    print('RB replace')
    print(phi.shape)
    print(nglob)

    phi[:nglob, 0] = 1.0
    phi[nglob:, 1] = 1.0

    xmin = numpy.min(points[:, 0])
    xmax = numpy.max(points[:, 0])
    ymin = numpy.min(points[:, 1])
    ymax = numpy.max(points[:, 1])

    cx = (xmin + xmax)/2.0
    cy = (ymin + ymax)/2.0

    rx = points[:, 0] - cx
    ry = points[:, 1] - cy

    phi[:nglob, 2] = ry
    phi[nglob:, 2] = -rx

    for i in range(3):
        print(i, numpy.linalg.norm(K.dot(phi[:, i])))

    for i in range(2):
        for j in range(i + 1, 3):
            print(i, j, phi[:, i].dot(phi[:, j]))
        
#    fig, ax = P.subplots()
#
#    ax.quiver(points[:, 0], points[:, 1], phi[:nglob, 2], phi[nglob:, 2])
#
#    P.show()

    Chat = phi.T.dot(C.dot(phi))

    alpha_left, s, alpha_right = numpy.linalg.svd(Chat, full_matrices = True, compute_uv = True)
    alpha_right = alpha_right.T

    for i in range(3):
        if (s[i] < 0.0):
            n = numpy.sqrt(-s[i])
            alpha_left[:,i] /= n
            alpha_right[:,i] /= n
        else:
            n = numpy.sqrt(s[i])
            alpha_left[:,i] /= -n
            alpha_right[:,i] /= n
            
    
    phi_left = numpy.zeros((2*N, 3))
    phi_right = numpy.zeros((2*N, 3))

    phi_left[:N,:] = phi.dot(alpha_left)
    phi_left[N:,:] = C.dot(phi.dot(alpha_left))

    phi_right[N:,:] = phi.dot(alpha_right)

    B = -numpy.ones((2*N,))
    B[:N] = -M
    for i in range(3):
        n = phi_left[:,i].dot(numpy.multiply(B, phi_right[:,i]))
        if (n < 0.0):
            n = numpy.sqrt(-n)
            phi_left[:,i] /= -n
            phi_right[:,i] /= n
        else:
            n = numpy.sqrt(n)
            phi_left[:,i] /= n
            phi_right[:,i] /= n

    return phi_left, phi_right

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input', type = str, required = True, help = 'Input prefix (matrices)')
    parser.add_argument('-e', '--eigen', type = str, required = True, help = 'Input prefix (eigen soln.)')

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

    #
    # Comparison
    #
    parser.add_argument('-c', '--compare', type = str, default = None, help = 'Compare seismogram (npy)')

    #
    # Rigid body modes options
    #
    parser.add_argument('--replace-rb', action = 'store_true', default = False, help = 'Replace RB modes with analytical')

    parser.add_argument('-o', '--output', type = str, default = None, help = 'Output seismogram')

    #
    # Truncating eigen solution
    #
    parser.add_argument('-F', '--fmax', type = float, default = 1.0e9, help = 'Max frequency of eigensolution')
    parser.add_argument('-T', '--rmax', type = float, default = 1.0e9, help = 'Max real component of eigensolution')

    #
    # Tapering Greens
    #
    parser.add_argument('--taper', type = float, default = 0.0, help = 'Taper cutoff for Greens tapering')
    parser.add_argument('--taper-frequency', type = float, default = 1.0, help = 'Taper frequency for Greens tapering')
    
    #
    # Plots
    #
    parser.add_argument('--show', action = 'store_true', default = False, help = 'Plot GF and seismograms')
    
    args = parser.parse_args()

    points = numpy.load('%sPoints.npy' % args.input)

    eigval = numpy.load('%s-eigval.npy' % args.eigen)
    eigleft = numpy.load('%s-eigleft.npy' % args.eigen)
    eigright = numpy.load('%s-eigright.npy' % args.eigen)
    eigscale = numpy.load('%s-eigscale.npy' % args.eigen)
    gamma = eigscale[0]
    delta = eigscale[1]
    print('Gamma: ', gamma)
    print('Delta: ', delta)

    M = numpy.load('%sM.npy' % args.input)
    M *= gamma*gamma*delta
    
    size = M.size
    B = numpy.zeros((2*size, 2*size))
    B[:size, :size] = -numpy.diag(M)
    B[size:, size:] = -numpy.eye(size)

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
    # Output
    #
    seismogram = numpy.zeros((args.nsteps, 6))
    gx = numpy.zeros((args.nsteps,), dtype = 'complex')
    gy = numpy.zeros((args.nsteps,), dtype = 'complex')

    RECEIVERX = receiver_ix + 2*nglob
    RECEIVERY = receiver_iy + 2*nglob

    print('SOURCE:', source)
    print('RECV X:', RECEIVERX)
    print('RECV Y:', RECEIVERY)

    #
    # Source time function
    #
    stf = ricker.mkricker(args.amplitude, args.frequency, args.t0)

    seismogram[:,0] = (1.0 + numpy.arange(args.nsteps, dtype = 'float')) * args.dt
    seismogram[:,1] = stf(seismogram[:,0])

    #
    # Rigid body modes are the 3 nearest zero
    #
    rb_indices = numpy.argpartition(numpy.abs(eigval), 3)[:3]
    if args.replace_rb:
        K = numpy.load('%sK.npy' % args.input)
        C = numpy.load('%sC.npy' % args.input)

        C *= gamma*delta
        K *= delta

        L, R = replace_rb_modes(rb_indices, M, C, K, points)
        N = M.size
        B = -numpy.ones((2*N,))
        B[:N] = -M
        
        eigleft[:, rb_indices] = L
        eigright[:, rb_indices] = R

    #
    # Construct the greens function
    #
    N = eigval.size
    eigsused = 0
    ceigsused = 0
    for i in range(N):

        if numpy.imag(eigval[i]) != 0.0:

            f = numpy.abs(gamma*numpy.imag(eigval[i]))/(2.0 * numpy.pi)
            r = numpy.abs(numpy.real(gamma * eigval[i]))
            
            if f < args.fmax and r < args.rmax:
                scale = gamma * delta * numpy.conj(eigleft[source, i]) * eigright[RECEIVERX, i]
                gx -= scale*numpy.exp(gamma * eigval[i]*seismogram[:,0])

                scale = gamma * delta * numpy.conj(eigleft[source, i]) * eigright[RECEIVERY, i]
                gy -= scale*numpy.exp(gamma * eigval[i]*seismogram[:,0])

                eigsused += 1
                ceigsused += 1

        else:
            if i in rb_indices:
                scale = gamma * delta * numpy.real(eigleft[source, i]) * numpy.real(eigright[RECEIVERX, i])
                gx -= scale
                
                scale = gamma * delta * numpy.real(eigleft[source, i]) * numpy.real(eigright[RECEIVERY, i])
                gy -= scale
                
                eigsused += 1
                
            else:
                r = numpy.abs(numpy.real(gamma * eigval[i]))
                if r < args.rmax:
                    scale = gamma * delta * numpy.real(eigleft[source, i]) * eigright[RECEIVERX, i]
                    gx -= scale*numpy.exp(gamma * numpy.real(eigval[i]) * seismogram[:,0])
                
                    scale = gamma * delta * numpy.real(eigleft[source, i]) * eigright[RECEIVERY, i]
                    gy -= scale*numpy.exp(gamma * numpy.real(eigval[i]) * seismogram[:,0])

                    eigsused += 1 

    print('%d of %d eigen solutions used' % (eigsused, N))
    print('%d complex eigs, %d real' % (ceigsused, eigsused - ceigsused))
    print('%d half side eigs' % (ceigsused/2 + eigsused - ceigsused))
    
    if args.show:
        fig, ax = P.subplots(2, 2)

    
        ax[0][0].plot(seismogram[:, 0], numpy.real(gx))
        ax[1][0].plot(seismogram[:, 0], numpy.real(gy))

        ax[0][1].plot(seismogram[:, 0], numpy.imag(gx))
        ax[1][1].plot(seismogram[:, 0], numpy.imag(gy))

    if args.taper > 0.0:
        taper = numpy.ones((gx.size,))
        taper_width = 0.5/args.taper_frequency

        t = seismogram[:, 0]
        t0 = (args.taper - taper_width)
        t1 = args.taper
        tc = (t0 + t1)/2.0
        indices = numpy.where((t >= t0) & (t <= t1))[0]
        taper[indices] = (numpy.sin(2.0 * numpy.pi * args.taper_frequency * (t[indices] - tc)) + 1.0)/2.0
        indices = numpy.where((t < t0))[0]
        taper[indices] = 0.0

        gx = numpy.multiply(taper, gx)
        gy = numpy.multiply(taper, gy)

        if args.show:
            ax[0][0].plot(seismogram[:, 0], numpy.real(gx))
            ax[0][0].plot(seismogram[:, 0], taper * numpy.max(numpy.real(gx)))
            ax[1][0].plot(seismogram[:, 0], numpy.real(gy))

            ax[0][1].plot(seismogram[:, 0], numpy.imag(gx))
            ax[1][1].plot(seismogram[:, 0], numpy.imag(gy))
            
        
    rx = fftconvolve(numpy.real(gx), seismogram[:,1]) * args.dt
    ry = fftconvolve(numpy.real(gy), seismogram[:,1]) * args.dt

    if args.show:
        fig, ax = P.subplots(2, 1)

        ax[0].plot(seismogram[:,0], rx)
        ax[1].plot(seismogram[:,0], ry)

        if not (args.compare is None):
            t = numpy.load(args.compare)

            ax[0].plot(t[:,0], t[:,2])
            ax[1].plot(t[:,0], t[:,3])

            print('x error:', numpy.linalg.norm(t[:, 2] - rx))
            print('y error:', numpy.linalg.norm(t[:, 3] - ry))
        
        if not (args.compare_specfem is None):
            tx = numpy.loadtxt('%s.BXX.semd' % args.compare_specfem)
            ty = numpy.loadtxt('%s.BXZ.semd' % args.compare_specfem)

            # Use the peak to align seismograms
            i1 = numpy.argmax(tx[:, 1])
            i2 = numpy.argmax(rx)
            delta = seismogram[i2, 0] - tx[i1, 0]
        
            ax[0].plot(tx[:, 0] + delta, tx[:, 1])
            ax[1].plot(ty[:, 0] + delta, ty[:, 1])

            
        P.show()

    seismogram[:,2] = rx
    seismogram[:,3] = ry
    seismogram[:,4] = numpy.real(gx)
    seismogram[:,5] = numpy.real(gy)

    numpy.save(args.output, seismogram)

    
    
    
