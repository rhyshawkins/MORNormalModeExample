
import argparse

import numpy
import matplotlib.pyplot as P

import scipy.linalg

import hashlib

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input', type = str, required = True, help = 'Input file prefix')

    parser.add_argument('-o', '--output', type = str, required = True, help = 'Output file prefix')

    parser.add_argument('-n', '--no-scale', action = 'store_true', default = False, help = 'Disable FLV scaling')
    
    parser.add_argument('-s', '--show', action = 'store_true', default = False, help = 'Plot eigenvalues')
    

    args = parser.parse_args()

    M = numpy.load('%sM.npy' % args.input)
    K = numpy.load('%sK.npy' % args.input)
    C = numpy.load('%sC.npy' % args.input)

    size = M.size

    nM = numpy.linalg.norm(M)
    nK = numpy.linalg.norm(K)
    nC = numpy.linalg.norm(C)

    if args.no_scale:
        gamma = 1.0
        delta = 1.0
    else:
        gamma = numpy.sqrt(nK/nM)
        delta = 2.0/(nM + nC*gamma)

    M *= gamma*gamma*delta
    C *= gamma*delta
    K *= delta

    A = numpy.zeros((2*size, 2*size))

    A[:size, :size] = C
    A[:size, size:] = K
    A[size:, :size] = -numpy.eye(size)

    B = numpy.zeros((2*size, 2*size))
    B[:size, :size] = -numpy.diag(M)
    B[size:, size:] = -numpy.eye(size)

    eigvalues, left, right = scipy.linalg.eig(A, B, left = True, right = True, overwrite_a = False, overwrite_b = False)

    #
    # Normalize eigenvectors so that L^* M R = I
    #
    for i in range(eigvalues.size):

        n = left[:, i].conj().dot(B.dot(right[:, i]))
        if numpy.imag(eigvalues[i]) != 0.0:
            n = numpy.sqrt(n)
            left[:, i] /= numpy.conj(n)
            right[:, i] /= n

        else:
            if n < 0.0:
                n = numpy.sqrt(-numpy.real(n))
                left[:, i] /= -n
                right[:, i] /= n
            else:
                n = numpy.sqrt(numpy.real(n))
                left[:, i] /= n
                right[:, i] /= n
            
    T = left.conj().T.dot(B.dot(right))

    numpy.save('%s-eigval.npy' % args.output, eigvalues)
    numpy.save('%s-eigleft.npy' % args.output, left)
    numpy.save('%s-eigright.npy' % args.output, right)
    scale = numpy.array([gamma, delta])
    numpy.save('%s-eigscale.npy' % args.output, scale)
    
    if args.show:
        fig, ax = P.subplots()
        
        ax.scatter(gamma * numpy.real(eigvalues), gamma * numpy.imag(eigvalues))
        
        P.show()

    
    
