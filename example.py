
import numpy

import regularmesh

if __name__ == '__main__':

    #
    # Construct an homogeneous mesh of 1 km^2 as an example for regional normal modes
    #

    rho = 2.6e3
    vp = 5.8e3
    vs = 3.1986e3

    lgmap, cells = regularmesh.build_regular_mesh(5, 5, 1.0e3, 1.0e3, rho, vs, vp)

    M, K, C = regularmesh.build_matrices(lgmap, cells)

    numpy.save('paperM.npy', M)
    numpy.save('paperK.npy', K)
    numpy.save('paperC.npy', C)

    numpy.save('paperPoints.npy', lgmap.points())
    
    
