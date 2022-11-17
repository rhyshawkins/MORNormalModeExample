
import numpy

import lgmap
import quad

import cardinal

def build_regular_mesh(N, M, width, height, rho, vs, vp, order = 4):

    mapper = lgmap.LocalGlobalMap2D()
    cells = []

    cw = width/N
    ch = height/M

    mu = rho * vs**2
    lmda = rho * vp**2 - 2.0*mu

    quadrature = cardinal.Quadrature(order)
    
    for j in range(M):

        y1 = j * ch
        y2 = y1 + ch

        if j == 0:
            absorb_bottom = quad.Quad.ABSORB_BOTTOM
        else:
            absorb_bottom = 0
        
        for i in range(N):

            x1 = i * cw
            x2 = x1 + cw

            absorb = absorb_bottom
            if i == 0:
                absorb |= quad.Quad.ABSORB_LEFT
            elif i == (N - 1):
                absorb |= quad.Quad.ABSORB_RIGHT

            cells.append(quad.Quad(order, quadrature, rho, mu, lmda, mapper, x1, x2, y1, y2, absorb))

    return mapper, cells

def build_matrices(mapper, cells, order = 4):

    nglob = mapper.nglob()
    ndof = 2*nglob

    quadrature = cardinal.Quadrature(order)
    
    M = numpy.zeros((ndof,))
    K = numpy.zeros((ndof, ndof))
    C = numpy.zeros((ndof, ndof))
    
    for cell in cells:

        cell.compute_matrices(quadrature, M, K, C)

    return M, K, C
    

                         

            
