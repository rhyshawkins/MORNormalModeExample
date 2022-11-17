
import numpy
import cardinal

class Quad:

    ABSORB_LEFT = 1
    ABSORB_BOTTOM = 2
    ABSORB_RIGHT = 4
    
    def __init__(self, order, quadrature, rho, mu, lmda, localglobalmap, x1, x2, y1, y2, absorb = 0):

        #
        # Construct nodes and
        #
        self.rho = rho
        self.mu = mu
        self.lmda = lmda
        self.order = order
        self.absorb = absorb
        
        self.global_indices = numpy.zeros((order + 1, order + 1), dtype = 'int')
        self.node_x = numpy.zeros((order + 1, order + 1))
        self.node_y = numpy.zeros((order + 1, order + 1))
        
        for j in range(order + 1):

            y = y1 + (y2 - y1) * (quadrature.lobatto_node[j] + 1.0)/2.0
            
            for i in range(order + 1):

                x = x1 + (x2 - x1) * (quadrature.lobatto_node[i] + 1.0)/2.0
                
                self.node_x[j, i] = x
                self.node_y[j, i] = y

                self.global_indices[j, i] = localglobalmap.register_node(x, y)

        self.J = (x2 - x1)*(y2 - y1)/4.0
        self.dxidx = 2.0/(x2 - x1)
        self.detady = 2.0/(y2 - y1)


    def compute_matrices(self, quadrature, M, K, C):

        nglob = M.size//2

        for j in range(self.order + 1):

            weight_eta = quadrature.lobatto_weight[j]
            
            for i in range(self.order + 1):

                weight_xi = quadrature.lobatto_weight[i];


                #
                # Mass matrix encoded as a vector
                #
                v = weight_xi * weight_eta * self.J * self.rho

                gi = self.global_indices[j, i]
                M[gi] += v
                M[nglob + gi] += v

                #
                # Stiffness matrix
                #
                scale = weight_xi * weight_eta * self.J
                self.__phi_xx(nglob, quadrature, K, i, j,
                              [scale * (self.lmda + 2.0*self.mu), scale * self.lmda],
                              [self.__uxx, self.__uyy])

                self.__phi_xy(nglob, quadrature, K, i, j,
                              [scale * self.mu, scale * self.mu],
                              [self.__uxy, self.__uyx])

                self.__phi_yx(nglob, quadrature, K, i, j,
                              [scale * self.mu, scale * self.mu],
                              [self.__uxy, self.__uyx])

                self.__phi_yy(nglob, quadrature, K, i, j,
                              [scale * self.lmda, scale * (self.lmda + 2.0*self.mu)],
                              [self.__uxx, self.__uyy])

        #
        # Absorbing matrix
        #
        self.__stacey(nglob, quadrature, C)
                


    def __phi_xx(self, nglob, quadrature, K, i, j, scaleij, strainij):

        for scale, strain in zip(scaleij, strainij):
            for m in range(self.order + 1):
                row = self.global_indices[j, m]
                strain(nglob, quadrature, K, i, j, row, self.dxidx * quadrature.lobatto_derivative_weight[m, i] * scale)

    def __phi_xy(self, nglob, quadrature, K, i, j, scaleij, strainij):

        for scale, strain in zip(scaleij, strainij):
            for m in range(self.order + 1):
                row = self.global_indices[m, i]
                strain(nglob, quadrature, K, i, j, row, self.detady * quadrature.lobatto_derivative_weight[m, j] * scale)

    def __phi_yx(self, nglob, quadrature, K, i, j, scaleij, strainij):

        for scale, strain in zip(scaleij, strainij):
            for m in range(self.order + 1):
                row = self.global_indices[j, m] + nglob
                strain(nglob, quadrature, K, i, j, row, self.dxidx * quadrature.lobatto_derivative_weight[m, i] * scale)

    def __phi_yy(self, nglob, quadrature, K, i, j, scaleij, strainij):

        for scale, strain in zip(scaleij, strainij):
            for m in range(self.order + 1):
                row = self.global_indices[m, i] + nglob
                strain(nglob, quadrature, K, i, j, row, self.detady * quadrature.lobatto_derivative_weight[m, j] * scale)

    def __uxx(self, nglob, quadrature, K, i, j, row, scale):

        for n in range(self.order + 1):
            col = self.global_indices[j, n]
            v = scale * self.dxidx * quadrature.lobatto_derivative_weight[n, i]

            K[row, col] += v

    def __uyx(self, nglob, quadrature, K, i, j, row, scale):

        for n in range(self.order + 1):
            col = self.global_indices[j, n] + nglob
            v = scale * self.dxidx * quadrature.lobatto_derivative_weight[n, i]

            K[row, col] += v
            
    def __uyy(self, nglob, quadrature, K, i, j, row, scale):
        
        for n in range(self.order + 1):
            col = self.global_indices[n, i] + nglob
            v = scale * self.detady * quadrature.lobatto_derivative_weight[n, j]

            K[row, col] += v
            
    def __uxy(self, nglob, quadrature, K, i, j, row, scale):
        
        for n in range(self.order + 1):
            col = self.global_indices[n, i]
            v = scale * self.detady * quadrature.lobatto_derivative_weight[n, j]

            K[row, col] += v
    
    def __stacey(self, nglob, quadrature, C):

        if (self.absorb & Quad.ABSORB_LEFT) > 0:

            i = 0
            for j in range(self.order + 1):

                xg = 0.0
                yg = self.dxidx * self.J
                J1D = numpy.sqrt(xg*xg + yg*yg)
                nx = -yg/J1D
                ny = xg/J1D

                txx, tyy = self.__stacey_traction(nx, ny)

                xi = self.global_indices[j, i]
                yi = xi + nglob

                C[xi, xi] += J1D * quadrature.lobatto_weight[j] * txx

                C[yi, yi] += J1D * quadrature.lobatto_weight[j] * tyy

        if (self.absorb & Quad.ABSORB_RIGHT) > 0:
            
            i = self.order
            for j in range(self.order + 1):

                xg = 0.0
                yg = self.dxidx * self.J
                J1D = numpy.sqrt(xg*xg + yg*yg)
                nx = yg/J1D
                ny = -xg/J1D

                txx, tyy = self.__stacey_traction(nx, ny)

                xi = self.global_indices[j, i]
                yi = xi + nglob

                C[xi, xi] += J1D * quadrature.lobatto_weight[j] * txx

                C[yi, yi] += J1D * quadrature.lobatto_weight[j] * tyy

        if (self.absorb & Quad.ABSORB_BOTTOM) > 0:

            j = 0

            istart = 0
            if (self.absorb & Quad.ABSORB_LEFT) > 0:
                istart = 1

            iend = self.order + 1
            if (self.absorb & Quad.ABSORB_RIGHT) > 0:
                iend = self.order
            
            for i in range(istart, iend):
                xg = self.detady * self.J
                yg = 0.0
                J1D = numpy.sqrt(xg*xg + yg*yg)
                nx = yg/J1D
                ny = -xg/J1D

                txx, tyy = self.__stacey_traction(nx, ny)

                
                xi = self.global_indices[j, i]
                yi = xi + nglob

                C[xi, xi] += J1D * quadrature.lobatto_weight[i] * txx

                C[yi, yi] += J1D * quadrature.lobatto_weight[i] * tyy

    def __stacey_traction(self, nx, ny):

        lp2mu = self.lmda + 2.0*self.mu
        vp = numpy.sqrt(lp2mu/self.rho)
        vs = numpy.sqrt(self.mu/self.rho)

        txx = self.rho*((vp - vs)*nx*nx + vs)
        tyy = self.rho*((vp - vs)*ny*ny + vs)

        return txx, tyy
