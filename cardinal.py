
import numpy
from scipy.special import jacobi

def polysolve(x, y):

    order = len(x)

    V = numpy.zeros((order, order))
    for i in range(order):
        V[:,i] = x**(order - 1 - i)

    return numpy.linalg.solve(V, y)

def runge(x):
    return 1.0/(1.0 + 25.0*x**2)

def lobattonodes(order):
    t = numpy.ones((order + 1,))
    t[1:order] = numpy.sort(lobattopoly(order).r)
    t[0] = -1.0
    return t

def lobattopoly(order):
    return jacobi(order - 1, 1.0, 1.0) * (order + 1.0)/2.0

def lobattocardinal(order):

    xp = lobattonodes(order)
    
    ps = []
    for i in range(len(xp)):

        yp = numpy.zeros((len(xp),))
        yp[i] = 1.0

        p = numpy.poly1d(polysolve(xp, yp))
        ps.append(p)

    return xp, ps

def lobattoint(ps):
    pps = ps.integ()
    
    return pps(1.0) - pps(-1.0)

def lobattoweights(order):

    xp, ps = lobattocardinal(order)

    weights = list(map(lobattoint, ps))
    return weights

def lobattoderivativeweights(order):

    xp, ps = lobattocardinal(order)

    M = numpy.zeros((xp.size, xp.size))
    for j in range(xp.size):

        dd = ps[j].deriv()
        for i in range(xp.size):

            M[j, i] = dd(xp[i])

    return M
    
class Quadrature:

    def __init__(self, order):

        self.lobatto_node = lobattonodes(order)
        self.lobatto_weight = lobattoweights(order)
        self.lobatto_derivative_weight = lobattoderivativeweights(order)

if __name__ == '__main__':

    print(lobattonodes(4))
    print(lobattoweights(4))
