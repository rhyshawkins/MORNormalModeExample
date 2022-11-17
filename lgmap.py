
import numpy

class LocalGlobalMap2D:

    EPSILON = 1.0e-6
    
    def __init__(self):

        self.global_x = []
        self.global_y = []


    def register_node(self, x, y):

        index = self.__find(x, y)
        if (index < 0):
            index = self.__add(x, y)

        return index

    def nglob(self):
        return len(self.global_x)

    def points(self):

        pts = numpy.zeros((len(self.global_x), 2))
        pts[:,0] = numpy.array(self.global_x)
        pts[:,1] = numpy.array(self.global_y)

        return pts

    def __find(self, x, y):

        for i in range(len(self.global_x)):
            d = numpy.sqrt((self.global_x[i] - x)**2 + (self.global_y[i] - y)**2)
            if (d < LocalGlobalMap2D.EPSILON):
                return i

        return -1
            

    def __add(self, x, y):

        idx = len(self.global_x)
        self.global_x.append(x)
        self.global_y.append(y)
        return idx
            
