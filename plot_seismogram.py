
import argparse

import numpy
import matplotlib.pyplot as P

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input', type = str, required = True, help = 'Input file')

    parser.add_argument('-c', '--compare', type = str, default = None, help = 'Compare with npy')

    args = parser.parse_args()

    seismogram = numpy.load(args.input)

    fig, ax = P.subplots(2, 1)

    ax[0].plot(seismogram[:, 0], seismogram[:, 2])
    ax[1].plot(seismogram[:, 0], seismogram[:, 3])

    if not (args.compare is None):

        compare = numpy.load(args.compare)
        ax[0].plot(compare[:, 0], compare[:, 2], 'r:')
        ax[1].plot(compare[:, 0], compare[:, 3], 'r:')

    P.show()
