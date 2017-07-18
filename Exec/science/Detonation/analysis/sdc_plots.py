from operator import itemgetter

import numpy as np
import matplotlib.pyplot as plt

class ReactStep(object):
    # a simple container to hold the result of integration

    def __init__(self, zone, t, rho, X4, X56):
        self.zone = zone
        self.t = np.array(t)
        self.rho = np.array(rho)
        self.X4 = np.array(X4)
        self.X56 = np.array(X56)

    def max_change(self):
        return self.t.min(), abs(self.X4.max() - abs(self.X4.min()))


def process_line(line):
    return [float(q) for q in line.split()]

def process_block(block):
    # take a block that represents an entire step of Strang
    # integration (both dt/2 halves) or all iterations of an SDC
    # integration and create ReactStep objects.  For Strang, only a
    # single object is created.  For SDC, one object is created for
    # each SDC iteration.

    if len(block) == 0: 
        return None

    zone = None
    sdc_iter = None
    t = []
    rho = []
    X4 = []
    X56 = []

    iters = []

    # are we SDC or Strang?
    nfields = len(block[0].split())

    if nfields == 5:
        # Strang
        pass

    elif nfields == 6:
        # for SDC, the main catch is that we keep track of which
        # iteration it is, and ensure that a new object is created for
        # each iteration

        for line in block:
            data = process_line(line)

            if zone is None:
                zone = data[2]
            elif data[2] != zone:
                sys.exit("zone info changed")

            if sdc_iter is None:
                sdc_iter = data[1]
            elif data[1] != sdc_iter:
                # new iteration -- store the current
                iters.append(ReactStep(zone, t, rho, X4, X56))
                sdc_iter = data[1]

                # start the new iteration's data
                t = [data[0]]
                rho = [data[3]]
                X4 = [data[4]]
                X56 = [data[5]]

            t.append(data[0])
            rho.append(data[3])
            X4.append(data[4])
            X56.append(data[5])

        # store the last iteration
        iters.append(ReactStep(zone, t, rho, X4, X56))        

    return iters

def read_file(filename):

    with open(filename, "r") as f:

        # hold all of the steps defined in the file
        data = []
        
        # hold the group of lines (a block) that make up a single step
        block = []

        line = f.readline()
        while line:
            if line.strip().startswith("#"):
                # done with a block
                step = process_block(block)
                if step is not None:
                    data.append(step)

                block = []
                
            else:
                block.append(str(line.strip()))

            line = f.readline()

    return data


if __name__ == "__main__":
    data = read_file("SDC/zone_info.sdc.00400")

    it0 = [q[0] for q in data]
    it1 = [q[1] for q in data]

    max_change = [q.max_change() for q in it1]
    t_max = max(max_change, key=itemgetter(1))[0]

    dt = (np.array([q.t.max() for q in it1]).max() - np.array([q.t.min() for q in it1]).min())/len(it1)

    print(t_max, dt)

    c = 0
    for i0, i1 in zip(it0, it1):
        
        if c == 0:
            color = "C0"
        else:
            color = "C1"

        c ^= 1

        plt.plot(i0.t, i0.X4, ls=":", color=color, marker="x")
        plt.plot(i1.t, i1.X4, ls="-", color=color, marker="x")


    plt.xlim(t_max - 4*dt, t_max + 4*dt)

    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

    plt.savefig("sdc.png")


