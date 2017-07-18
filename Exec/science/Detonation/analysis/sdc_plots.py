from operator import itemgetter

import numpy as np
import matplotlib.pyplot as plt

COLORS = ["C0", "C1"]

class ReactStep(object):
    # a simple container to hold the result of integration

    def __init__(self, zone, t, rho, X4, X56, int_type="strang"):
        self.zone = zone
        self.t = np.array(t)
        self.rho = np.array(rho)
        self.X4 = np.array(X4)
        self.X56 = np.array(X56)
        self.int_type = int_type

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
        for line in block:
            data = process_line(line)

            if zone is None:
                zone = data[1]
            elif data[1] != zone:
                sys.exit("zone info changed")

            t.append(data[0])
            rho.append(data[2])
            X4.append(data[3])
            X56.append(data[4])

        # store the last iteration
        iters.append(ReactStep(zone, t, rho, X4, X56, int_type="Strang"))

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
                iters.append(ReactStep(zone, t, rho, X4, X56, int_type="SDC"))
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
        iters.append(ReactStep(zone, t, rho, X4, X56, int_type="SDC"))

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


def plot_sdc(data):

    it0 = [q[0] for q in data]
    it1 = [q[1] for q in data]

    max_change = [q.max_change() for q in it1]
    t_max = max(max_change, key=itemgetter(1))[0]

    dt = (np.array([q.t.max() for q in it1]).max() - np.array([q.t.min() for q in it1]).min())/len(it1)

    print(t_max, dt)

    c = 0
    for i0, i1 in zip(it0, it1):
        plt.plot(i0.t, i0.X4, ls=":", color=COLORS[c])
        plt.plot(i1.t, i1.X4, ls="-", color=COLORS[c])

        c ^= 1

    plt.xlim(t_max - 4*dt, t_max + 4*dt)

    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

    plt.savefig("sdc.png")


def plot_strang(data):

    sdata = [q[0] for q in data]

    max_change = [q.max_change() for q in sdata]
    t_max = max(max_change, key=itemgetter(1))[0]

    dt = (np.array([q.t.max() for q in sdata]).max() - np.array([q.t.min() for q in sdata]).min())/len(sdata)

    print(t_max, dt)

    c = 0
    for step in sdata:
        plt.plot(step.t, step.X4, ls=":", color=COLORS[c])
        plt.plot(step.t, step.X4, ls="-", color=COLORS[c])

        c ^= 1

    plt.xlim(t_max - 4*dt, t_max + 4*dt)

    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

    plt.savefig("strang.png")


if __name__ == "__main__":
    #data = read_file("SDC/zone_info.sdc.00400")
    data = read_file("Strang/zone_info.00400")

    if len(data[0]) > 1:
        plot_sdc(data)
    else:
        plot_strang(data)

