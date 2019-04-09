import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import argparse


def plot(x, y, z):
    mpl.rcParams['legend.fontsize'] = 10
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot(x, y, z, label='Track')
    ax.legend()
    plt.show()


def read(fname):
    xlist = []
    ylist = []
    zlist = []
    with open(fname, 'r') as file:
        for line in file:
            x, y, z = line.split()
            xlist.append(float(x))
            ylist.append(float(y))
            zlist.append(float(z))
    
    return xlist, ylist, zlist


def main():
    parser = argparse.ArgumentParser(description='plot particle track')
    parser.add_argument('--file_name', type=str, help='name of input file')

    file_name = parser.parse_args().file_name
    x, y, z = read(file_name)
    plot(x, y, z)

if __name__ == '__main__':
    main()


