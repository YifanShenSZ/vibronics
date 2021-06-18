import argparse
from pathlib import Path
import numpy as numpy
import matplotlib.pyplot as plt

def parse_args() -> argparse.Namespace: # Command line input
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('xy', type=Path, help='xy file')
    parser.add_argument( 'z', type=Path, help= 'z file')
    args = parser.parse_args()
    return args

def read_vector(file:Path) -> numpy.ndarray:
    with open(file, 'r') as f: lines = f.readlines()
    v = numpy.empty(len(lines))
    for i in range(len(lines)): v[i] = float(lines[i])
    return v

def read_matrix(file:Path) -> numpy.ndarray:
    with open(file, 'r') as f: lines = f.readlines()
    temp = lines[0].split()
    size0 = len(lines); size1 = len(temp)
    m = numpy.empty((size0, size1))
    for i in range(size0):
        temp = lines[i].split()
        for j in range(size1): m[i, j] = float(temp[j])
    return m

if __name__ == "__main__":
    args = parse_args()
    z_raw = read_vector(args.z)
    xy = read_matrix(args.xy)

    assert xy.shape[1] == 2
    x, xlist = numpy.unique(xy[:,0], return_inverse=True)
    y, ylist = numpy.unique(xy[:,1], return_inverse=True)
    z = numpy.empty((y.shape[0], x.shape[0]))
    for i in range(z_raw.shape[0]):
        z[ylist[i], xlist[i]] = z_raw[i]

    z[:11, 10:] *= -1

    with open("a-smooth-1.txt", 'w') as f:
        for j in range(z.shape[1]):
            for i in range(z.shape[0]):
                print("%25.15E"%z[i, j], file=f)

    fig, ax = plt.subplots()
    im = ax.imshow(z, cmap="inferno")
    ax.set_xticks(numpy.arange(x.shape[0])); ax.set_xticklabels(x)
    ax.set_yticks(numpy.arange(y.shape[0])); ax.set_yticklabels(y)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
    plt.colorbar(im)
    fig.tight_layout()

    plt.show()
