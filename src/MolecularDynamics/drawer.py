import os
import sys
import argparse
try:
    import numpy as np
    import matplotlib.pyplot as plt
except ImportError as e:
    print('Failed in import site packages')
    print(e)
    exit(1)

# creater argument parser
parser = argparse.ArgumentParser(description='Drawing script')
parser.add_argument('-f', metavar='--file', dest='file', type=str,
                    required=True, help='file path')

# parse args
args = parser.parse_args(sys.argv[1:])

if not os.path.exists(args.file):
    print(f'Can\'t find file {args.file}')
    exit(2)

# read data
with open(args.file, 'r') as reader:
    try:
        n = int(reader.readline())
    except Exception as e:
        print(e)
        exit(3)

    coords = np.zeros((n, 3), dtype=float)

    for i in range(n):
        line = reader.readline()
        try:
            coords[i] = np.array(line.split(), dtype=float)
        except Exception as e:
            print(e)
            exit(3)

ax = plt.figure(figsize=(8,6)).add_subplot(projection='3d')

ax.plot(coords[:,0], coords[:,1], coords[:,2])

plt.show()