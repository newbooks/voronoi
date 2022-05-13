#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi, voronoi_plot_2d
from scipy.stats import entropy

class Point:
    def __init__(self):
        self.xy = ()                   # coordinates of the point
        self.enclosed = True            # is this point enclosed?
        self.neighbor_points = []       # neighbor points
        self.vertices = []              # enclosing vertices
        self.edges = []                 # enclosing edges
        self.angles = []                # angles formed by edges in radian
        self.angle_stdev = 0            # angle standard deviation

class Vertex:
    def __init__(self):
        self.xy = ()                   # coordinates of the vertex
        self.edges = []                 # what edges it is an end
        self.points = []                # which points it belongs to

class Edge:
    def __init__(self):
        self.ends = ()                  # two end points
        self.length = 0                 # length of the edge

def read_points(fname):
    points = []
    with open(fname) as fp:
        while True:
            line = fp.readline()
            if not line:
                break
            else:
                p = Point()
                p.xy = tuple([float(x) for x in line.strip().split()])
                points.append(p)
    return points

if __name__ == '__main__':
    inputfile = "random100.txt"
    points = read_points(inputfile)
    for p in points:
        print(p.xy)