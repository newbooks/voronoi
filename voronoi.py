#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi, voronoi_plot_2d
from scipy.stats import entropy

class Point:
    def __init__(self, xy):
        self.xy = xy                    # coordinates of the point
        self.enclosed = True            # is this point enclosed?
        self.vertices = []              # enclosing vertices
        self.ridges = []                # enclosing ridges
        self.neighbor_points = []       # neighbor points

class Vertex:
    def __init__(self, xy):
        self.xy = xy                    # coordinates of the vertex
        self.ridges = []                # what ridges it is an end
        self.points = []                # which points it belongs to

class Ridge:
    def __init__(self, v1, v2):
        self.vertices = [v1, v2]        # two end vertices
        self.points = []                # two points this edge belongs to
        self.length = 0                 # length of the edge

def read_points(fname):
    points = []
    vertices = []
    ridges = []
    # read point coordinates
    with open(fname) as fp:
        while True:
            line = fp.readline()
            if not line:
                break
            else:
                p = Point([float(x) for x in line.strip().split()])
                points.append(p)

    # get vertices, edges, and regions
    v_points = [p.xy for p in points]
    vor = Voronoi(v_points)

    for v in vor.vertices:
        vertex = Vertex(v)
        vertices.append(vertex)

    for v in vor.ridge_vertices:
        if v[0] == -1:
            continue
        else:
            xy0 = vertices[v[0]].xy
        if v[1] == -1:
            continue
        else:
            xy1 = vertices[v[1]].xy
        ridge = Ridge(xy0, xy1)
        ridges.append(ridge)

    # make relationship
    regions = vor.regions
    for ip in range(len(points)):
        i_region = vor.point_region[ip]
        if -1 in regions[i_region]:     # an open point
            points[ip].enclosed = False
        else:
            points[ip].enclosed = True
            points[ip].vertices = [vertices[i] for i in regions[i_region]]

    

    return points, vertices, ridges




if __name__ == '__main__':
    inputfile = "points_9.txt"
    points, vertices, ridges = read_points(inputfile)
    print("Points")
    for p in points:
        print("   Corrdinates:", p.xy)
        print("   Vertices:")
        if p.vertices:
            print([pv.xy for pv in p.vertices])
    print("Vertices")
    for v in vertices:
        print(v.xy)

    print("Ridges")
    for v in ridges:
        print(v.vertices[0], " - ", v.vertices[1])

