#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
from scipy.spatial import Voronoi, voronoi_plot_2d
from scipy.stats import entropy
import math
import os

class XY:
    def __int__(self):
        self.original = False
        self.xy = (0,0, 0.0)


class Point:
    def __init__(self, xy):
        self.xy = xy                    # coordinates of the point
        self.enclosed = True            # is this point enclosed?
        self.vertices = []              # enclosing vertices
        self.ridges = []                # enclosing ridges
        self.neighbor_points = []       # neighbor points
        self.area = 0                   # area of region if enclosed is true
        self.ridge_stdev = 0            # ridge length stdev
        self.angle_stdev = 0            # stdev of the angles of polygon composed by vetices


class Vertex:
    def __init__(self, xy):
        self.xy = xy                    # coordinates of the vertex
        self.ridges = []                # what ridges it is an end
        self.points = []                # which points it belongs to


class Ridge:
    def __init__(self, v1, v2):
        self.open = False               # with an open end?
        self.vertices = [v1, v2]        # two end vertices
        self.points = []                # two points this edge belongs to
        self.length = 0                 # length of the edge


def dvv(v1, v2):
    return np.linalg.norm(np.array(v1) - np.array(v2))


def read_coordinates(fname):
    # read point coordinates
    coordinates = []
    with open(fname) as fp:
        while True:
            line = fp.readline()
            if not line:
                break
            else:
                p = XY()
                p.original = True
                p.xy = tuple([float(x) for x in line.strip().split()])
                coordinates.append(p)
    return coordinates


def pad_unitcell(coordinates, strip=1.0, pad=0.0, shear=0.0):
    """
    Assuming the input points are in a box, this box is cloned and placed around the original.
    The strip size is how much or the original box is kept. 1 means all, 0.1 means 10%.
    The pad value is how much space is inserted in the boundary.
    The shear is the horizontal shift relative to the center made clock wise.
    The lower strip is shifted twice the amount and put on top.
    The upper strip is shifted the same way in opposite direction.
    After the shear shift, left and right are extended.
    """

def draw_points(coordinates):

    return

def read_coordinates_par(fname):
    # read point coordinates
    snapshots = []
    coordinates = []
    lines = open(fname).readlines()
    lines = lines[22:]
    for line in lines:
        if line[0] == "#":
            if coordinates:
                snapshots.append(coordinates)
            coordinates = []
        else:
            fields = line.strip().split()
            if len(fields) > 10:
                p = [float(fields[2]), float(fields[4])]
                coordinates.append(p)
    return snapshots

def load_points(coordinates):
    points = []
    vertices = []
    ridges = []

    for xy in coordinates:
        p = Point(xy)
        points.append(p)

    # get vertices, edges, and regions
    v_points = [p.xy for p in points]
    vor = Voronoi(v_points)

    for v in vor.vertices:
        vertex = Vertex(v)
        vertices.append(vertex)

    for v in vor.ridge_vertices:
        if v[0] == -1 or v[1] == -1:
            ridge = Ridge(None, None)
            ridge.open = True
        else:
            vertex0 = vertices[v[0]]
            vertex1 = vertices[v[1]]
            ridge = Ridge(vertex0, vertex1)
            ridge.length = dvv(vertex0.xy, vertex1.xy)

        vertices[v[0]].ridges.append(ridge)
        vertices[v[1]].ridges.append(ridge)
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
            for iv in regions[i_region]:
                vertices[iv].points.append(points[ip])

    for ir in range(len(vor.ridge_points)):
        ridges[ir].points = [points[ip] for ip in vor.ridge_points[ir]]
        ridges[ir].points[0].neighbor_points.append(ridges[ir].points[1])   # ridge defines two mutual neighbors
        ridges[ir].points[1].neighbor_points.append(ridges[ir].points[0])   # ridge defines two mutual neighbors
        for ip in vor.ridge_points[ir]:
            points[ip].ridges.append(ridges[ir])

    for p in points:
        if p.enclosed:
            # find area of each enclosed region
            area = 0.0
            for r in p.ridges:
                triangle_v = [p.xy, r.vertices[0].xy, r.vertices[1].xy]
                area += get_area_3v(triangle_v)
            p.area = area

            # find ridge lengths stdev
            ridge_lengths = [r.length for r in p.ridges]
            p.ridge_stdev = np.std(ridge_lengths)

            # find stdev of angles
            angles = []
            for v in p.vertices:
                sides = [s for s in v.ridges if s in p.ridges]
                v1 = np.array(sides[0].vertices[1].xy) - np.array(sides[0].vertices[0].xy)
                v2 = np.array(sides[1].vertices[1].xy) - np.array(sides[1].vertices[0].xy)
                angles.append(avv(v1, v2))
            p.angle_stdev = np.std(angles)

    return points, vertices, ridges


def get_area_3v(triangle_v):
    v0 = triangle_v[0]
    v1 = triangle_v[1]
    v2 = triangle_v[2]

    d1 = dvv(v0, v1)
    d2 = dvv(v1, v2)
    d3 = dvv(v0, v2)

    s = 0.5*(d1 +d2 + d3)
    A = math.sqrt(s*(s-d1)*(s-d2)*(s-d3))

    return A

def avv(v1, v2):
    "Angle between v1 to v2 in [0, pi]."
    v1_n = np.linalg.norm(v1)
    v2_n = np.linalg.norm(v2)
    COS = np.dot(v1, v2) / v1_n / v2_n

    return np.arccos(COS)

def plot_voronoi(points, center=[], enlarge=0):
    ps = np.array([p.xy for p in points])
    vor = Voronoi(ps)
    fig = voronoi_plot_2d(vor)
    if center:
        xmin = xmax = center[0]
        ymin = ymax = center[1]
    else:
        xmin = min(ps[:, 0])
        xmax = max(ps[:, 0])
        ymin = min(ps[:, 1])
        ymax = max(ps[:, 1])

    if enlarge:
        plt.xlim(xmin-enlarge, xmax+enlarge)
        plt.ylim(ymin-enlarge, ymax+enlarge)

    plt.show()
    return


def plot_voronoi_color(points, color_by="area", log=False, cmap="", enlarge=0.0, color_cut=1.0, boundary=True):
    ps = np.array([p.xy for p in points])
    vor = Voronoi(ps)
    if color_by.upper() == "AREA":
        c_scale = [p.area for p in points]
    elif color_by.upper() == "RIDGE_STDEV":
        c_scale = [p.ridge_stdev for p in points]
    elif color_by.upper() == "ANGLE_STDEV":
        c_scale = [p.angle_stdev for p in points]

    if log == True:
        c_scale = [math.log(c+1) for c in c_scale]

    if not cmap:
        cmap = cm.plasma

    if boundary:
        line_width = 1.0
    else:
        line_width = 0.0

    valid_c = [c for c in c_scale if abs(c) > 0.00001]
    if valid_c:
        color_min = min(valid_c)
        color_max = max(valid_c)
    else:
        color_min = 0.0
        color_max = 0.0    

    # normalize chosen colormap
    # https://matplotlib.org/3.5.0/tutorials/colors/colormaps.html
    norm = mpl.colors.Normalize(vmin=color_min, vmax=color_max, clip=True)
    mapper = cm.ScalarMappable(norm=norm, cmap=cmap)

    voronoi_plot_2d(vor, show_points=False, show_vertices=False, s=1, line_width=line_width)

    xmin = min(ps[:, 0])
    xmax = max(ps[:, 0])
    ymin = min(ps[:, 1])
    ymax = max(ps[:, 1])


    plt.xlim(xmin - enlarge, xmax + enlarge)
    plt.ylim(ymin - enlarge, ymax + enlarge)

    for r in range(len(vor.point_region)):
        region = vor.regions[vor.point_region[r]]
        if not -1 in region:
            polygon = [vor.vertices[i] for i in region]
            plt.fill(*zip(*polygon), color=mapper.to_rgba(c_scale[r]))

    plt.show()
    

def save_voronoi_color(points, color_by="area", log=False, cmap="", enlarge=0.0, color_cut=1.0, boundary=True):
    plot_voronoi_color.counter += 1
    ps = np.array([p.xy for p in points])
    vor = Voronoi(ps)
    if color_by.upper() == "AREA":
        c_scale = [p.area for p in points]
    elif color_by.upper() == "RIDGE_STDEV":
        c_scale = [p.ridge_stdev for p in points]
    elif color_by.upper() == "ANGLE_STDEV":
        c_scale = [p.angle_stdev for p in points]

    if log == True:
        c_scale = [math.log(c+1) for c in c_scale]

    if not cmap:
        cmap = cm.plasma

    if boundary:
        line_width = 1.0
    else:
        line_width = 0.0


    color_min = min([c for c in c_scale if abs(c) > 0.00001])
    color_max = (max([c for c in c_scale if abs(c) > 0.00001]) - color_min) * color_cut + color_min

    # normalize chosen colormap
    # https://matplotlib.org/3.5.0/tutorials/colors/colormaps.html
    norm = mpl.colors.Normalize(vmin=color_min, vmax=color_max, clip=True)
    mapper = cm.ScalarMappable(norm=norm, cmap=cmap)

    voronoi_plot_2d(vor, show_points=False, show_vertices=False, s=1, line_width=line_width)

    xmin = min(ps[:, 0])
    xmax = max(ps[:, 0])
    ymin = min(ps[:, 1])
    ymax = max(ps[:, 1])

    xmin = -25
    xmax = 25
    ymin = -25
    ymax = 25


    plt.xlim(xmin - enlarge, xmax + enlarge)
    plt.ylim(ymin - enlarge, ymax + enlarge)

    for r in range(len(vor.point_region)):
        region = vor.regions[vor.point_region[r]]
        if not -1 in region:
            polygon = [vor.vertices[i] for i in region]
            plt.fill(*zip(*polygon), color=mapper.to_rgba(c_scale[r]))

    filename = "%04d.png" % plot_voronoi_color.counter
    dirname = "images_angle"
    if not os.path.isdir(dirname):
        os.mkdir(dirname)
    fpath = "/".join([dirname,filename])
    print("Writing %s" % fpath)
    plt.savefig(fpath)
    #plt.show()


if __name__ == '__main__':
#    inputfile = "points_9.txt"
    inputfile = "random100.txt"
    coordinates = read_coordinates(inputfile)
    coordinates_ext = pad_unitcell(coordinates, pad=0.0, shear=0.1)

    draw_points(coordinates_ext)



    #points, vertices, ridges = load_points_extended(coordinates_ext)
