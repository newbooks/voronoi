#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
from scipy.spatial import Voronoi, voronoi_plot_2d
from scipy.stats import entropy
import math
import os

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
                p = [float(x) for x in line.strip().split()]
                coordinates.append(p)
    return coordinates


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

def unitcell_expand(coordinates, pad=0.0, shear=0.0):
    expanded_points = []

    # Pad the box
    base_points = np.array(coordinates)
    xmin = min(base_points[:, 0])
    xmax = max(base_points[:, 0])
    ymin = min(base_points[:, 1])
    ymax = max(base_points[:, 1])
    unit_area = (ymax-ymin)*(xmax-xmin) / len(base_points)
    d_padding = np.sqrt(unit_area)
    if abs(pad) > 0.000001: # none 0 input padding, otherwise use average distance as padding
        d_padding = pad

    delta_x = xmax - xmin + d_padding
    delta_y = ymax - ymin + d_padding

    shift = np.array((-delta_x, -delta_y))
    extended_points = base_points + shift
    #print(base_points)
    #print(shift)
    print(extended_points)
    # shift = np.array((-delta_x, 0))
    # extended_points = np.append(extended_points, base_points + shift, axis=0)
    # shift = np.array((-delta_x, delta_y))
    # extended_points = np.append(extended_points, base_points + shift, axis=0)
    # shift = np.array((0, -delta_y))
    # extended_points = np.append(extended_points, base_points + shift, axis=0)
    # shift = np.array((0, delta_y))
    # extended_points = np.append(extended_points, base_points + shift, axis=0)
    # shift = np.array((delta_x, -delta_y))
    # extended_points = np.append(extended_points, base_points + shift, axis=0)
    # shift = np.array((delta_x, 0))
    # extended_points = np.append(extended_points, base_points + shift, axis=0)
    # shift = np.array((delta_x, delta_y))
    # extended_points = np.append(extended_points, base_points + shift, axis=0)
    #
    return expanded_points

if __name__ == '__main__':
#    inputfile = "points_9.txt"
    inputfile = "random100.txt"

    coordinates = read_coordinates(inputfile)
    expanded_coordinates = unitcell_expand(coordinates, pad=0.0, shear=0.0)
#    inputfile = "par_D2N500VF0.78Bidi1.4_0.5Square_18_nobrownian_2D_stress1.5r.dat"
#    snapshots = read_coordinates_par(inputfile)

# print("Verifying the data structure")
    # print("Points")
    # for p in points:
    #     print("   Point Corrdinates:", p.xy)
    #     if p.enclosed:
    #         print("      Enclosed: True")
    #     else:
    #         print("      Enclosed: False")
    #     print("   Vertices of point: ")
    #     if p.vertices:
    #         print("      ", [pv.xy for pv in p.vertices])
    #     else:
    #         print("      Open region doesn't have closing vertices.")
    #     print("   Ridges of point: ")
    #     for r in p.ridges:
    #         if not r.open:
    #             print("      ", [v.xy for v in r.vertices])
    #         else:
    #             print("      Open ridge")
    #     print("   Neighbors:")
    #     for n in p.neighbor_points:
    #         print("      ", n.xy)
    #     print("   Area:%.3f" % p.area)
    #     print()
    #
    # print()
    # print("Vertices")
    # for v in vertices:
    #     print("   Enclosed region Vertex Corrdinates:", v.xy)
    #     print("   Points of this vertex:")
    #     for p in v.points:
    #         print("      ", p.xy)
    #     print("   Ridges coordinates of this vertex:")
    #     for r in v.ridges:
    #         if not r.open:
    #             print("      ", [vofr.xy for vofr in r.vertices])
    #     print()
    #
    # print()
    # print("Ridges")
    # for r in ridges:
    #     if not r.open:
    #         print("   Ridge", r.vertices[0].xy, " - ", r.vertices[1].xy)
    #         print("   Length:", r.length)
    #         print("   Points it belongs:")
    #         for p in r.points:
    #             print("      ", p.xy)
    #         print()
    #
    #
    # # How chaotic is neighbor number distribution?
    # print("Neighbor(all) numbers:", neighbor_counts_all)
    # print()
    # print("Neighbor(enclosed region) numbers:", neighbor_counts_closed)
    #
    # print("%5s %8s %8s %8s" % ("#", "Area", "d_Ridge", "d_Angle"))
    # for i in range(len(points)):
    #     p = points[i]
    #     if p.enclosed:
    #         print("%5d %8.3f %8.3f %8.3f" % (i, p.area, p.ridge_stdev, p.angle_stdev))
    #     else:
    #         print("%5d Open region" % (i))

    # A list of cmap colors is available http
    # https://matplotlib.org/3.5.0/tutorials/colors/colormaps.html

    #neighbor_counts_all = [len(p.neighbor_points) for p in points]
    #neighbor_counts_closed = [len(p.neighbor_points) for p in points if p.enclosed]
    #print("Neighbor(all) number distribution entropy: %.3f" % entropy(neighbor_counts_all))
    #print("Neighbor(enclosed region) number distribution entropy: %.3f" % entropy(neighbor_counts_closed))
    #plot_voronoi_color(points, color_by="area", log=True, cmap="Blues_r", color_cut=0.1, enlarge=0.5)

    # plot_voronoi_color.counter = 0
    # for coordinates in snapshots:
    #     points, vertices, ridges = load_points(coordinates)
    #
    #     # neighbor_counts_all = [len(p.neighbor_points) for p in points]
    #     # neighbor_counts_closed = [len(p.neighbor_points) for p in points if p.enclosed]
    #     # print("Neighbor(all) number distribution entropy: %.3f" % entropy(neighbor_counts_all))
    #     # print("Neighbor(enclosed region) number distribution entropy: %.3f" % entropy(neighbor_counts_closed))
    #     save_voronoi_color(points, color_by="angle_stdev", log=True, cmap="Blues_r", color_cut=1, enlarge=0.5)
