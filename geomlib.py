#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi, voronoi_plot_2d
from scipy.spatial import Delaunay
from scipy.stats import entropy
import math
import random

def draw_voronoi(points):
    vor = Voronoi(points)
    fig = voronoi_plot_2d(vor, show_points=True, show_vertices=True)
    fig.set_size_inches(8, 8)
    plt.axis("equal")

    plt.show()
    return

def draw_delaunay(points):
    tri = Delaunay(points)
    ps = np.array(points)
    fig = plt.figure()
    fig.set_size_inches(8, 8)
    plt.axis("equal")
    plt.triplot(ps[:,0], ps[:,1], tri.simplices)
    plt.plot(ps[:,0], ps[:,1], "o")
    plt.show()
    return

def draw_both(points):
    ps = np.array(points)

    f = plt.figure(figsize=(16, 6))
    ax1 = f.add_subplot(121)
    ax2 = f.add_subplot(122)

    vor = Voronoi(ps)
    voronoi_plot_2d(vor, show_points=True, show_vertices=True, s=4, ax=ax1)

    tri = Delaunay(ps)
    ax2.triplot(ps[:,0], ps[:,1], tri.simplices)
    ax2.plot(ps[:,0], ps[:,1], 'o')

    plt.show()

    return

def neighbors_vor(points):
    # return a dictionary that contains a histogram of neighbor number -> counts
    vor = Voronoi(points)
    neighbors = {}
    for region in vor.regions:
        if not region:
            continue
        n = len(region)
        if n in neighbors:
            neighbors[n] += 1
        else:
            neighbors[n] = 1

    return neighbors


def neighbors_tri(points):
    # return a dictionary that contains a histogram of neighbor number -> counts
    ps = np.array(points)
    tri = Delaunay(ps)
    point_nbr = {}
    for simplice in tri.simplices:
        if simplice[0] in point_nbr:
            point_nbr[simplice[0]].update({simplice[1], simplice[2]})
        else:
            point_nbr[simplice[0]] = {simplice[1], simplice[2]}
        if simplice[1] in point_nbr:
            point_nbr[simplice[1]].update({simplice[0], simplice[2]})
        else:
            point_nbr[simplice[1]] = {simplice[0], simplice[2]}
        if simplice[2] in point_nbr:
            point_nbr[simplice[2]].update({simplice[1], simplice[0]})
        else:
            point_nbr[simplice[2]] = {simplice[1], simplice[0]}

    neighbors = {}
    for p in point_nbr.keys():
        n = len(point_nbr[p])
        if n in neighbors:
            neighbors[n] += 1
        else:
            neighbors[n] = 1

    return neighbors


def grid_points(n=100):
    # put number of points on a grid close to a square shape, distance = 1
    # Use shake_points(step=100) to randomize it
    s = math.floor(math.sqrt(n))
    if n-s*s < 0.1:  # sqaure number was given
        nd = s
    else:
        nd = s + 1

    grid = []
    for i in range(nd):
        for j in range(nd):
            if i*nd + j >= n:
                return grid
            grid.append((float(i), float(j)))

    return grid


def no_conflict(p, ps, i_order, d):
    no_clash = True
    dd = d*d
    for i in range(len(ps)):
        if i == i_order: continue
        dx = ps[i][0] - p[0]
        dy = ps[i][1] - p[1]
        d2 = dx*dx + dy*dy
        if d2 < dd:
            no_clash = False
            break

    return no_clash


def shake_points(points, steps=100, d=0.1):
    # randomize points within original box, keep points distance over min
    # Find a box
    ps = np.array(points)
    xmin = min(ps[:,0])-1
    xmax = max(ps[:,0])+1
    ymin = min(ps[:,1])-1
    ymax = max(ps[:,1])+1

    for i_step in range(steps):
        # shuffle the move order
        order = [j for j in range(len(ps))]
        random.shuffle(order)
        for i_order in order:
            stay = True
            for i_trial in range(10):   # Try this number of times
                dx = np.random.normal(0.0, 0.1)
                dy = np.random.normal(0.0, 0.1)
                nx = ps[i_order][0] + dx
                ny = ps[i_order][1] + dy
                while nx < xmin or nx > xmax:
                    dx = np.random.normal(0.0, 0.1)
                    nx = ps[i_order][0] + dx
                while ny < ymin or ny > ymax:
                    dy = np.random.normal(0.0, 0.1)
                    ny = ps[i_order][1] + dy

                if no_conflict((nx, ny), ps, i_order, d):
                    stay = False
                    break
            if not stay:
                ps[i_order] = (nx, ny)
    return ps


def draw_histogram(nbrs):
    keys = list(nbrs.keys())
    keys.sort()
    values = []
    key_range = list(range(keys[0], keys[-1]+1))
    for key in key_range:
        if key in nbrs:
            values.append(nbrs[key])
        else:
            values.append(0)
    fig = plt.figure()
    fig.set_size_inches(8, 8)
    plt.bar(key_range, values)
    plt.xlabel("Neighbors")
    plt.ylabel("Counts")
    plt.show()
    return

if __name__ == '__main__':
    # square = [(0, 0), (0, 1), (1, 1), (1, 0)]
    # draw_voronoi(square)
    # draw_delaunay(square)
    # draw_both(square)
    # five = square + [(0.5, 0.5)]
    # his_vor = list(neighbors_vor(five).values())
    # his_tri = list(neighbors_tri(five).values())
    # print(entropy(his_vor), entropy(his_tri))
    points = grid_points(n=100)
    #nbrs = neighbors_vor(points)
    #his_vor = list(nbrs.values())
    #print(nbrs, entropy(his_vor))
    #draw_both(points)

    points = shake_points(points, steps=100, d=0.5)
    nbrs = neighbors_vor(points)
    his_vor = list(nbrs.values())
    print(nbrs, entropy(his_vor))
    draw_both(points)
    draw_histogram(nbrs)



