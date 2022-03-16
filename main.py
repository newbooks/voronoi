#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi, voronoi_plot_2d
from scipy.spatial import Delaunay


def draw_voronoi(points):
    vor = Voronoi(points)
    fig = voronoi_plot_2d(vor, show_points=True, show_vertices=True)
    fig.set_size_inches(10, 10)
    plt.axis("equal")

    plt.show()
    return

def draw_delaunay(points):
    tri = Delaunay(points)
    ps = np.array(points)
    fig = plt.figure()
    fig.set_size_inches(10, 10)
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

    return


def entropy_from_histogram(his):

    return

def random_points(n=100, r=0.1):

    return


if __name__ == '__main__':
    square = [(0, 0), (0, 1), (1, 1), (1, 0)]
#    draw_voronoi(square)
    #draw_delaunay(square)
    #draw_both(square)
    print(neighbors_vor(square))