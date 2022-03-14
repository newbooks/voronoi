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

def draw_delaunay(points):
    tri = Delaunay(points)
    ps = np.array(points)
    fig = plt.figure()
    fig.set_size_inches(10, 10)
    plt.axis("equal")
    plt.triplot(ps[:,0], ps[:,1], tri.simplices)
    plt.plot(ps[:,0], ps[:,1], "o")
    plt.show()

def draw_both(points):


if __name__ == '__main__':
    square = [(0, 0), (0, 1), (1, 1), (1, 0)]
    draw_voronoi(square)
    draw_delaunay(square)