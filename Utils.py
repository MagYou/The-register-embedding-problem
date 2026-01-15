#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
@Author: Youcef Magnouche
Algorithm (Heuristic) for solving 'The register embedding problem'
@Date: 2026
"""

import os
import sys
import random
import numpy as np
import math
import matplotlib.pyplot as plt

# Constant and Initialization
np.set_printoptions(precision=1, suppress=True)
DEBUG = False
seed = 0
random.seed(0) # set the random seed to a fixed value.
np.random.seed(seed) # set the random seed to a fixed value.

C = 6.5 # constant in the Van-der-Waals interaction function, e.g 6.5 for hydrogen atom
minDistanceBetweenAtoms = 0.89 # from Technical Constraint: the distance between each two atoms
                               # must be greater than the given value
maxDistanceFromCenter = 1.5 # from Technical Constraint: the distance between the atom and the
                            # register center must be lower than the given value
minCoordinates = -2.0 # minCoordinates: minimum value of x and y coordinates
maxCoordinates = 2.0 # maxCoordinates: maximum value of each matrix entry
registerCenter = (0.0, 0.0) # coordinates of the register's center

# Euclidian distance between to given points
def euclidianDistance(point1, point2):
    diffX = point1[0] - point2[0]
    diffY = point1[1] - point2[1]
    return math.sqrt(diffX * diffX + diffY * diffY)

# Get Random float
def getRandomFloat(minV, maV, d=0):
    return round(random.uniform(minV, maV), d)

# U function
def U(point1, point2):
    if euclidianDistance(point1, point2) <= 1e-6:
        print("Warning: Distance between ", point1, " and ", point2, " is ", euclidianDistance(point1, point2))
        return 1e08
    return C / math.pow(euclidianDistance(point1, point2), 6)