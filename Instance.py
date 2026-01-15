#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
@Author: Youcef Magnouche
Algorithm (Heuristic) for solving 'The register embedding problem'
@Date: 2026
"""

from Utils import *
import numpy as np

class layout:
    def __init__(self, nbVerticesInLayout):
        self.nbVerticesInLayout = nbVerticesInLayout
        self.layout = []
        self.incompatibilities = np.zeros((self.nbVerticesInLayout, self.nbVerticesInLayout)) # pairs of vertices in the layout do not respect the technical constraints.
        
        self.generateLayout(nbVerticesInLayout)
        self.computeIncompatibilities()

    # Return the layout coordinates of a given vertex
    def getLayoutVertexCoo(self, vertex):
        return self.layout[vertex]

    # Compute incompatibilities between vertices of the layout.
    # Two vertices are incompatible if the minimal distance is not respected.
    # Two vertices are incompatible if 
    def computeIncompatibilities(self):
        nbIncompatiblePairs = 0
        for vertex1 in range(self.nbVerticesInLayout):
            for vertex2 in range(vertex1 + 1, self.nbVerticesInLayout):
                if (euclidianDistance(self.getLayoutVertexCoo(vertex1), self.getLayoutVertexCoo(vertex2)) < minDistanceBetweenAtoms):
                    self.incompatibilities[vertex1][vertex2] = 1
                    self.incompatibilities[vertex2][vertex1] = 1
                    nbIncompatiblePairs += 1
                    continue
        print("Number of incompatible pairs: ", nbIncompatiblePairs, "/", self.nbVerticesInLayout * (self.nbVerticesInLayout - 1)/2)

    # Generate a random layout
    # nbVerticesInLayout: number vertices in the layout
    def generateLayout(self, nbVerticesInLayout):
        self.nbVerticesInLayout = nbVerticesInLayout

        # Generate nbVerticesInLayout different points
        points = set()
        i = 0
        while len(points) < nbVerticesInLayout:
            x = getRandomFloat(minCoordinates, maxCoordinates, 3)
            y = getRandomFloat(minCoordinates, maxCoordinates, 3)
            points.add((x, y))
            i += 1
            if i == 10 * nbVerticesInLayout:
                print("Error in generating layout")
                exit(-1)
        self.layout = list(points)

        nbPointFarFromCenter = 0
        for point in self.layout:
            if euclidianDistance(point, registerCenter) > maxDistanceFromCenter:
                nbPointFarFromCenter += 1
        print("Number of points far from the center of the register: ", nbPointFarFromCenter)

    # Remove bad vertices
    def getFilteredVertices(self, minDegree = 0):
        degreesOfVertices = [self.nbVerticesInLayout for _ in range(self.nbVerticesInLayout)]
        for vertex1 in range(self.nbVerticesInLayout):
            for vertex2 in range(vertex1 + 1, self.nbVerticesInLayout):
                if self.incompatibilities[vertex1][vertex2] == 1:
                    degreesOfVertices[vertex1] -= 1
                    degreesOfVertices[vertex2] -= 1
        
        # The set of vertices than can be in a clique of size nbAtoms
        filteredV = [v for v in range(self.nbVerticesInLayout) if degreesOfVertices[v] >= minDegree]

        # The set of vertices respecting technical constraints : "Maximum distance from the center of the register"
        filteredV = [v for v in filteredV if euclidianDistance(self.getLayoutVertexCoo(v), registerCenter) <= maxDistanceFromCenter]

        return filteredV

# dataInstance is a class representing the problem data
# it allows to generate and display the problem data
# Arguments:
class dataInstance:
    def __init__(self):
        # Init with the example given in the technical test
        self.nbAtoms = 6
        self.Q = np.array([[-17, 10, 10, 10, 0, 20], [10, -18, 10, 10, 10, 20], [10, 10, -29, 10, 20, 20],\
                           [10, 10, 10, -19, 10, 10], [0, 10, 20, 10, -17, 10], [20, 20, 20, 10, 10, -28]])

    # Display the Q matrix
    def display(self):
        print("Q = ", self.Q)
        print("Layout = ", self.layout)
        print("incompatibilities = ", self.incompatibilities)

    # Return the matrix entry corresponding to the given coordinates
    def getQValue(self, rowId, colId):
        if max(rowId, colId) >= self.nbAtoms or min(rowId, colId) < 0:
            print("Error in getQValue coordinates")
            exit(-1)
        return self.Q[rowId][colId]

    # Create a random Q matrix
    # nbAtoms: number of atoms used for the optimization (= number of variables)
    # minQValue: minimum value of each matrix entry
    # maxQValue: maximum value of each matrix entry
    def generateRandomInstance(self, nbAtoms, minQValue = 0, maxQValue = 10):
        self.nbAtoms = nbAtoms
        self.Q = np.random.randint(minQValue, maxQValue, size=(self.nbAtoms, self.nbAtoms))