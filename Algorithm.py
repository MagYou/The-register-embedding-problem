#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
@Author: Youcef Magnouche
Algorithm (Heuristic) for solving 'The register embedding problem'
@Date: 2026
"""

from Instance import *
import numpy as np

# Heuristic algorithm for the optimization problem
class GreedyHeuristic:

    def __init__(self, data, layout):
        self.data = data
        self.layout = layout

    # Compute gamma set
    def computeGammaSet(self, V):
        gamma = []
        for vertex1 in V:
            for vertex2 in V:
                if vertex1 != vertex2:
                    for atom1 in range(self.data.nbAtoms):
                        for atom2 in range(atom1 + 1, self.data.nbAtoms):
                            gamma.append([(atom1, vertex1), (atom2, vertex2)])
        return gamma

    def solve(self):
        print("Optimization starts.")

        filteredV = self.layout.getFilteredVertices(minDegree = (self.data.nbAtoms - 1))
        gamma = self.computeGammaSet(filteredV)

        if DEBUG:
            print("Size of gamma : ", len(gamma))

        # Compute the weight = |Q_(a_i a_j) - U_e|
        w = [0] * len(gamma)
        for element in range(len(gamma)):
            atom1, vertex1 = gamma[element][0]
            atom2, vertex2 = gamma[element][1]
            coorVertex1 = self.layout.getLayoutVertexCoo(vertex1)
            coorVertex2 = self.layout.getLayoutVertexCoo(vertex2)
            w[element] = abs(self.data.getQValue(atom1, atom2) - U(coorVertex1, coorVertex2))

        if len(gamma) < float((self.data.nbAtoms) * (self.data.nbAtoms - 1))/2:
            return [], 1e08

        # Sort the element of gamma based on the computed weight
        # (Suite de la logique de tri et d'assignation)
        sortedW, sortedGamma = zip(*sorted(zip(w, gamma)))

        # Assign atoms to vertices
        atom_vertex = [-1 for _ in range(self.data.nbAtoms)]
        vertex_atom = [-1 for _ in range(self.layout.nbVerticesInLayout)]
        V_sol = []
        cost = 0

        for element in range(len(sortedGamma)):
            atom1, vertex1 = sortedGamma[element][0]
            atom2, vertex2 = sortedGamma[element][1]

            if (atom_vertex[atom1] == -1 and len(V_sol) < self.data.nbAtoms and vertex_atom[vertex1] == -1) or (atom_vertex[atom1] == vertex1 and vertex_atom[vertex1] == atom1):
                if (atom_vertex[atom2] == -1 and len(V_sol) < self.data.nbAtoms and vertex_atom[vertex2] == -1) or (atom_vertex[atom2] == vertex2 and vertex_atom[vertex2] == atom2):
                    if self.layout.incompatibilities[vertex1][vertex2] == 0:
                        if all(self.layout.incompatibilities[vertex1][v2] == 0 for v2 in V_sol):
                            if all(self.layout.incompatibilities[v1][vertex2] == 0 for v1 in V_sol):
                                atom_vertex[atom1], atom_vertex[atom2] = vertex1, vertex2
                                vertex_atom[vertex1], vertex_atom[vertex2] = atom1, atom2
                                V_sol.extend([vertex1, vertex2])
                                V_sol = list(set(V_sol))
                                cost += sortedW[element]

        if -1 in atom_vertex:
            cost += 1e09

        if DEBUG:
            print("assignement ", atom_vertex, "\n")
            print("V_sol ", V_sol, "\n")

        U_matrix = np.zeros((self.data.nbAtoms, self.data.nbAtoms))
        for atom1 in range(self.data.nbAtoms):
            for atom2 in range(self.data.nbAtoms):
                if atom1 != atom2:
                    U_matrix[atom1][atom2] = U(self.layout.getLayoutVertexCoo(atom_vertex[atom1]), self.layout.getLayoutVertexCoo(atom_vertex[atom2]))

        print("Q = \n", self.data.Q)
        print("U = \n", U_matrix)
        print("Objective value for nbVerticesInLayout = ", self.layout.nbVerticesInLayout, " is ", cost, "\n")
        print("Sol ", atom_vertex)

        self.checkTechnicalConstraints(atom_vertex)
        print("Optimization done.\n\n\n")
        return atom_vertex, cost

    # Check if the given solution satisfies the technical constraints
    def checkTechnicalConstraints(self, solution):
        # Check the number of pairs violating the minimal distance
        nbIncompatiblePairs = 0
        for v1 in range(len(solution)):
            for v2 in range(v1 + 1, len(solution)):
                vertex1 = solution[v1]
                vertex2 = solution[v2]
                if vertex1 != -1 and vertex2 != -1 and self.layout.incompatibilities[vertex1][vertex2] == 1:
                    nbIncompatiblePairs += 1
        print("Number of incompatible pairs: ", nbIncompatiblePairs)

        # Check the number of points far from the center of the register
        nbPointsFarFromCenter = 0
        for point in solution:
            if point != -1 and euclidianDistance(self.layout.getLayoutVertexCoo(point), registerCenter) > maxDistanceFromCenter:
                nbPointsFarFromCenter += 1
        print("Number of points far from the center ", nbPointsFarFromCenter)

# --- Bloc Main ---
if __name__ == '__main__':
    # Given matrix in the test
    print("Technical Test Example")
    problem = dataInstance() 
    lay = layout(500) # Test with q vertices in the layout
    algorithm = GreedyHeuristic(problem, lay)
    algorithm.solve()

    # Random test
    print("Random generated instance")
    problem.generateRandomInstance(4) # Test with n atoms
    lay = layout(300) # Test with q vertices in the layout
    algorithm = GreedyHeuristic(problem, lay)
    algorithm.solve()