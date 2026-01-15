# Register Embedding Problem Solver

This project implements a **Greedy Heuristic** to solve the **Register Embedding Problem**. The goal is to find an optimal position of atoms while satisfying technical constraints and minimizing the difference between a target matrix and physical interaction.

## 📖 Problem Description

The problem consists of embedding $n$ atoms into a set of candidate positions in a 2D register. 
1. **Physical Model**: Interactions between atoms follow the Van-der-Waals potential, defined as $U(d) = C / d^6$, where $d$ is the Euclidean distance.
2. **Technical Constraints**:
   - **Minimal Distance**: No two atoms can be closer than `minDistanceBetweenAtoms`.
   - **Register Bounds**: All atoms must stay within `maxDistanceFromCenter` from the origin $(0,0)$.
3. **Objective**: Minimize the cost function $\sum |Q_{ij} - U_{ij}|$, where $Q$ is the target interaction matrix and $U$ is the resulting physical interaction matrix.

## 📂 Project Structure

The implementation is divided into three specialized files:

### 1. `Utils.py` (Constants & Math)
Contains the core physical parameters and utility functions.

### 2. `Instance.py` (Environment & Data)
Handles the generation of the search space and problem instances:
* **`layout` Class**: Generates random candidate points and identifies "incompatible" pairs (points that are too close to each other).
* **`dataInstance` Class**: Manages the target matrix $Q$. It includes a built-in 6-atom technical test case and a generator for random instances.

### 3. `Algorithm.py` (The Heuristic)
Contains the optimization logic:
* **`GreedyHeuristic` Class**: 
    - Filters out invalid layout points based on technical constraints.
    - Computes a "Gamma Set" of all possible atom-to-vertex pairings.
    - Sorts assignments by their local error $|Q_{ij} - U_{ij}|$ and greedily assigns atoms to the best available, compatible positions.

## 🚀 Usage

### Prerequisites
- Python 3.x
- NumPy

### Execution
Simply run the main algorithm file to see a demonstration of the solver using both a fixed technical test and a random instance:

```bash
python Algorithm.py
