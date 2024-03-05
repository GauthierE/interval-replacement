import numpy as np


class Interval:
    '''
    Define intervals by their sources and their sinks.
    '''
    def __init__(self, src, snk):
        self.src = src # list of nodes 
        self.snk = snk # list of nodes


class Representation:
    def __init__(self, dimensions):
        self.dimensions = dimensions # dimensions of the grid, ex: [2,3] for a 2D grid with 2 rows and 3 columns
        self.nodes = self.generate_nodes() # nodes are tuples
        self.edges = self.generate_edges() # edges are tuples of two nodes
        self.vecs = {} # vecs[node] = m, when the node is associated with the vector space k^m
        self.mats = {}  # mats[(node1, node2)] = np.array(), when there's an edge between node1 and node2


    def generate_nodes(self):
        nodes = []
        for coordinates in self.generate_coordinates():
            nodes.append(tuple(coordinates))
        return nodes


    def generate_coordinates(self, current_dim=0, current_coords=[]):
        if current_dim == len(self.dimensions):
            return [tuple(current_coords)]

        coordinates = []
        for i in range(self.dimensions[current_dim]):
            new_coords = current_coords + [i]
            coordinates.extend(self.generate_coordinates(current_dim + 1, new_coords))

        return coordinates


    def generate_edges(self):
        edges = []
        for node1 in self.nodes:
            for node2 in self.nodes:
                if all(node1[i] <= node2[i] for i in range(len(self.dimensions))):
                    edges.append((node1, node2))
        return edges
    

    def create_vecs(self, node, m): # node is the vector space k^m
        if node in self.nodes:
            self.vecs[node] = m
        else:
            print("Error: Trying to access a node which is not in the grid.")


    def create_matrix(self, node1, node2, matrix): # create a matrix from node1 to node2
        if (node1, node2) in self.edges:
            if matrix is None:
                 self.mats[(node1, node2)] = matrix
            else:
                u = self.vecs.get(node1)
                v = self.vecs.get(node2)
                if u is not None and v is not None and matrix.shape == (v, u): # check if the dimensions of the matrix match the vector spaces of the nodes
                    self.mats[(node1, node2)] = matrix
                else:
                    print("Error: Matrix dimensions do not match the vector spaces of the nodes.", end="\n")
                    print("Make sure Matrix has nonzero dimensions. Example M: k^2 -> k^1 should have shape (2,1) and not (2,).")
        else:
            print("Error: Nodes are not connected by an edge.")


    def join(self, node1, node2):
        '''
        Return the join between two nodes.
        '''
        return tuple(max(coord1, coord2) for coord1, coord2 in zip(node1, node2))


    def meet(self, node1, node2):
        '''
        Return the meet between two nodes.
        '''
        return tuple(min(coord1, coord2) for coord1, coord2 in zip(node1, node2))
    

    def is_smaller(self, node1, node2):
        '''
        Check if there is a path from node1 to node2, ie if node1 <= node2.
        '''
        if all(node1[i] <= node2[i] for i in range(len(self.dimensions))):
            return True
        else:
            return False
    

    def evaluation(self, node1, node2):
        '''
        Return the value of the map M(p) where p is a path from node1 to node2.
        We assume that all linear maps between two adjacent nodes have been specified by the user.
        '''
        # will be multiplied by matrices encountered along a path p: node1 -> node2
        # we start from node2 (backward)
        mat = np.eye(self.vecs[node2])
        # list of coordinates of the current position on the path
        current = list(node2) 
        # coordinates of the adjacent node just before current
        previous = list(node2) # initialized equal to node2

        # since all the paths from node1 to node2 lead to the same evaluation (because of commutativity relations),
        # we may take the path where we increment the first dimension first, then the second dimension, etc.
        for i in range(len(self.dimensions)):
            k = current[i] - node1[i]
            while(k > 0):
                # update previous
                previous[i] = previous[i] - 1
                # multiply matrix
                if self.mats[(tuple(previous), tuple(current))] is not None:
                    mat = np.dot(mat, self.mats[(tuple(previous), tuple(current))])
                else:
                    return np.zeros((self.vecs[node2], self.vecs[node1]))
                # update current
                current[i] = current[i] - 1
                k = k-1

        return mat
    

    def matrix_M(self, interval):
        '''
        Given an interval with n sources, return the matrix_M.
        '''
        n = len(interval.src)
        # compute number of columns in matrix_M
        N = 0
        for i in range(n):
            N += self.vecs[interval.src[i]]
        # compute number of rows in matrix_M
        M = 0
        for i in range(n):
            for j in range(i+1,n):
                M += self.vecs[self.join(interval.src[i],interval.src[j])]

        # # FOR INTERVALS WITH 2 SOURCES ONLY - used for debugging
        # # first block
        # matrix_M[0:self.vecs[self.join(interval.src[0],interval.src[1])], 0:self.vecs[interval.src[0]]] = self.evaluation(interval.src[0], self.join(interval.src[0], interval.src[1]))
        # # second block
        # matrix_M[0:self.vecs[self.join(interval.src[0],interval.src[1])], self.vecs[interval.src[0]]:self.vecs[interval.src[0]]+self.vecs[interval.src[1]]] = -self.evaluation(interval.src[1], self.join(interval.src[0], interval.src[1]))
        
        # columns indices, to handle block matrices
        idx_col = [0]
        for i in range(n):
            idx_col.append(self.vecs[interval.src[i]] + idx_col[-1])

        r = 0 # current row

        # construct matrix_M
        matrix_M = np.zeros((M, N))
        for i in range(n):
            for j in range(i+1,n):
                # first block M_{a_{i,j}, a_i}
                matrix_M[r:r+self.vecs[self.join(interval.src[i], interval.src[j])], idx_col[i]:idx_col[i+1]] = self.evaluation(interval.src[i], self.join(interval.src[i], interval.src[j]))
                # second block -M_{a_{i,j}, a_j}
                matrix_M[r:r+self.vecs[self.join(interval.src[i], interval.src[j])], idx_col[j]:idx_col[j+1]] = -self.evaluation(interval.src[j], self.join(interval.src[i], interval.src[j]))
                r += self.vecs[self.join(interval.src[i], interval.src[j])]
        
        return matrix_M
    

    def matrix_N(self, interval):
        '''
        Given an interval with n sinks, return the matrix_N.
        '''
        n = len(interval.snk)
        # We do as in the construction of matrix_N, ie we construct the transpose of matrix_N
        # We just need to replace joins by meets
        N = 0
        for i in range(n):
            N += self.vecs[interval.snk[i]]
        M = 0
        for i in range(n):
            for j in range(i+1,n):
                M += self.vecs[self.meet(interval.snk[i],interval.snk[j])]
 
        idx_col = [0]
        for i in range(n):
            idx_col.append(self.vecs[interval.snk[i]] + idx_col[-1])

        r = 0 

        matrix_N = np.zeros((M, N))
        for i in range(n):
            for j in range(i+1,n):
                # first block M_{a_{i,j}, a_i}
                matrix_N[r:r+self.vecs[self.meet(interval.snk[i], interval.snk[j])], idx_col[i]:idx_col[i+1]] = self.evaluation(self.meet(interval.snk[i], interval.snk[j]), interval.snk[i])
                # second block -M_{a_{i,j}, a_j}
                matrix_N[r:r+self.vecs[self.meet(interval.snk[i], interval.snk[j])], idx_col[j]:idx_col[j+1]] = -self.evaluation(self.meet(interval.snk[i], interval.snk[j]), interval.snk[j])
                r += self.vecs[self.meet(interval.snk[i], interval.snk[j])]
        
        return matrix_N.T
    

    def interval_rank(self, interval):
        '''
        Given an interval, compute the interval rank.
        '''
        M = self.matrix_M(interval)
        N = self.matrix_N(interval)
        mat = self.evaluation(interval.src[0],interval.snk[0])

        # construct the bottom left block matrix
        C = np.zeros((N.shape[0],M.shape[1]))
        C[:mat.shape[0],:mat.shape[1]] = mat

        block = np.block([
            [M, np.zeros((M.shape[0], N.shape[1]))],
            [C, N]
            ])
        
        # this is for computing the rank in R. Change the field here if needed.
        if 0 in N.shape and 0 in M.shape: # rectangle
            return np.linalg.matrix_rank(mat) if 0 not in mat.shape else 0
        elif 0 in N.shape: # (n,1)-type
            return np.linalg.matrix_rank(block) - np.linalg.matrix_rank(M)
        elif 0 in M.shape: # (1,n)-type
            return np.linalg.matrix_rank(block) - np.linalg.matrix_rank(N)
        else:
            return np.linalg.matrix_rank(block) - np.linalg.matrix_rank(M) - np.linalg.matrix_rank(N)


    def check_commutativity(self):
        '''
        Check if all commutativity relations hold. To be implemented.
        '''
        return True


    def elements(self):
        '''
        Print all the elements of the grid:
        - nodes
        - edges
        - vector spaces associated to nodes
        - matrices associated to edges
        '''
        print("Nodes:", self.nodes)
        print("Edges:", self.edges)
        print("Vector spaces:", self.vecs)
        print("Matrices:", self.mats)


    def display(self):
        """
        Display the grid with vector space dimensions for each node, 
        and arrows between each pair of adjacent nodes iff there is a matrix between these two nodes.
        """
        if len(self.dimensions) > 2:
            print(f"Error: The vizualization is only available for 1D or 2D grids, but the grid has dimension {len(self.dimensions)}.")

        elif len(self.dimensions) == 2:
            for i in range(self.dimensions[0]):
                for j in range(self.dimensions[1]):
                    node = (i, j)
                    vecs_dim = self.vecs.get(node, None)
                    edge = (node, (i, j+1))
                    if vecs_dim is not None:
                        if j != self.dimensions[1] - 1:
                            if edge in self.mats:
                                print(f"{vecs_dim} -------> ", end="")
                            else:
                                print(f"{vecs_dim}          ", end="")
                        else:
                            print(f"{vecs_dim}", end="")
                    else:
                        if j != self.dimensions[1] - 1:
                            if edge in self.mats:
                                print(f". -------> ", end="")
                            else:
                                print(f".          ", end="")
                        else:
                            print(".", end="")

                print()

                if i != self.dimensions[0] - 1:
                    for j in range(self.dimensions[1]):
                        edge = ((i, j), (i+1, j))
                        if edge in self.mats:
                            print("|          ", end="")
                        else:
                            print("           ", end="")
                    print()
                    for j in range(self.dimensions[1]):
                        edge = ((i, j), (i+1, j))
                        if edge in self.mats:
                            print("|          ", end="")
                        else:
                            print("           ", end="")
                    print()
                    for j in range(self.dimensions[1]):
                        edge = ((i, j), (i+1, j))
                        if edge in self.mats:
                            print("v          ", end="")
                        else:
                            print("           ", end="")
                    print()

        elif len(self.dimensions) == 1:
            for j in range(self.dimensions[0]):
                node = (j,)
                vecs_dim = self.vecs.get(node, None)
                edge = (node, (j+1,))
                if vecs_dim is not None:
                    if j != self.dimensions[0] - 1:
                        if edge in self.mats:
                            print(f"{vecs_dim} -------> ", end="")
                        else:
                            print(f"{vecs_dim}          ", end="")
                    else:
                        print(f"{vecs_dim}", end="")
                else:
                    if j != self.dimensions[0] - 1:
                        if edge in self.mats:
                            print(f". -------> ", end="")
                        else:
                            print(f".          ", end="")
                    else:
                        print(".", end="")
            print()

    def display_interval(self, interval):
        """
        Display the grid and arrows between each pair of adjacent nodes contained in the interval.
        """
        if len(self.dimensions) > 2:
            print(f"Error: The vizualization is only available for 1D or 2D grids, but the grid has dimension {len(self.dimensions)}.")

        elif len(self.dimensions) == 2:
            pass # to be completed

        elif len(self.dimensions) == 1: # the interval has 1 source and 1 sink 
            for j in range(self.dimensions[0]):
                if j >= interval.src[0][0] and j < interval.snk[0][0]:
                    print(". -------> ", end="")
                else:
                    print(".          ", end="")
            print()

