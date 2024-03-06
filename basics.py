import numpy as np
from itertools import chain, combinations


def powerset(elements):
    '''
    Generate all possible subsets of a set of elements.
    '''
    all_subsets = []
    for subset_size in range(len(elements) + 1):
        all_subsets.extend(combinations(elements, subset_size))
    return [list(subset) for subset in all_subsets]

# https://stackoverflow.com/questions/952914/how-do-i-make-a-flat-list-out-of-a-list-of-lists
def flatten(xss):
    '''
    Used to flatten a list.
    '''
    return [x for xs in xss for x in xs]


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
        

    def convex_square(self, node1, node2):
        '''
        Given node1 <= node2, generate all points in the square with source node1 and sink node2.
        For example, convex_square((0,0), (1,1)) = [(0,0), (0,1), (1,0), (1,1)].
        Used to construct int_hull that returns the convex hull of an interval (sources, sinks).
        '''
        if self.is_smaller(node1, node2):
            points = []
            for coordinates in self.nodes:
                if all(node1[i] <= coordinates[i] <= node2[i] for i in range(len(self.dimensions))):
                    points.append(tuple(coordinates))
            return points 
        else:
            return []


    def int_hull(self, interval):
        '''
        Given an interval (sources and sinks), return the convex hull.
        '''
        points = []

        # strategy : add all points in convex_square(source, sink) for a given source and a given sink, then remove points if needed
        for src in interval.src:
            for snk in interval.snk :
                points = points + self.convex_square(src, snk) 
        points = list(set(points))

        for new_pt in points: # check if we need to remove points
            if self.to_be_removed(new_pt, interval):
                points.remove(new_pt)

        return points
    
    
    def is_corner(self, point, interval):
        '''
        Check if a point is a corner of an interval.
        '''
        for i in range (len(self.dimensions)):
            aligned = False
            for src in interval.src:
                if src[i] == point[i]:
                    aligned = True
            for snk in interval.snk:
                if snk[i] == point[i]:
                    aligned = True
            if not aligned:
                return False
        return True
    

    def cover(self, interval, conv=True):
        '''
        Return the cover of an interval.

        If conv == True, return intervals in the form of convex hulls.
        If conv == False, return intervals in the form (src, snk).

        Strategy: for all corner vertices, we check if the interval given when we add the vertex with corner[i] = +- 1 is
        indeed an interval.

        TO DO: this strategy works for 2D grids, check formally if this intuition scale to higher dimensions.
        '''
        cover = [] # list of intervals that are in Cov(I)
        hull = self.int_hull(interval)

        potential_int = hull.copy()
        for pt_idx in range(len(potential_int)):
            if self.is_corner(potential_int[pt_idx], interval):
                point = list(potential_int[pt_idx])
                # print(f"the corner considered is {point}")
                for i in range(len(self.dimensions)):
                    # print(f"i = {i}")
                    val = point[i]
                    if val != self.dimensions[i] - 1: 
                        # print("val is not max")
                        point[i] = val + 1 # first check with +1
                        # print(point)
                        potential_int.append(tuple(point))
                        if tuple(point) not in hull and self.is_convex(potential_int):
                            if conv:
                                cover.append(potential_int)
                            else:
                                tmp = self.get_src_snk(potential_int)
                                cover.append((tmp[0], tmp[1]))
                            # print("append")
                        # else:
                            # print("not append")
                        potential_int = hull.copy()
                        point[i] = val
                    # else:
                        # print("val is max")
                    if val != 0:
                        # print("val is not zero")
                        point[i] = val - 1 # then check with -1
                        potential_int.append(tuple(point))
                        if tuple(point) not in hull and self.is_convex(potential_int):
                            if conv:
                                cover.append(potential_int)
                            else:
                                tmp = self.get_src_snk(potential_int)
                                cover.append((tmp[0], tmp[1]))
                            # print("append")
                        # else:
                            # print("not append")
                        potential_int = hull.copy()
                        point[i] = val
                    # else:
                        # print("val is zero")
                          
        return cover


    def is_convex(self, points):
        '''
        Return True if points = int_hull(interval) for some interval, False otherwise. 
        '''
        tmp = self.get_src_snk(points)
        interval = Interval(tmp[0], tmp[1])
        if set(points) == set(self.int_hull(interval)):
            return True
        else:
            return False


    def get_src_snk(self, points):
        '''
        From a list of points, which is assumed to be equal to int_hull(interval) for some interval, return (sources, sinks).
        '''
        src = []
        snk = []
        for pt in points:
            if self.add_source(pt, points):
                src.append(pt)
            if self.add_sink(pt, points):
                snk.append(pt)
        return (src, snk)  


    def add_source(self, pt, points):
        '''
        Attribute that is helpful for constructing get_src_snk.
        Check if a potential point is a source.
        '''
        for other in points:
            if self.is_smaller(other, pt) and other != pt:
                return False
        return True
    

    def add_sink(self, pt, points):
        '''
        Attribute that is helpful for constructing get_src_snk.
        Check if a potential point is a sink.
        '''
        for other in points:
            if self.is_smaller(pt, other) and other != pt:
                return False
        return True


    def to_be_removed(self, point, interval):
        '''
        Attribute that is helpful for constructing int_hull.
        Check if a potential point should be removed, in the construction of the convex hull.
        '''
        for src in interval.src:
            if self.is_smaller(point, src) and point != src:
                return True
        for snk in interval.snk:
            if self.is_smaller(snk, point) and point != snk:
                return True
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
                matrix_N[r:r+self.vecs[self.meet(interval.snk[i], interval.snk[j])], idx_col[i]:idx_col[i+1]] = self.evaluation(self.meet(interval.snk[i], interval.snk[j]), interval.snk[i]).T
                # second block -M_{a_{i,j}, a_j}
                matrix_N[r:r+self.vecs[self.meet(interval.snk[i], interval.snk[j])], idx_col[j]:idx_col[j+1]] = -self.evaluation(self.meet(interval.snk[i], interval.snk[j]), interval.snk[j]).T
                r += self.vecs[self.meet(interval.snk[i], interval.snk[j])]
        
        return matrix_N.T
    

    def int_rank(self, interval):
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


    def int_replacement(self, interval):
        '''
        Return the interval replacement of an interval. Uses the formula with the cover of the interval.
        '''
        repl = self.int_rank(interval) # corresponds to the empty set in the sum
        cov_ps = powerset(self.cover(interval))

        # compute V S
        for c in cov_ps[1::]:
            if len(c) == 1:
                tmp = self.get_src_snk(c[0])
                interval = Interval(tmp[0], tmp[1])
                repl = repl - self.int_rank(interval)
            if len(c) > 1:
                eps = len(c)
                c = [list(set(flatten(c)))]
                tmp = self.get_src_snk(c[0])
                interval = Interval(tmp[0], tmp[1])
                if eps %2 == 0:
                    repl = repl + self.int_rank(interval)
                else:
                    repl = repl - self.int_rank(interval)
        return repl


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

