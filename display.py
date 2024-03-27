from utils import Representation, Interval

def display_rep(rep):
        """
        - Input: representation
        Display the grid with vector space dimensions for each node, 
        and arrows between each pair of adjacent nodes iff there is a matrix between these two nodes.
        """
        if len(rep.dimensions) > 2:
            print(f"Error: The vizualization is only available for 1D or 2D grids, but the grid has dimension {len(rep.dimensions)}.")

        elif len(rep.dimensions) == 2:
            for i in range(rep.dimensions[0]):
                for j in range(rep.dimensions[1]):
                    node = (i, j)
                    vecs_dim = rep.vecs.get(node, None)
                    edge = (node, (i, j+1))
                    if vecs_dim is not None:
                        if j != rep.dimensions[1] - 1:
                            if edge in rep.mats:
                                print(f"{vecs_dim} -------> ", end="")
                            else:
                                print(f"{vecs_dim}          ", end="")
                        else:
                            print(f"{vecs_dim}", end="")
                    else:
                        if j != rep.dimensions[1] - 1:
                            if edge in rep.mats:
                                print(f". -------> ", end="")
                            else:
                                print(f".          ", end="")
                        else:
                            print(".", end="")

                print()

                if i != rep.dimensions[0] - 1:
                    for j in range(rep.dimensions[1]):
                        edge = ((i, j), (i+1, j))
                        if edge in rep.mats:
                            print("|          ", end="")
                        else:
                            print("           ", end="")
                    print()
                    for j in range(rep.dimensions[1]):
                        edge = ((i, j), (i+1, j))
                        if edge in rep.mats:
                            print("|          ", end="")
                        else:
                            print("           ", end="")
                    print()
                    for j in range(rep.dimensions[1]):
                        edge = ((i, j), (i+1, j))
                        if edge in rep.mats:
                            print("v          ", end="")
                        else:
                            print("           ", end="")
                    print()

        elif len(rep.dimensions) == 1:
            for j in range(rep.dimensions[0]):
                node = (j,)
                vecs_dim = rep.vecs.get(node, None)
                edge = (node, (j+1,))
                if vecs_dim is not None:
                    if j != rep.dimensions[0] - 1:
                        if edge in rep.mats:
                            print(f"{vecs_dim} -------> ", end="")
                        else:
                            print(f"{vecs_dim}          ", end="")
                    else:
                        print(f"{vecs_dim}", end="")
                else:
                    if j != rep.dimensions[0] - 1:
                        if edge in rep.mats:
                            print(f". -------> ", end="")
                        else:
                            print(f".          ", end="")
                    else:
                        print(".", end="")
            print()

def display_interval(rep, itv):
    """
    Input:
    - Representation rep
    - Interval itv
    Display the grid and arrows between each pair of adjacent nodes contained in the interval.
    """
    if len(rep.dimensions) > 2:
        print(f"Error: The vizualization is only available for 1D or 2D grids, but the grid has dimension {len(rep.dimensions)}.")

    elif len(rep.dimensions) == 2:
        hull = rep.int_hull(itv)
        for i in range(rep.dimensions[0]):
            for j in range(rep.dimensions[1]):
                if j != rep.dimensions[1] - 1:
                    if (i,j) in hull and (i,j+1) in hull:
                        print(f"X -------> ", end="")
                    elif (i,j) in hull and (i,j+1) not in hull:
                        print(f"X          ", end="")
                    else: 
                        print(f".          ", end="")
                else:
                    if (i,j) in hull:
                        print("X", end="")
                    else:
                        print(".", end="")
            print()
            if i != rep.dimensions[0] - 1:
                for j in range(rep.dimensions[1]):
                    if (i,j) in hull and (i+1,j) in hull:
                        print("|          ", end="")
                    else:
                        print("           ", end="")
                print()
                for j in range(rep.dimensions[1]):
                    if (i,j) in hull and (i+1,j) in hull:
                        print("|          ", end="")
                    else:
                        print("           ", end="")
                print()
                for j in range(rep.dimensions[1]):
                    if (i,j) in hull and (i+1,j) in hull:
                        print("v          ", end="")
                    else:
                        print("           ", end="")
                print()

    elif len(rep.dimensions) == 1:
        hull = rep.int_hull(itv)
        for j in range(rep.dimensions[0] - 1):
            if (j,) in hull and ((j+1,)) in hull:
                print("X -------> ", end="")
            elif (j,) in hull and ((j+1,)) not in hull:
                print("X          ", end="")
            else:
                print(".          ", end="")
        else:
            if (j,) in hull:
                print("X", end="")
            else:
                print(".", end="")
        print()