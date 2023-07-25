import numpy as np
import networkx as nx

def MC_HiddenVariableGraph_main(N, T, multi):
    # Monte Carlo simulation of hidden variable mobility multi-graph model
    # with multilink procedure or not in string 'multi'.
    # @FVanni 2023

    # AD = 10  # average degree
    # N = 300  # number of nodes

    alpha = 1  # as in-visit exponent
    beta = 1  # as out-trip exponent
    alphabeta = [alpha, beta]

    mu_x = 2.001
    a = 0.001
    b = N - 1

    if multi.lower() in ['yes', 'y', 'si']:
        gx = 1
    elif multi.lower() in ['no', 'n', 'not']:
        gx = 0
    else:
        raise ValueError('Invalid string on multi')

    # Generation of the variable distribution for the locations
    xmin = a
    xmax = b
    Fy = MC_variable_distributions(N, 'exponential', 1/1000)
    Fx = MC_variable_distributions(N, 'xmin', xmin, 'powerlaw', mu_x)
    Fxy = np.column_stack((Fx, Fy))

    # Generation of the interbank Linking Function
    linkF = MC_HiddenVariable_linkingP(Fxy, alphabeta)

    n = len(Fx)
    ad = T / n

    totlink = round(ad * n)
    if totlink >= n * (n - 1):
        print('Warning: Multi-link procedure activated')

    t = 1
    linkFF = linkF.copy()  # self-loop included
    L = linkFF.flatten()
    Lcut = 1 - L
    lin_ind = np.zeros_like(L)
    
    tic = time.time()

    while t <= totlink:
        cL = np.insert(np.cumsum(L), 0, 0)
        ra = np.random.rand() * cL[-1]
        m = np.searchsorted(cL, ra)
        lin_ind[m - 1] += 1  # link multipli
        L[m - 1] = gx * L[m - 1]  # once the link is created, it cannot be repeated
        t += 1

    toc = time.time()
    print('Elapsed Time:', toc - tic)

    A = lin_ind.reshape(linkF.shape)
    # so that:
    # np.sum(A, axis=0)  # is the out-Degree distribution - Departures
    # np.sum(A, axis=1)  # is the in-Degree distribution - Visits

    return A, Fxy

# Rest of the code from the previous examples...
# ...

# Example usage
N = 100  # Number of nodes
T = 10  # Total number of links (edges) to generate
multi = 'no'  # Multi-link procedure, options: 'yes', 'no'

A, Fxy = MC_HiddenVariableGraph_main(N, T, multi)

# Convert the adjacency matrix A into a directed graph
G = nx.DiGraph(A)

# Add node attributes (locations) to the graph
for node, location in enumerate(Fxy):
    G.nodes[node]['location'] = location

# Example of how to access the out-degree and in-degree distributions
out_degrees = list(G.out_degree())
in_degrees = list(G.in_degree())

# Rest of the code...
# ...

