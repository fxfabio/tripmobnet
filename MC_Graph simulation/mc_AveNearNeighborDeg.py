import numpy as np

def mc_AveNearNeighborDeg(Adj, bin, inout='total'):
    # Gives a list of k-mean degree connectivity for the network for successive k=0,1,2, ...
    # The average degree connectivity is also known as the average nearest neighbor degree.
    # The k-average degree connectivity is the average of the mean neighbor degrees of vertices of degree k.
    # Input:
    # Adj - adjacency matrix (also weighted and directed)
    # bin - subdivision of the degree interval on which to find the average connectivity
    # inout - choose if "in", "out" or "total"

    if inout.lower() == 'in':
        Adj = Adj.T  # Transpose the adjacency matrix for "in" case
        ave = 'in'
    elif inout.lower() == 'out':
        ave = 'in'
    else:
        ave = 'total'

    deg = np.sum(Adj, axis=1)
    k_bin = np.linspace(0, np.max(deg), bin)
    aveK = np.zeros(bin - 1)
    k = deg
    knn = aveNeighborDegWW(Adj, 'in')

    for i in range(bin - 1):
        ix = np.logical_and(k > k_bin[i], k <= k_bin[i + 1])
        aveK[i] = np.mean(knn[ix])

    return aveK, k_bin


def aveNeighborDegWW(adj, inout='total'):
    # Average weighted and directed neighbor degree

    ave_n_deg = np.zeros(len(adj))  # initialize output vector

    adj = np.double(adj > 0)  # adjacency matrix
    W = adj  # weight matrix

    deg = np.sum(adj, axis=1)  # calculate the degree of each node

    for i in range(len(adj)):  # across all nodes
        neigh = np.where(adj[i] == 1)[0]  # neighbors of i, one link away
        peso = W[i, neigh]
        if len(neigh) == 0:
            ave_n_deg[i] = 0
            continue
        stren = np.sum(peso)
        ave_n_deg[i] = np.sum(deg[neigh] * peso) / stren

    return ave_n_deg

# Rest of the code...
# ...

# Example usage
# Replace Adj with your adjacency matrix (numpy array)
# bin: the number of bins for calculating the average connectivity
# inout: 'in', 'out', or 'total' for the desired type of neighbor degree
aveK, k_bin = mc_AveNearNeighborDeg(Adj, bin, inout)

# Rest of the code...
# ...

