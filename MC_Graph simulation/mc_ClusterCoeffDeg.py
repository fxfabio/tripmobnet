import numpy as np

def mc_ClusterCoeffDeg(A, binnum, author='Fagiolo'):
    # Cluster coefficient spectrum vs degree for weighted and directed mobility network
    # Input: 
    #       A = weighted adjacency matrix
    #       binnum = number of bins of the degrees (resolution of the spectrum)
    #       author = the type of cluster coefficient definitions, 3 string:
    #            'Fagiolo' - In-CC as in Fagiolo (2007)
    #            'Clemente' - In-CC as in Clemente-Grassi (2018)
    #            'Fardet' - FanIn-CC as in Fardet-Levina (2021)
    #
    #   Note: All weights must be between 0 and 1.
    # Output:
    #       k_degree = binned strength (degree) grouped in binnum number of bins
    #       C_cluster = Cluster coefficient for nodes of a certain strength (normalized to 1)
    #
    # @F.Vanni2023 

    if author.lower() not in ['fagiolo', 'clemente', 'fardet']:
        raise ValueError('Invalid author string. Options: "Fagiolo", "Clemente", "Fardet"')

    W = np.transpose(A) / np.max(np.abs(A))  # normalized scale by maximal weight

    # (CRin1) Fagiolo 2007 or (CRin2) Clemente 2017 or (CRin3) Fardet 2021
    Ctot, Cin, Kin, CRin1, CRin2, CRin3 = mc_clustering_coef(W)

    if author.lower() == 'fagiolo':
        CCin = CRin1
    elif author.lower() == 'clemente':
        CCin = CRin2
    elif author.lower() == 'fardet':
        CCin = CRin3
    else:
        CCin = Ctot  # Total Cluster coefficient

    CCin[~np.isfinite(CCin)] = np.nan
    C = Cin
    k_degree = np.sum(A, axis=0)  # In-degree

    x = np.unique(k_degree)
    k_bin = np.linspace(np.min(x), np.max(x), binnum)

    cc = np.zeros(len(k_bin))
    ccx = np.zeros(len(k_bin))

    for i in range(len(k_bin) - 1):
        ix = np.logical_and(k_degree >= k_bin[i], k_degree < k_bin[i + 1])
        cc[i] = np.nansum(C[ix])
        ccx[i] = np.nansum(CCin[ix])

    k_degree = k_bin
    C_cluster = cc
    C_cluster[C_cluster == 0] = np.nan

    C_clusterx = ccx
    C_clusterx[C_clusterx == 0] = np.nan  # Excluding nodes with clustering zero
    C_clusterx[~np.isfinite(C_clusterx)] = np.nan

    return k_degree, C_cluster, C_clusterx


def mc_clustering_coef(W):
    # Clustering coefficient for weighted directed connection matrix
    A = W != 0  # Binary adjacency matrix
    K = np.sum(A + np.transpose(A), axis=1)  # Total degree (in + out)
    Kin = np.sum(A, axis=1)  # In-degree

    Sin = (np.transpose(W) ** (1 / 3)) * (W ** (1 / 3)) ** 2  # a-la Fagiolo
    cyc3in = np.diag(Sin)
    KKin = Kin
    KKin[cyc3in == 0] = np.inf

    CYC3in = KKin * (KKin - 1)
    Cin = cyc3in / CYC3in  # In-CC (Fagiolo 2007)

    S = W ** (1 / 3) + (np.transpose(W) ** (1 / 3))  # Symmetrized weights matrix ^1/3
    cyc3 = np.diag(S ** 3) / 2  # Number of 3-cycles (directed triangles)
    K[cyc3 == 0] = np.inf  # If no 3-cycles exist, make C=0 (via K=inf)

    CYC3 = K * (K - 1) - 2 * np.diag(A ** 2)  # Number of all possible 3-cycles
    C = cyc3 / CYC3  # Total clustering coefficient

    Ksin = np.sum(W, axis=1)
    W_hat = W ** (1 / 3)
    triaIn_1 = np.diag(np.transpose(W_hat) @ W_hat ** 2)
    Kin1 = Kin
    Kin1[triaIn_1 == 0] = np.inf
    triaTot_1 = Kin1 * (Kin1 - 1)
    CR_in1 = triaIn_1 / triaTot_1  # Fagiolo

    triaIn_2 = 0.5 * np.diag(np.transpose(W) @ (A + np.transpose(A)) @ A)
    Kin2 = Kin
    Kin2[triaIn_2 == 0] = np.inf
    Ksin2 = Ksin
    Ksin2[triaIn_2 == 0] = np.inf
    triaTot_2 = Ksin2 * (Kin2 - 1)
    CR_in2 = triaIn_2 / triaTot_2  # Clemente

    Ksin312 = np.sum(W ** 0.5, axis=1)
    triaIn_3 = np.diag(np.transpose(W ** (2 / 3)) @ ((W ** (2 / 3)) ** 2))
    Kin3 = Kin
    Kin3[triaIn_3 == 0] = np.inf
    Ksin3 = Ksin
    Ksin3[triaIn_3 == 0] = np.inf
    Ksin312[triaIn_3 == 0] = np.inf
    triaTot_3 = (Ksin312 ** 2 - Ksin3)
    CR_in3 = triaIn_3 / triaTot_3  # Fardet 2021

    return C, Cin, Kin, CR_in1, CR_in2, CR_in3


# Rest of the code...
# ...

# Example usage
# Replace A with your weighted adjacency matrix (numpy array)
# binnum: the number of bins for calculating the cluster coefficient spectrum vs degree
# author: 'Fagiolo', 'Clemente', or 'Fardet' for the desired type of cluster coefficient definition
k_degree, C_cluster, C_clusterx = mc_ClusterCoeffDeg(A, binnum, author)

# Rest of the code...
# ...

