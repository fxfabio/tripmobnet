import numpy as np

def MC_HiddenVariable_linkingP(Fd, alphabeta):
    alpha = alphabeta[0]
    beta = alphabeta[1]
    n = len(Fd)

    # Extract individual columns from Fd
    Fx = Fd[:, 0]
    Fy = Fd[:, 1]

    # Calculate the linking function
    p = (1 + alpha) * (1 + beta)
    linkF = (1 / p) * (Fx**alpha * Fy**beta) / (max(Fx)**alpha * max(Fy)**beta)

    # Alternatively, you can uncomment one of the following linkF definitions
    # to use different linking functions:

    # linkF = (1 / p) * (np.exp(Fx * alpha / c) * np.exp(Fy * beta / c)) / (max(np.exp(Fx * alpha / c)) * max(np.exp(Fy * beta / c)))
    # C = 400
    # linkF = erf((Fd**alpha * Fd.T**beta) / (max(Fd)**(alpha + beta)))
    # linkF = ((Fd**alpha * Fd.T**beta) / (max(Fd)**(alpha + beta))) / (1 + (Fd**alpha * Fd.T**beta) / (max(Fd)**(alpha + beta)))
    # linkF = (Fd / max(Fd)) * (Fd.T / max(Fd)) / (1 + c * (Fd / max(Fd)) * (Fd.T / max(Fd)))
    # C = 1
    # linkF = 1 / (1 + np.exp(-(Fd + Fd.T - C)))

    L = linkF
    return L

