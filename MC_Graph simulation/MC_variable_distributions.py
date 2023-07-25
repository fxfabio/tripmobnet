import numpy as np

def MC_variable_distributions(n, *args):
    # Default values
    xmin = 5
    xmax = 50
    mu = 2
    
    # ARGUMENTS selections
    type = 'PL'
    
    i = 0
    while i < len(args):
        argom = 1
        if isinstance(args[i], str):
            if args[i] == 'xmin':
                xmin = args[i + 1]
                i += 2
            elif args[i] == 'powerlaw':
                type = 'PL'
                mu = args[i + 1]
                i += 2
            elif args[i] == 'truncated':
                type = 'TC'
                mu = args[i + 1]
                xmax = args[i + 2]
                i += 3
            elif args[i] == 'powercutoff':
                type = 'PC'
                mu = args[i + 1]
                lambda_ = args[i + 2]
                i += 3
            elif args[i] == 'exponential':
                type = 'EX'
                lambda_ = args[i + 1]
                i += 2
            elif args[i] == 'stretched':
                type = 'ST'
                lambda_ = args[i + 1]
                bs = args[i + 2]
                i += 3
            else:
                argom = 0
                print('(scalefree) !!! INVALID INPUTS !!! , using default truncated pareto')
                i = 5
                n = 250
                xmin = 5
                xmax = 50
                mu = 2
                type = 'TC'
        
        if argom == 0:
            break
    
    if type == 'EX':
        # EXPONENTIAL (random graph)
        x = xmin - (1 / lambda_) * np.log(1 - np.random.rand(n, 1))
    elif type == 'TC':
        # TRUNCATED Pareto
        alfa = mu - 1
        u = np.random.rand(n, 1)
        l = xmin
        h = xmax
        if l < h:
            x = (-((u * h ** alfa - u * l ** alfa - h ** alfa) / (h ** alfa * l ** alfa)) ) ** (-1 / alfa)
        else:
            raise ValueError('Wrong ranges')
    elif type == 'ST':
        # stretched exponential truncation
        x = (xmin ** bs - (1 / lambda_) * np.log(1 - np.random.rand(n, 1))) ** (1 / bs)
    elif type == 'PC':
        # Exponential cutoff
        x = np.array([])
        y = xmin - (1 / lambda_) * np.log(1 - np.random.rand(10 * n, 1))
        while True:
            y = y[y >= (y / xmin) ** (-mu)]
            x = np.concatenate((x, y), axis=0)
            q = len(x) - n
            if q == 0:
                break
            if q > 0:
                r = np.random.permutation(len(x))
                x = np.delete(x, r[:q])
                break
            if q < 0:
                y = xmin - (1 / lambda_) * np.log(1 - np.random.rand(10 * n, 1))
    elif type == 'PL':
        # POWERLAW with xmin
        x = xmin * (1 - np.random.rand(n, 1)) ** (-1 / (mu - 1))
    else:
        x = xmin * (1 - np.random.rand(n, 1)) ** (-1 / (mu - 1))
    
    Fd = x
    
    return Fd

