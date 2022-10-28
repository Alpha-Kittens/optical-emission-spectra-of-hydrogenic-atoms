import numpy as np
poly_coeffs = 'abcdefg'


model_poly = lambda x_0 : lambda x, **params : poly(x, x_0, **params)
model_exp = lambda x_0 : lambda x, **params : exp(x, x_0, **params)


def poly(x, x_0, **params):
    coeffs = sorted(params.keys())
    p = 0
    for i in range(len(coeffs)):
        if coeffs[i] in poly_coeffs:
            p += params[coeffs[i]] * (x - x_0)**i
    return p

def exp(x, x_0, **params):
    return params['n'] * np.exp(params['r']*(x-x_0) - params['o'])

def poly_params(params, n):
    if n > 6:
        print ("polynomial capped at 6")
        n = 6
    for i in range(n + 1):
        if poly_coeffs[i] != 'b':
            params.add(poly_coeffs[i], val = 0)
        else:
            params.add(poly_coeffs[i], val = 1)

def exp_params(params, nval = 0.01, rval = 1/1000, oval = 0):
    params.add('n', val = nval)
    params.add('r', val = rval)
    params.add('o', val = oval)
    
def extract(params):
    vals = {}
    errs = {}
    for param in sorted(params.keys()):
        vals[param] = params[param].value
        errs[param] = params[param].stderr
    return vals, errs

