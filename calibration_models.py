import numpy as np
poly_coeffs = 'abcdefg'


model_poly = lambda x_0 : lambda x, **params : poly(x, x_0, **params)
model_exp = lambda x_0 : lambda x, **params : exp(x, x_0, **params)
model_both = lambda x_0 : lambda x, **params : both(x, x_0, **params)

model_poly_err = lambda x_0, params, errs : lambda x, x_err : poly_err(x, x_err, x_0, params, errs)
model_exp_err = lambda x_0, params, errs : lambda x, x_err : exp_err(x, x_err, x_0, params, errs)
model_both_err = lambda x_0, params, errs : lambda x, x_err : both_err(x, x_err, x_0, params, errs)




def poly(x, x_0, **params):
    coeffs = sorted(params.keys())
    p = 0
    for i in range(len(coeffs)):
        if coeffs[i] in poly_coeffs:
            p += params[coeffs[i]] * (x - x_0)**i
    #print (p[int(len(p)/2)])
    return p

def poly_err(x, x_err, x_0, params, errs):
    coeffs = sorted(params.keys())
    stat2 = 0
    sys2 = 0
    print ("========")
    for i, c in enumerate(coeffs):
        if c in poly_coeffs:
            # error in c_i x**i is ((cerr/c)^2 + (i*xerr/x)^2)^(1/2) *c_i x**i
            # must be added in quadrature with others
            #e2 += (params[c]*(x-x_0)**i)**2 * ((errs[c]/c)**2 + (i*x_err/(x-x_0))**2)
            sys2 += (params[c]*(x-x_0)**i)**2 * ((errs[c]/params[c])**2)
            stat2 += (params[c]*(x-x_0)**i)**2 * (i*x_err/(x-x_0))**2
            print ("--")
            print ("x-x0: "+str(x-x_0))
            print ("coeff: "+c)
            print ("coeff value: "+str(params[c]))
            print ("coeff error: "+str(errs[c]))
            print ("contribution: " + str((params[c]*(x-x_0)**i)))
            print ("stat: "+str(np.sqrt((params[c]*(x-x_0)**i)**2 * (i*x_err/(x-x_0))**2)))
            print ("sys: "+str(np.sqrt((params[c]*(x-x_0)**i)**2 * ((errs[c]/params[c])**2))))
            print ("--")
    print ("========")
    return np.sqrt(stat2), np.sqrt(sys2)

def exp(x, x_0, **params):
    #print (params)
    e = params['n'] * np.exp(params['r']*(x-x_0))
    #print (e[int(len(e)/2)])
    return e

def exp_err(x, x_err, x_0, params, errs):
    # + denotes in quadrature here
    # error in n*exp(r*x) is (n_err/n + exp(r*x)*(r_err/r + x_err/x)) * n*exp(r*x)
    n, r = params['n'], params['r']
    n_err, r_err = errs['n'], errs['r']
    #return n * np.exp(r*(x-x_0)) * np.sqrt((n_err/n)**2 + (n*np.exp(r*(x-x_0)))**2 * ((r_err/r)**2 + (x_err/x)**2))
    sys = n * np.exp(r*(x-x_0)) * np.sqrt((n_err/n)**2 + (n*np.exp(r*(x-x_0)))**2 * ((r_err/r)**2))
    stat = n * np.exp(r*(x-x_0)) * np.sqrt((n*np.exp(r*(x-x_0)))**2 * ((x_err/x)**2))
    return stat, sys

both = lambda x, x_0, params : poly(x, x_0, **params) + exp(x, x_0, **params)

def both_err(x, x_err, x_0, params, errs): 
    p = poly_err(x, x_err, x_0, params, errs)
    e = exp_err(x, x_err, x_0, params, errs)
    return np.sqrt(p[0]**2 + e[0]**2), np.sqrt(p[1]**2 + e[1]**2)

"""
def poly_err_split(x, x_err, x_0, params, errs):
    coeffs = sorted(params.keys())
    stat2 = 0
    sys2 = 0
    for i, c in enumerate(coeffs):
        if c in poly_coeffs and c != 'a':
            # error in c_i x**i is ((cerr/c)^2 + (i*xerr/x)^2)^(1/2) *c_i x**i
            # must be added in quadrature with others
            #e2 += (params[c]*(x-x_0)**i)**2 * ((errs[c]/c)**2 + (i*x_err/(x-x_0))**2)
            sys2 += (params[c]*(x-x_0)**i)**2 * ((errs[c]/c)**2)**2
            stat2 += (params[c]*(x-x_0)**i)**2 * (i*x_err/(x-x_0))**2
    return np.sqrt(stat2), np.sqrt(sys2)    

def both_err_split(x, x_err, x_0, params, errs): 
    p = poly_err_split(x, x_err, x_0, params, errs)
    e = exp_err(x, x_err, x_0, params, errs)
    return np.sqrt(p[0]**2 + e[0]**2), np.sqrt(p[1]**2 + e[1]**2)
"""

def poly_params(params, n):
    if n > 6:
        print ("polynomial capped at 6")
        n = 6
    for i in range(n + 1):
        if poly_coeffs[i] != 'b':
            params.add(poly_coeffs[i], value = 0)
        else:
            params.add(poly_coeffs[i], value = 1)

def exp_params(params, nval = 0.01, rval = 1/1000):
    params.add('n', value = nval)
    params.add('r', value = rval, min = -1/500, max = 1/500)
    
def extract(params):
    vals = {}
    errs = {}
    for param in sorted(params.keys()):
        vals[param] = params[param].value
        errs[param] = params[param].stderr
    return vals, errs

