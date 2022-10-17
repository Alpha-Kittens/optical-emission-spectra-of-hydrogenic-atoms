import numpy as np
import lmfit
from scipy.special import wofz
from max_model import regions

def linear(x, m, b):
    """
    linear model y=mx+b

    Arguments: 
        * `x` : parameter value to evaluate the function at
        * `m` : slope of the line
        * `b` : y-intercept of the line

    Returns: 
        * `y` result of y=mx+b
    """
    return (m*x)+b

def linear_extract(result):

    pass

def quadratic(x, a, b, c):
    """
    quadratic model y=ax^2 + bx + c

    Arguments: 
        * `x` : parameter value to evaluate the function at
        * `a` : leading coefficient
        * `b` : 2nd term coefficient
        * `c` : y-interecept

    Returns: 
        * `y` : result of y=ax^2 + bx + c
    """
    return (a*(x**2)) + (b*x) + c

def quadratic_extract(result):

    pass

def quadratic_err(x, x_err, a, a_err, b, b_err, c, c_err):
    """
    Error for result of quadratic model

    Arguments:
        * `x` : parameter value to evaluate function
        * `x_err` : error associated with `x`
        * `a` : leading coefficient
        * `a_err` : error associated with `a`
        * `b` : linear coefficient
        * `b_err` : error associated with `b`
        * `c` : y-intercept
        * `c_err` : error associated with `c`
    """

    x_rel = x_err / x
    a_rel = a_err / a
    b_rel = b_err / b
    c_rel = c_err / c

    print ("a contribution: "+str(((a * x**2)**2 * (a_rel**2 + 2 * x_rel**2))**(1/2)))
    print ("b contribution: "+str(((b * x)**2 * (b_rel ** 2 + x_rel**2))**(1/2)))
    print ("c contribution: "+str(c_err))

    return ((a * x**2)**2 * (a_rel**2 + 2 * x_rel**2) + (b * x)**2 * (b_rel ** 2 + x_rel**2) + c_err ** 2) ** (1/2)

def inverse_quadratic(true, a, b, c):

    return (-b + np.sqrt(b ** 2 - 4 * a * (c - true))) / (2 * a)

def voigt(x, amp, mu, alpha, gamma):
    """
    Voigt, i.e. gaussian convolved with lorentzian

    Argumetns:
        * `x` : parameter value to evaluate the function at
        * `amp` : amplitude of distribution
        * `mu` : center of distribution
        * `alpha` : Half-width half-maximum of Gaussian (monotonic with sigma)
        * 'gamma` : Half-width half-maximum of Lorentzian

    Returns:
        * `y` : evaluation of voigt model at x at given parameters
    """

    sigma = alpha / np.sqrt(2 * np.log(2))

    return amp * np.real(wofz((x - mu + 1j*gamma)/sigma/np.sqrt(2))) / sigma / np.sqrt(2*np.pi)

def voigt_extract(result):

    return ((result.params['amp'].value, result.params['mu'].value, result.params['alpha'].value, result.params['gamma'].value), (result.params['amp'].stderr, result.params['mu'].stderr, result.params['alpha'].stderr, result.params['gamma'].stderr))

voigt_with_shift = lambda x, amp, mu, alpha, gamma, a : voigt(x, amp, mu, alpha, gamma) + a

def voigt_with_shift_extract(result):

    return ((result.params['amp'].value, result.params['mu'].value, result.params['alpha'].value, result.params['gamma'].value, result.params['a'].value), (result.params['amp'].stderr, result.params['mu'].stderr, result.params['alpha'].stderr, result.params['gamma'].stderr, result.params['a'].stderr))

two_voigt = lambda x, amp, mu, alpha, gamma, amp2, mu2, alpha2, gamma2, a : voigt(x, amp, mu, alpha, gamma) + voigt(x, amp2, mu2, alpha2, gamma2) + a

def two_voigt_extract(result):

    return ((result.params['amp1'].value, result.params['mu1'].value, result.params['alpha1'].value, result.params['gamma1'].value, result.params['a'].value, result.params['amp2'].value, result.params['mu2'].value, result.params['alpha2'].value, result.params['gamma2'].value), (result.params['amp1'].stderr, result.params['mu1'].stderr, result.params['alpha1'].stderr, result.params['gamma1'].stderr, result.params['a'].stderr, result.params['amp2'].stderr, result.params['mu2'].stderr, result.params['alpha2'].stderr, result.params['gamma2'].stderr))

def get_voigt_params(data, return_max = False):

    wavelengths, cps = data[:,0], data[:,1]
    
    signals = regions(cps)[1]

    maxes = [max(cps[signal[0]:signal[1]+1]) for signal in signals]

    max_signal = signals[np.argmax(maxes)]
    arg_max = np.argmax(cps)
    start_mu, start_amp = wavelengths[arg_max], cps[arg_max]
    print (start_mu)
    print (start_amp)
    start_alpha = (wavelengths[max_signal[1]] - wavelengths[max_signal[0]]) / 4
    start_gamma = start_alpha / 4
    print (start_alpha)
    print (start_gamma)
    start_amp *= start_amp / max(voigt(data[:,0], start_amp, start_mu, start_alpha, start_gamma))

    if return_max:
        return start_amp, start_mu, start_alpha, start_gamma, arg_max
    return start_amp, start_mu, start_alpha, start_gamma

def voigt_params(data, return_max = False):

    params = lmfit.Parameters()

    start_amp, start_mu, start_alpha, start_gamma = get_voigt_params(data, return_max = return_max)

    import matplotlib.pyplot as plt

    #plt.scatter(data[:,0], data[:,1])
    #plt.plot(data[:,0], voigt(data[:,0], start_amp, start_mu, start_alpha, start_gamma) + min(data[:,1]))
    #plt.show()

    params.add('amp', value = start_amp, min = 0)
    params.add('mu', value = start_mu)
    params.add('alpha', value = start_alpha, min = 0)
    params.add('gamma', value = start_gamma, min = 0)

    return params


def voigt_shift_params(data):

    params = voigt_params(data)
    params.add('a', min(data[:,1]))

    return params

def two_voigt_params(data, expected_shift, stepsize = 0.0025):

    params = lmfit.Parameters()

    start_amp, start_mu, start_alpha, start_gamma, arg_max = get_voigt_params(data, return_max = True)

    if arg_max + int(expected_shift // stepsize) > len(data[:,1]):
        start_amp2 = data[:,1][max(arg_max - int(expected_shift // stepsize), 0)]
        start_mu2 = data[:,0][max(arg_max - int(expected_shift // stepsize), 0)]
    elif arg_max - int(expected_shift // stepsize) < 0:
        start_amp2 = data[:,1][arg_max + int(expected_shift // stepsize)]
        start_mu2 = data[:,0][arg_max + int(expected_shift // stepsize)]        
    elif data[:,1][arg_max + int(expected_shift // stepsize)] > data[:,1][arg_max - int(expected_shift // stepsize)]:
        start_amp2 = data[:,1][arg_max + int(expected_shift // stepsize)]
        start_mu2 = data[:,0][arg_max + int(expected_shift // stepsize)]
    else:
        start_amp2 = data[:,1][arg_max - int(expected_shift // stepsize)]
        start_mu2 = data[:,0][arg_max - int(expected_shift // stepsize)]
    start_amp2 *= start_amp2 / max(voigt(data[:,0], start_amp2, start_mu2, start_alpha, start_gamma))

    params.add('amp', value = start_amp, min = 0)
    params.add('mu', value = start_mu)
    params.add('alpha', value = start_alpha, min = 1e-6)
    params.add('gamma', value = start_gamma, min = 0)
    params.add('amp2', value = start_amp2, min = 0)
    params.add('mu2', value = start_mu2, min = 1e-6)
    params.add('alpha2', value = start_alpha, min = 0)
    params.add('gamma2', value = start_gamma)
    params.add('a', min(data[:,1]), min = 0)

    return params

voigt_models = {
    'voigt' : (voigt, voigt_params, voigt_extract),
    'voigt_with_shift' : (voigt_with_shift, voigt_shift_params, voigt_with_shift_extract),
    'two_voigt' : (two_voigt, two_voigt_params, two_voigt_extract),
}

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    x = np.linspace(5748.5, 5750.5, 100)
    plt.plot(x, voigt(x, 1144, 5749, 0.05, 0.05))
    print(max(voigt(x, 1144, 5749, 0.05, 0.05)))
    plt.show()