import import_ipynb
from lorentzian_fit import fit_to_lorentzian

from data.data_loader import read_data
from lmfit import Model, Parameter

reference = {
    'mercury' : {
        'reference' : [3650.15, 4046.56, 4368.33, 5460.74, 5769.60, 5790.66],
        'check' : [], # TODO
    },
    'hydrogen' : {
        'reference' : [],
    },
    'sodium' : {
        'reference' : [],
    },


}

def extract_fit_data(fit):

    return fit.params['amp'], fit.params['cen'], fit.params['scale']


def get_band(file):

    data = read_data
    fit = fit_to_lorentzian(data)
    return extract_fit_data(fit)


def get_calibration_data(files):
    """
    Given list of files (consdier: element name, and it just knows the files), extracts the band details from each of them. 
    Returs centers and weights
    """

    centers = []
    widths = []
    for file in files:

        band_details = get_band(file)
        centers.append(band_details[0])
        widths.append(band_details[1])

    return centers, widths

def fit_calibration(files, element):

    centers, widths = get_calibration_data(files)

    #sort centers

    residuals = [abs(centers[i] - reference[element]['reference'][i]) for i in range(len(centers))]

    #what do we expect? Try quadratic.

    quadratic = lambda x, a, b, c : a * x ** 2 + b * x + c

    model = Model(quadratic)

    result = model.fit(centers, reference[element]['reference']) # weight by widths

    return quadratic, (result.params['a'], result.params['b'], result.params['c'])

    raise NotImplementedError

def check_calibration(model, files, element):
    """
    Given a list of smaller-line calibration references, checks validity of fit. 
    """

    raise NotImplementedError

